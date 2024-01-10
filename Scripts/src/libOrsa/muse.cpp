// Adapted from https://github.com/vxl/vxl/tree/c6c899aaf9cbf7ccabfca1ba35b2248fb611ffbc
// Algorithm from MUSE: Robust Surface Fitting using Unbiased Scale Estimates, James V. Miller and Charles V. Stewart
// By Clément Riu
// 2021

#include "muse.hpp"
#include "sampling.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

namespace orsa {

/// Constructor
    Muse::Muse(const ModelEstimator *estimator)
            : RansacAlgorithm(estimator) {
        setHyperParameters();
        _table = new MuseTable(_model->NbData());
    }

/// Setters for RANSAC hyperparameters.
/// \param precision: max error to consider a point inlier
/// \param cpII: stopping criterion confidence
/// \param nModelMin: Number of non-contaminated model to see before stopping
    void Muse::setHyperParameters(const double cpIIT,
                                  const double min_frac,
                                  const double max_frac,
                                  const double frac_inc,
                                  const bool use_sk_refine) {
        setCpIIT(cpIIT);
        setMinFrac(min_frac);
        setMaxFrac(max_frac);
        setFracInc(frac_inc);
        setSKRefine(use_sk_refine);
    }

    void Muse::setCpIIT(const double cpIIT) {
        _cpIIT = cpIIT;
    }

    void Muse::setMinFrac(const double min_frac) {
        _min_frac = min_frac;
    }

    void Muse::setMaxFrac(const double max_frac) {
        _max_frac = max_frac;
    }

    void Muse::setFracInc(const double frac_inc) {
        _frac_inc = frac_inc;
    }

    void Muse::setSKRefine(const bool use_sk_refine) {
        _use_sk_refine = use_sk_refine;
    }

/// Generic implementation of RANSAC
    double Muse::run(RunResult &res, int nIterMax, bool verbose) const {
        double log_pII = log(1 - _cpIIT);

        const int nData = _model->NbData();
        const int sizeSample = _model->SizeSample();

        if (nData <= sizeSample)
            return 0;

        res.vInliers.clear();
        res.vpm = nData;

        double minObjective = std::numeric_limits<double>::infinity();
        bool modelFound = false;

        for (res.T = 0; res.T < nIterMax; res.T++) {
            std::vector<int> vSample(sizeSample); // Sample indices
            UniformSample(sizeSample, nData, &vSample); // Get random sample

            std::vector<ModelEstimator::Model> vModels;
            _model->Fit(vSample, &vModels);

            std::vector<ModelEstimator::Model>::const_iterator it;
            for (it = vModels.begin(); it != vModels.end(); ++it) {
                std::vector<double> residuals;
                _model->computeResiduals(*it, residuals);
                double objective = computeObjective(residuals);

                if (objective < minObjective) {
                    modelFound = true;
                    res.model = *it;
                    minObjective = objective;

                    std::vector<int> vInliers;
                    _model->FindInliers(res.model, objective, vInliers);
                    double pIn = pow(vInliers.size() / (double) nData, sizeSample);
                    double denom = log(1 - pIn);

                    if (denom < 0) { // Protect against 1-eps==1
                        double newIter = log_pII / denom;
                        nIterMax = (int) std::min((double) nIterMax, ceil(newIter));
                    }
                }
            }
        }
        if (modelFound) {
            std::vector<double> residuals;
            _model->computeResiduals(res.model, residuals);
            res.sigma = computeObjective(residuals);
            _model->FindInliers(res.model, res.sigma, res.vInliers);
            return res.sigma;
        }
        return 0;
    }

    double Muse::computeObjective(std::vector<double> &residuals) const {
        double sigma_est;
        int best_k = 0;

        std::sort(residuals.begin(), residuals.end());

        unsigned int num_residuals = residuals.size();
        bool at_start = true;
        double best_sk = 0;
        double best_objective = 0;

        constexpr double min_exp_kth_to_stddev_ratio = 3.0;
        static bool notwarned = true;

        double sum_residuals = 0;
        double best_sum = 0;
        int prev_k = 0;

        //  Find the best k
        for (double frac = _min_frac; frac <= _max_frac + 0.00001; frac += _frac_inc) {
            unsigned int k = std::round(frac * num_residuals);
            if (k > num_residuals) k = num_residuals;
            if (k <= 0) k = 1;
            if (_table->expected_kth(k, num_residuals) / _table->standard_dev_kth(k, num_residuals)
                < min_exp_kth_to_stddev_ratio) {
                if (notwarned) {
                    std::cerr << "WARNING: attempted evaluation at value of k that lead to unstable estimates\n";
                    notwarned = false;
                }
                continue;
            }

            for (unsigned int i = prev_k; i < k; ++i) {
                sum_residuals += std::sqrt(residuals[i]);
            }
            double sk = sum_residuals / _table->muset_divisor(k, num_residuals);
            double objective = sk * _table->standard_dev_kth(k, num_residuals) /
                               _table->expected_kth(k, num_residuals);

            if (at_start || objective < best_objective) {
                at_start = false;
                best_k = k;
                best_sk = sk;
                best_sum = sum_residuals;
                best_objective = objective;
            }
            prev_k = k;
        }
        if (at_start) {
            std::cerr << "WARNING:  There were NO values of k with stable estimates.\n"
                      << "          Setting sigma = +Infinity\n";
            return std::numeric_limits<double>::infinity();
        }

        //  Refine the scale estimate
        if (!_use_sk_refine) {
            sigma_est = best_sk;
        } else {
            unsigned int new_n = best_k;
            while (new_n < num_residuals && std::sqrt(residuals[new_n]) < 2.5 * best_sk)
                ++new_n;
            sigma_est = best_sum / _table->muset_divisor(best_k, new_n);
        }
        return sigma_est;
    }

    // Verifies the runOutput metric is good enough.
    bool Muse::satisfyingRun(const double runOutput) const {
        return runOutput > 0.0;
    }

}  // namespace orsa
