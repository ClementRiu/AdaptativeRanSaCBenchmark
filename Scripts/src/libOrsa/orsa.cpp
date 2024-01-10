/**
* @file orsa.cpp
* @brief Model estimation by ORSA (aka AC-RANSAC) algorithm.
* @author Lionel Moisan, Pascal Monasse, Pierre Moulon
*
* Copyright (c) 2007 Lionel Moisan
* Copyright (c) 2010-2011,2020 Pascal Monasse
* Copyright (c) 2010-2011 Pierre Moulon
* All rights reserved.
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "orsa.hpp"
#include "sampling.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace orsa {

/// The class does not take ownership of the estimator instance but depends on
/// it. Be careful that it is still valid during the lifetime of Orsa object.
Orsa::Orsa(const ModelEstimator *estimator)
: RansacAlgorithm(estimator) {
    setHyperParameters();
}

//// Setters for RANSAC hyperparameters.
/// \param precision : Parameter controlling the maximum error to consider a possible threshold.
/// \param alpha0Right : Parameter controlling TODO
/// \param alpha0Left : Parameter controlling TODO
/// \param bConvergence : Parameter controlling TODO
void Orsa::setHyperParameters(const double precision,
                              const bool bConvergence) {
    setPrecision(precision);
    setAlpha0Right();
    setAlpha0Left();
    setBConvergence(bConvergence);
}

//// Sets the value of precision.
/// \param precision : Parameter controlling the maximum error to consider a point inlier.
void Orsa::setPrecision(const double precision) {
    assert(precision >= 0);
    _precision = precision;
}

//// Sets the value of alpha0Right.
/// \param alpha0Right : Parameter controlling TODO.
void Orsa::setAlpha0Right() {
    _alpha0Right = _model->pSigma(1, false);
    logalpha0_[1] = log10(_alpha0Right);
}

//// Sets the value of alpha0Left.
/// \param alpha0Left : Parameter controlling TODO.
void Orsa::setAlpha0Left() {
    _alpha0Left = _model->pSigma(1, true);
    logalpha0_[0] = log10(_alpha0Left);
}

//// Sets the value of bConvergence.
/// \param bConvergence : Parameter controlling TODO.
void Orsa::setBConvergence(const bool bConvergence) {
    _bConvergence = bConvergence;
}

/// logarithm (base 10) of binomial coefficient
static float logcombi(int k, int n) {
    if (k >= n || k <= 0) return (0.0);
    if (n - k < k) k = n - k;
    double r = 0.0;
    for (int i = 1; i <= k; i++)
        r += log10((double) (n - i + 1)) - log10((double) i);

    return static_cast<float>(r);
}

/// tabulate logcombi(.,n)
void makelogcombi_n(int n, std::vector<float> &l) {
    l.resize(n + 1);
    for (int k = 0; k <= n; k++)
        l[k] = logcombi(k, n);
}

/// tabulate logcombi(k,.)
void makelogcombi_k(int k, int nmax, std::vector<float> &l) {
    l.resize(nmax + 1);
    for (int n = 0; n <= nmax; n++)
        l[n] = logcombi(k, n);
}

/// Find best NFA and number of inliers wrt square error threshold in e.
Orsa::ErrorIndex Orsa::bestNFA(const std::vector<ErrorIndex> &e,
                               double loge0,
                               double maxThreshold,
                               const std::vector<float> &logc_n,
                               const std::vector<float> &logc_k) const {
    const int startIndex = _model->SizeSample();
    const double multError = (_model->DistToPoint() ? 1.0 : 0.5);

    ErrorIndex bestIndex(std::numeric_limits<double>::infinity(),
                         startIndex,
                         0);

    const int n = static_cast<int>(e.size());
    for (int k=startIndex+1; k <= n && e[k - 1].error <= maxThreshold; ++k) {
        double logalpha = logalpha0_[e[k - 1].side]
            + multError * log10(e[k - 1].error + std::numeric_limits<double>::epsilon());

        ErrorIndex index(loge0 + logalpha * (double) (k - startIndex) + logc_n[k] + logc_k[k],
                         k, e[k - 1].side);
        if (index.error < bestIndex.error)
            bestIndex = index;
    }
    return bestIndex;
}

/// Generic implementation of 'ORSA':
/// A Probabilistic Criterion to Detect Rigid Point Matches
///    Between Two Images and Estimate the Fundamental Matrix.
/// Bibtex :
/// @article{DBLP:journals/ijcv/MoisanS04,
///  author    = {Lionel Moisan and B{\'e}renger Stival},
///  title     = {A Probabilistic Criterion to Detect Rigid Point Matches
///    Between Two Images and Estimate the Fundamental Matrix},
///  journal   = {International Journal of Computer Vision},
///  volume    = {57},
///  number    = {3},
///  year      = {2004},
///  pages     = {201-218},
///  ee        = {http://dx.doi.org/10.1023/B:VISI.0000013094.38752.54},
///  bibsource = {DBLP, http://dblp.uni-trier.de}
///}
///
/// ORSA is based on an a contrario criterion of inlier/outlier discrimination,
/// is parameter free and relies on an optimized
/// random sampling procedure. It returns the log of NFA and
/// the best estimated model.
/// \param res Output results
/// \param verbose Display optimization statistics
double Orsa::run(RunResult& res, int nIterMax, bool verbose) const {
    const int nData = _model->NbData();
    const int sizeSample = _model->SizeSample();
    if (nData <= sizeSample)
        return std::numeric_limits<double>::infinity();

    const double maxThreshold = (_precision > 0) ?
        _precision * _precision : // Square max error
        std::numeric_limits<double>::infinity();

    std::vector<ErrorIndex> vResiduals(nData); // [residual,index]
    std::vector<int> vSample(sizeSample); // Sample indices

    // Possible sampling indices (could change in the optimization phase)
    std::vector<int> vIndex(nData);
    for (int i = 0; i < nData; ++i)
        vIndex[i] = i;

    // Precompute log combi
    double loge0 = log10((double) _model->NbModels() * (nData - sizeSample));
    std::vector<float> vLogc_n, vLogc_k;
    makelogcombi_n(nData, vLogc_n);
    makelogcombi_k(sizeSample, nData, vLogc_k);

    // Reserve 10% of iterations for focused sampling
    int nIter = nIterMax;
    int nIterReserve = nIter / 10;
    nIter -= nIterReserve;

    // Output parameters
    double minNFA = std::numeric_limits<double>::infinity();
    double errorMax = 0;
    int side = 0;
    res.vInliers.clear();
    res.vpm = nData;

    // Main estimation loop.
    for (res.T = 0; res.T < nIter && vIndex.size() > sizeSample; res.T++) {
        UniformSample(sizeSample, vIndex, &vSample); // Get random sample

        // Evaluate models
        bool better = false;
        std::vector<ModelEstimator::Model> vModels;
        _model->Fit(vSample, &vModels);
        std::vector<ModelEstimator::Model>::const_iterator it;
        for (it=vModels.begin(); it!=vModels.end(); ++it) {
            // Residuals computation and ordering
            for (int i = 0; i < nData; ++i) {
                int s;
                double error = _model->Error(*it, i, &s);
                vResiduals[i] = ErrorIndex(error, i, s);
            }

            // Most meaningful discrimination inliers/outliers
            std::sort(vResiduals.begin(), vResiduals.end());
            ErrorIndex best = bestNFA(vResiduals, loge0, maxThreshold,
                                      vLogc_n, vLogc_k);

            if (best.error < minNFA) {// A better model was found
                if (best.error < 0) res.model = *it;
                better = true;
                minNFA = best.error; // best stores nfa in error field.
                side = best.side;
                res.vInliers.resize(best.index);
                for (int i = 0; i < best.index; ++i)
                    res.vInliers[i] = vResiduals[i].index;
                errorMax = vResiduals[best.index - 1].error; // Error threshold
                if (verbose) {
                    std::cout << "  nfa=" << minNFA
                              << " inliers=" << res.vInliers.size()
                              << " precision=" << sqrt(errorMax)
                              << " im" << side + 1
                              << " (iter=" << res.T;
                    if (best.error < 0) {
                        std::cout << ",sample=" << vSample.front();
                        std::vector<int>::const_iterator it = vSample.begin();
                        for (++it; it != vSample.end(); ++it)
                            std::cout << ',' << *it;
                    }
                    std::cout << ")" << std::endl;
                }
            }
        }
        // ORSA optimization: draw samples among best set of inliers so far
        if ((better && minNFA < 0) || (res.T + 1 == nIter && nIterReserve > 0)) {
            if (res.vInliers.empty()) { // No model found at all so far
                nIter++; // Continue to look for any model, even not meaningful
                nIterReserve--;
            } else {
                vIndex = res.vInliers;
                if (nIterReserve) {
                    nIter = res.T + 1 + nIterReserve;
                    nIterReserve = 0;
                }
            }
        }
    }

    if (minNFA >= 0)
        res.vInliers.clear();

    if (_bConvergence)
        res.T+=refineUntilConvergence(vLogc_n, vLogc_k, loge0,
                                      maxThreshold, minNFA, &res.model, verbose,
                                      res.vInliers, errorMax, side);

    res.sigma = sqrt(errorMax);
    return minNFA;
}

/// Refine the model on all the inliers with the "a contrario" model
/// The model is refined while the NFA threshold is not stable.
/// Return the number of iterations.
int Orsa::refineUntilConvergence(const std::vector<float> &vLogc_n,
                                 const std::vector<float> &vLogc_k,
                                 double loge0,
                                 double maxThreshold,
                                 double &minNFA,
                                 ModelEstimator::Model *model,
                                 bool bVerbose,
                                 std::vector<int> &vInliers,
                                 double &errorMax,
                                 int &side) const {
    std::cout << "\n\n Orsa::refineUntilConvergence(...)\n" << std::endl;
    const int nData = _model->NbData();
    std::vector<ErrorIndex> vResiduals(nData); // [residual,index]

    bool bContinue = true;
    int iter = 0;
    do {
        std::vector<ModelEstimator::Model> vModels;
        _model->Fit(vInliers, &vModels);

        // Evaluate models
        std::vector<ModelEstimator::Model>::const_iterator it;
        for (it=vModels.begin(); it!=vModels.end(); ++it) {
            // Residuals computation and ordering
            for (int i = 0; i < nData; ++i) {
                double error = _model->Error(*it, i);
                vResiduals[i] = ErrorIndex(error, i);
            }

            // Most meaningful discrimination inliers/outliers
            std::sort(vResiduals.begin(), vResiduals.end());
            ErrorIndex best = bestNFA(vResiduals, loge0, maxThreshold,
                                      vLogc_n, vLogc_k);

            if (best.error < 0 && best.error < minNFA) { // better model found
                minNFA = best.error;
                side = best.side;
                vInliers.resize(best.index);
                for (int i = 0; i < best.index; ++i)
                    vInliers[i] = vResiduals[i].index;
                errorMax = vResiduals[best.index - 1].error; // Error threshold
                *model = *it;

                if (bVerbose)
                    std::cout << "  nfa=" << minNFA
                              << " inliers=" << vInliers.size()
                              << " precision=" << errorMax
                              << " (iter=" << iter << ")\n";

            } else
                bContinue = false;
        }
        if (vModels.empty())
            bContinue = false;
        ++iter;
    } while (bContinue);
    return iter;
}

/// Toggle iterative refinement NFA/RMSE.
void Orsa::setRefineUntilConvergence(bool value) {
    _bConvergence = value;
}

/// Iterative refinement NFA/RMSE.
bool Orsa::getRefineUntilConvergence() const {
    return _bConvergence;
}

// Verifies the runOutput metric is good enough.
bool Orsa::satisfyingRun(const double runOutput) const {
    return runOutput < 0.0;
}
} // namespace orsa
