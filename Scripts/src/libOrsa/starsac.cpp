// @INPROCEEDINGS{5206678,
//  author={Jongmoo Choi and Medioni, Gerard},
//  booktitle={2009 IEEE Conference on Computer Vision and Pattern Recognition},
//  title={StaRSaC: Stable random sample consensus for parameter estimation},
//  year={2009},
//  volume={},
//  number={},
//  pages={675-682},
//  doi={10.1109/CVPR.2009.5206678}}
//
// Created by Clément Riu.
//

#include <limits>

#include "starsac.hpp"
#include "libNumerics/numerics.h"
#include "sampling.hpp"

namespace orsa {

StarSac::StarSac(const ModelEstimator *estimator) : LRTSac(estimator) {
}

/// Setters for RANSAC hyperparameters.
/// \param precision: max error to consider a point inlier
/// \param cpII: stopping criterion confidence
void StarSac::setHyperParameters(const double cpIIT, const double sigmaMax) {
    setCpIIT(cpIIT);
    setSigmaMax(sigmaMax);
}

double StarSac::computeVar(const std::vector<RunResult> &estimatedResults) const {
    const int modelDimension = (estimatedResults[0].model).ncol() * (estimatedResults[0].model).nrow();
    libNumerics::matrix<double> estimatedResultsAsMatrix(modelDimension, _numRansacRun);
    for (size_t iModel = 0; iModel < estimatedResults.size(); iModel++) {
        for (int iParam = 0; iParam < modelDimension; iParam++) {
            estimatedResultsAsMatrix(iParam, iModel) = (estimatedResults[iModel].model)(iParam);
        }
    }
    libNumerics::matrix<double> covMatrix = estimatedResultsAsMatrix * estimatedResultsAsMatrix.t();
    libNumerics::SVD svd(covMatrix);
    return svd.D(0);
}

double StarSac::run(RunResult &res, int nIterMax, bool verbose) const {
    const int nData = _model->NbData();
    const int sizeSample = _model->SizeSample();
    res.vInliers.clear();
    res.vpm = nData;
    res.T = 0;
    res.sigma = _sigmaMax;
    double bestVar = std::numeric_limits<double>::infinity();

    if (nData <= sizeSample)
        return bestVar;

    if (nIterMax <= 0) // Authorize "infinite" number of iterations
        nIterMax = std::numeric_limits<int>::max();

    // Computation of array of values for sigma
    std::vector<double> Sigma;
    initSigma(Sigma);

    RunResult bestResultOverall;
    std::vector<std::vector<RunResult>> results(Sigma.size());

    std::vector<RunResult> bestRestultLocal(Sigma.size());
    std::vector<RunResult>::iterator itBestRestultLocal = bestRestultLocal.begin();
    for (; itBestRestultLocal != bestRestultLocal.end(); itBestRestultLocal++) {
        (*itBestRestultLocal).vInliers.clear();
    }
    for (int nRun = 0; nRun < _numRansacRun; nRun++) {
        if (verbose) {
            std::cout << "Run " << nRun << " of " << _numRansacRun << std::endl;
        }
        std::vector<RunResult> tempRes;
        res.T += ransac(Sigma, tempRes, nIterMax);
        for (size_t iSigma = 0; iSigma < Sigma.size(); iSigma++) {
            results[iSigma].push_back(tempRes[iSigma]);
            if (tempRes[iSigma].vInliers.size() > bestRestultLocal[iSigma].vInliers.size()) {
                bestRestultLocal[iSigma] = tempRes[iSigma];
            }
        }
    }
    for (size_t iSigma = 0; iSigma < Sigma.size(); iSigma++){
        double var = computeVar(results[iSigma]);
        if (var < bestVar) {
            bestVar = var;
            res.sigma = Sigma[iSigma];
            bestResultOverall = bestRestultLocal[iSigma];
        }
    }

    if (verbose) {
        std::cout << "Optimal threshold: " << res.sigma << std::endl;
    }
    res.vInliers = bestResultOverall.vInliers;
    res.model = bestResultOverall.model;
    return bestVar;
}

/// Find inliers of \a model passed as parameter (within \a precision).
void StarSac::FindAllInliers(const ModelEstimator::Model &model,
                             const std::vector<double> &precisions,
                             std::vector<std::vector<int>> &inliers) const {
    const int nData = _model->NbData();
    for (int i = 0; i < nData; i++) {
        double error = _model->Error(model, i);
        std::vector<double>::const_iterator itPrec = precisions.begin();
        for (int j = 0; itPrec != precisions.end(); itPrec++, j++) {
            if (error <= (*itPrec) * (*itPrec)) inliers[j].push_back(i);
        }
    }
}

int StarSac::ransac(const std::vector<double> &precisions, std::vector<RunResult> &res, int nIterMax) const {
    double log_pII = log(1 - _cpIIT);

    const int nData = _model->NbData();
    const int sizeSample = _model->SizeSample();
    res.clear();
    for (size_t i = 0; i < precisions.size(); i++) {
        RunResult result;
        result.vInliers.clear();
        result.sigma = precisions[i];
        result.vpm = nData;
        res.push_back(result);
    }
    std::vector<int> nIterMaxes(precisions.size(), nIterMax);
    int T = 0;
    for (; T < nIterMax; T++) {
        std::vector<int> vSample(sizeSample); // Sample indices
        UniformSample(sizeSample, nData, &vSample); // Get random sample
        std::vector<ModelEstimator::Model> vModels;
        _model->Fit(vSample, &vModels);
        std::vector<ModelEstimator::Model>::const_iterator it;
        for (it = vModels.begin(); it != vModels.end(); ++it) {
            std::vector<std::vector<int>> inliers(precisions.size());
            FindAllInliers(*it, precisions, inliers);
            std::vector<double>::const_iterator itPrec = precisions.begin();
            bool change = false;
            for (int j = 0; itPrec != precisions.end(); itPrec++, j++) {
                if (nIterMaxes[j] > T) {
                    if (res[j].vInliers.size() < inliers[j].size()) {
                        res[j].model = *it;

                        std::swap(inliers[j], res[j].vInliers); // Avoid copy
                        double pIn = pow(res[j].vInliers.size() / (double) nData, sizeSample);
                        double denom = log(1 - pIn);

                        if (denom < 0) { // Protect against 1-eps==1
                            change = true;
                            double newIter = log_pII / denom;
                            nIterMaxes[j] = (int) std::min((double) nIterMaxes[j], ceil(newIter));
                        }
                    }
                }
            }
            if (change) {
                int maxIterMax = 0;
                for (size_t j = 0; j < precisions.size(); j++) {
                    if (maxIterMax < nIterMaxes[j]) maxIterMax = nIterMaxes[j];
                }
                nIterMax = maxIterMax;
            }
        }
    }
    return T;
}

/// Verify the runOutput metric is good enough.
bool StarSac::satisfyingRun(double runOutput) const {
    return true;
}

} // namespace orsa
