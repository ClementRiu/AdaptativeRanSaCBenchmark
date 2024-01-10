/**
* @file ransac.cpp
* @brief Model estimation by classical RANSAC algorithm.
* @author Pascal Monasse
*
* Copyright (c) 2020-2021 Pascal Monasse
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

#include "ransac.hpp"
#include "sampling.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

namespace orsa {

/// Constructor
Ransac::Ransac(const ModelEstimator *estimator)
: RansacAlgorithm(estimator) {
    setHyperParameters();
}

/// Setters for RANSAC hyperparameters.
/// \param precision: max error to consider a point inlier
/// \param cpII: stopping criterion confidence
/// \param nModelMin: Number of non-contaminated model to see before stopping
void Ransac::setHyperParameters(double precision, double cpII, int nModelMin) {
    setPrecision(precision);
    setCpII(cpII);
    setNModelMin(nModelMin);
}

/// Sets the value of precision.
/// \param precision: max error to consider a point inlier
void Ransac::setPrecision(double precision) {
    assert(precision > 0);
    _precision = precision;
}

/// Sets the value of cpII, the confidence probability wrt type II error.
/// This is used in the dynamic adjustment of the number of iterations.
/// \param cpII Stopping criterion confidence
void Ransac::setCpII(double cpII) {
    assert(cpII > 0);
    _cpII = cpII;
    if (_cpII > 1.0) {
        std::cerr<<"RANSAC cpII parameter adjusted to not exceed 1"<<std::endl;
        _cpII = 1.0;
    }
}

//// Sets the value of nModelMin.
/// \param nModelMin: Number of non-contaminated model to see before stopping
void Ransac::setNModelMin(int nModelMin) {
    assert(nModelMin > 0);
    _nModelMin = nModelMin;
}

/// Some information output
void display_info(size_t nInliers, size_t iter, size_t nIterMax,
                  const std::vector<int> &vSample) {
    std::cout << " inliers=" << nInliers
              << " (iter=" << iter
              << ",iterMax=" << nIterMax;
    std::cout << ",sample=" << vSample.front();
    std::vector<int>::const_iterator it = vSample.begin();
    for (++it; it != vSample.end(); ++it)
        std::cout << ',' << *it;
    std::cout << ")" << std::endl;
}

/// Generic implementation of RANSAC
double Ransac::run(RunResult& res, int nIterMax, bool verbose) const {
    double log_pII = log(1 - _cpII);

    const int nData = _model->NbData();
    const int sizeSample = _model->SizeSample();
    res.vInliers.clear();
    res.sigma = _precision;
    res.vpm = nData;
    for (res.T = 0; res.T < nIterMax; res.T++) {
        std::vector<int> vSample(sizeSample); // Sample indices
        UniformSample(sizeSample, nData, &vSample); // Get random sample
        std::vector<ModelEstimator::Model> vModels;
        _model->Fit(vSample, &vModels);
        std::vector<ModelEstimator::Model>::const_iterator it;
        for (it=vModels.begin(); it!=vModels.end(); ++it) {
            std::vector<int> inliers;
            _model->FindInliers(*it, _precision, inliers);
            if (res.vInliers.size() < inliers.size()) {
                res.model = *it;

                std::swap(inliers, res.vInliers); // Avoid copy
                double pIn = pow(res.vInliers.size()/(double)nData, sizeSample);
                double denom = log(1 - pIn);

                if (denom < 0) { // Protect against 1-eps==1
                    double newIter = log_pII / denom;
                    if (_nModelMin > 1 && pIn < 1) {
                        double iterUpdate = 0;
                        for (int nModel = 0; nModel < _nModelMin; nModel++)
                            iterUpdate += std::pow(pIn / (1 - pIn),
                                                   (double)nModel);
                        newIter -= log(iterUpdate) / denom;
                    }

                    if (_nModelMin > 1 && pIn == 1)
                        newIter = _nModelMin;
                    nIterMax = (size_t)std::min((double)nIterMax,ceil(newIter));
                }
                if (verbose)
                    display_info(res.vInliers.size(), res.T, nIterMax, vSample);
            }
        }
    }
    return (double) res.vInliers.size();
}

bool Ransac::satisfyingRun(double runOutput) const {
    return runOutput > 0;
}

}  // namespace orsa
