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

#include "orsa_fast.hpp"
#include "sampling.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include "utilities/histogram.hpp"

namespace orsa {

/// The class does not take ownership of the estimator instance but depends on
/// it. Be careful that it is still valid during the lifetime of OrsaFast object.
OrsaFast::OrsaFast(const ModelEstimator *estimator)
        : Orsa(estimator) {
    setHyperParameters();
}

/// Find best NFA and number of inliers wrt square error threshold in e.
std::pair<double, double> OrsaFast::bestNFA(const std::vector<double> &e,
                                            double loge0,
                                            double maxThreshold,
                                            const std::vector<float> &logc_n,
                                            const std::vector<float> &logc_k) const {
    const double multError = (_model->DistToPoint() ? 1.0 : 0.5);

    const int nBins = 20;
    utility::Histogram<double> histo(0.0f, maxThreshold, nBins);
    histo.Add(e.cbegin(), e.cend());

    std::pair<double, double> current_best_nfa(std::numeric_limits<double>::infinity(), 0.0);
    unsigned int cumulative_count = 0;
    const std::vector<size_t> &frequencies = histo.GetHist();
    const std::vector<double> residual_val = histo.GetXbinsValue();
    for (int bin = 0; bin < nBins; ++bin) {
        cumulative_count += frequencies[bin];
        if (cumulative_count > _model->SizeSample()
            && residual_val[bin] > std::numeric_limits<float>::epsilon()) {
            const double logalpha = logalpha0_[0]
                                    + multError * log10(residual_val[bin]
                                                        + std::numeric_limits<float>::epsilon());
            const std::pair<double, double> current_nfa(loge0
                                                        + logalpha *
                                                          (double) (cumulative_count - _model->SizeSample())
                                                        + logc_n[cumulative_count]
                                                        + logc_k[cumulative_count],
                                                        residual_val[bin]);
            // Keep the best NFA iff it is meaningful ( NFA < 0 ) and better than the existing one
            if (current_nfa.first < current_best_nfa.first)
                current_best_nfa = current_nfa;
        }
    }
    return current_best_nfa;
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
double OrsaFast::run(RunResult &res, int nIterMax, bool verbose) const {
    const int nData = _model->NbData();
    const int sizeSample = _model->SizeSample();
    if (nData <= sizeSample)
        return std::numeric_limits<double>::infinity();

    const double maxThreshold = (_precision > 0) ?
                                _precision * _precision : // Square max error
                                std::numeric_limits<double>::infinity();

    std::vector<double> vResiduals(nData); // [residual,index]
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
        for (it = vModels.begin(); it != vModels.end(); ++it) {
            // Residuals computation and ordering
            for (int i = 0; i < nData; ++i) {
                int s;
                vResiduals[i] = _model->Error(*it, i, &s);
            }

            // Most meaningful discrimination inliers/outliers
            std::pair<double, double> bestnfa = bestNFA(vResiduals, loge0, maxThreshold,
                                                        vLogc_n, vLogc_k);

            if (bestnfa.first < minNFA) {// A better model was found
                res.model = *it;
                better = true;
                minNFA = bestnfa.first;
                side = 0;
                res.vInliers.clear();
                for (int i = 0; i < nData; ++i) {
                    if (vResiduals[i] <= bestnfa.second) {
                        res.vInliers.push_back(i);
                    }
                }
                errorMax = bestnfa.second; // Error threshold
                if (verbose) {
                    std::cout << "  nfa=" << minNFA
                              << " inliers=" << res.vInliers.size()
                              << " precision=" << sqrt(errorMax)
                              << " im" << side + 1
                              << " (iter=" << res.T;
                    if (bestnfa.first < 0) {
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
        res.T += refineUntilConvergence(vLogc_n, vLogc_k, loge0,
                                        maxThreshold, minNFA, &res.model, verbose,
                                        res.vInliers, errorMax, side);

    res.sigma = sqrt(errorMax);
    return minNFA;
}
} // namespace orsa
