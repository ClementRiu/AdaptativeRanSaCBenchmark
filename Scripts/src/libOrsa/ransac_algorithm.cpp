/**
* @file ransac_algorithm.cpp
* @brief Generic class for Ransac algorithms
* @author Clement Riu
*
* Copyright (c) 2021 Clement Riu
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

#include "ransac_algorithm.hpp"
#include <ctime>

namespace orsa {
/// Refine model based on all inliers, and display statistics.
static void refine(const ModelEstimator *model,
                   const std::vector<int> &vec_inliers,
                   ModelEstimator::Model *M,
                   bool verbose = true) {
    std::pair<double, double> err = model->ErrorStats(vec_inliers, *M);
    if (verbose)
        std::cout << "Before refinement: RMSE/max error: " << err.first
                  << "/" << err.second << std::endl;
    ModelEstimator::Model M2(*M);
    if (model->ComputeModel(vec_inliers, &M2)) { // Re-estimate with all inliers
        double maxBefore = err.second;
        std::pair<double, double> err2 = model->ErrorStats(vec_inliers, M2);
        if (verbose)
            std::cout << "After  refinement: RMSE/max error: " << err2.first
                      << "/" << err2.second << std::endl;
        if (err2.first <= maxBefore)
            *M = M2;
        else
            std::cerr<<"Warning: error after refinement too large and ignored"
                      << std::endl;
    } else
        std::cerr<<"Warning: error in refinement, result is suspect"<<std::endl;
}

/// Estimate a given model for a given algorithm, with refinement.
/// \param[in] model: Model estimator.
/// \param[in] algorithm: Algorithm to use.
/// \param[in] maxIter: Maximum number of iteration of the algorithm.
/// \param[out] modelParams: Estimated model params.
/// \param[out] vec_inliers: Index of inliers in model.data_
/// \param[out] T: Number of iterations.
/// \param[out] vpm: Number of verification per model.
/// \param[out] runtime: Runtime of the estimation.
/// \param[out] computedSigma: Estimated sigma for adaptative methods.
/// \param[in] verbose: Set the run information verbose.
bool RansacAlgorithm::evalModel(RunResult& res, double& runtime,
                                int maxIter, bool verbose) const {
    if (_model->NbData() >= _model->SizeSample()) {
        clock_t begin = std::clock();
        double runOutput = run(res, maxIter, verbose);
        clock_t end = std::clock();
        runtime = double(end - begin) / CLOCKS_PER_SEC;

        if (! satisfyingRun(runOutput))
            return false;
        refine(_model, res.vInliers, &res.model, verbose);
        return true;
    }
    std::cerr << "Error: The algorithm needs " << _model->SizeSample()
              << " matches or more to proceed" << std::endl;
    return false;
}
} // namespace orsa
