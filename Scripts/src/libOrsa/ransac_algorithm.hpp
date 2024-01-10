/**
* @file ransac_algorithm.hpp
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

#ifndef MMM_ORSA_RANSAC_ALGORITHM_HPP
#define MMM_ORSA_RANSAC_ALGORITHM_HPP

#include "model_estimator.hpp"

namespace orsa {

/// Generic class for Ransac algorithms
class RansacAlgorithm {
public:
    struct RunResult {
        std::vector<int> vInliers; ///< Index of inliers
        ModelEstimator::Model model; ///< Best model parameters
        double sigma; ///< Inlier/outlier threshold (may be estimated)
        int T; ///< Number of performed RANSAC iterations
        double vpm; ///< Mean numbero of verifications per model
    };

    explicit RansacAlgorithm(const ModelEstimator *model): _model(model) {}
    virtual ~RansacAlgorithm() {}

    /// Core Ransac run
    virtual double run(RunResult& res, int maxIter, bool verbose) const=0;
    /// Run + refinement with inliers
    bool evalModel(RunResult& res, double& runtime,
                   int maxIter, bool verbose) const;    
    
    /// Verify the runOutput metric is good enough
    virtual bool satisfyingRun(double runOutput) const { return true; }
protected:
    const ModelEstimator *_model;
};

} // namespace orsa

#endif //MMM_ORSA_RANSAC_ALGORITHM_HPP
