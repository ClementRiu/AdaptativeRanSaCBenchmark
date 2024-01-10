/**
* @file orsa.hpp
* @brief Model estimation by ORSA (aka AC-RANSAC) algorithm.
* @author Pascal Monasse, Pierre Moulon
* 
* Copyright (c) 2011,2020 Pascal Monasse
* Copyright (c) 2011 Pierre Moulon
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

#ifndef ORSA_H
#define ORSA_H

#include "ransac_algorithm.hpp"

namespace orsa {

/// Model estimation with ORSA algorithm.
class Orsa : public RansacAlgorithm {
public:
    /// Constructor
    explicit Orsa(const ModelEstimator *estimator);

    /// Generic implementation of ORSA (Optimized Random Sampling Algorithm)
    double run(RunResult& res, int nIterMax=1000, bool verbose=false) const;

    /// Enable the model refinement until convergence
    void setRefineUntilConvergence(bool value);

    /// Return if convergence check is on or off
    bool getRefineUntilConvergence() const;

    //// Setters for RANSAC hyperparameters
    void setHyperParameters(double precision = 0, bool bConvergence = false);
    void setPrecision(double precision);
    void setAlpha0Right();
    void setAlpha0Left();
    void setBConvergence(bool bConvergence = false);

    // Verifies the runOutput metric is good enough.
    bool satisfyingRun(double runOutput) const;

protected:
    double logalpha0_[2]; ///< Log probability of error<=1, set by subclass
    bool _bConvergence; ///< Refine NFA until it is stable
    double _precision; ///< Maximum inlier/outlier threshold
    /// Probabilities of having an error of 1 pixel in left or right image
    double _alpha0Left, _alpha0Right;

protected:
    /// Distance and associated index
    struct ErrorIndex {
        double error; ///< Square error
        int index; ///< Correspondence index
        int side;     ///< Error in image 1 (side=0) or 2 (side=1)?
        /// Constructor
        ErrorIndex(double e=0, int i=0, int s=0) : error(e),index(i),side(s) {}

        bool operator<(const ErrorIndex &e) const { return (error<e.error); }
    };

    ErrorIndex bestNFA(const std::vector<ErrorIndex> &e,
                       double loge0, double maxThreshold,
                       const std::vector<float> &vLogc_n,
                       const std::vector<float> &vLogc_k) const;

    /// Iterative minimization NFA/RMSE.
    int refineUntilConvergence(const std::vector<float> &vLogc_n,
                               const std::vector<float> &vLogc_k,
                               double loge0,
                               double maxThreshold,
                               double &minNFA,
                               ModelEstimator::Model *model,
                               bool bVerbose,
                               std::vector<int> &vInliers,
                               double &errorMax,
                               int &side) const;
};

void makelogcombi_n(int n, std::vector<float> &l);
void makelogcombi_k(int k, int nmax, std::vector<float> &l);

}  // namespace orsa

#endif
