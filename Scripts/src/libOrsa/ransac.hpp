/**
* @file ransac.hpp
* @brief Model estimation by classical RANSAC algorithm.
* @author Pascal Monasse
* 
* Copyright (c) 2020 Pascal Monasse
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

#ifndef RANSAC_H
#define RANSAC_H

#include "ransac_algorithm.hpp"
#include "model_estimator.hpp"

namespace orsa {

/// Some information output
void display_info(size_t nInliers, size_t iter, size_t nIterMax,
                  const std::vector<int> &vSample);

/// Model estimation with ORSA algorithm.
class Ransac : public RansacAlgorithm {
public:
    /// Constructor
    explicit Ransac(const ModelEstimator *estimator);

    //// Setters for RANSAC hyperparameters.
    void setHyperParameters(double precision=1, double cpII=0.99,
                            int nModelMin=1);
    void setPrecision(double precision);
    void setCpII(double cpII);
    void setNModelMin(int nModelMin);

    /// Generic implementation of RANSAC
    double run(RunResult& res, int nIterMax=1000, bool verbose=false) const;

protected:
    double _cpII;
    int _nModelMin;
    bool satisfyingRun(double runOutput) const;

private:
    // The <name>Ready values need to be set before calling the run method.
    double _precision;


    void find_inliers(const ModelEstimator::Model &model,
                      double precision, std::vector<int> &inliers) const;
};

}  // namespace orsa

#endif
