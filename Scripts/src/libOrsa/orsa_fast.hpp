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

#ifndef ORSA_FAST_H
#define ORSA_FAST_H

#include "orsa.hpp"

namespace orsa {

/// Model estimation with ORSA algorithm.
class OrsaFast : public Orsa {
public:
    /// Constructor
    explicit OrsaFast(const ModelEstimator *estimator);

    /// Generic implementation of ORSA (Optimized Random Sampling Algorithm)
    double run(RunResult& res, int nIterMax=1000, bool verbose=false) const;
protected:
    std::pair<double, double> bestNFA(const std::vector<double> &e,
                       double loge0, double maxThreshold,
                       const std::vector<float> &vLogc_n,
                       const std::vector<float> &vLogc_k) const;

};

}  // namespace orsa

#endif
