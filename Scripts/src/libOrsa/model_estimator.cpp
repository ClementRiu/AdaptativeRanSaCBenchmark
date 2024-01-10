/**
 * @file model_estimator.cpp
 * @brief Model estimation by ORSA (aka AC-RANSAC) algorithm
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 *
 * Copyright (c) 2010-2011,2020-2021 Pascal Monasse
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

#include "libOrsa/model_estimator.hpp"

namespace orsa {

/// Matrix \a data is mxn, representing n datapoints of dimension m.
ModelEstimator::ModelEstimator(const Mat &data, bool symmetric)
: symError(symmetric), data_(data) {}

/// If multiple solutions are possible, return false.
bool ModelEstimator::ComputeModel(const std::vector<int> &indices,
                                  Model *model) const {
    std::vector<Model> models;
    Fit(indices, &models);
    if (models.size() != 1)
        return false;
    *model = models.front();
    return true;
}

ModelEstimator::Mat ModelEstimator::data() const {
    return data_;
}

/// Find inliers of \a model passed as parameter (within \a precision).
void ModelEstimator::FindInliers(const Model &model, double precision,
                                 std::vector<int> &inliers) const {
    precision = precision*precision;
    const int nData = NbData();
    for (int i = 0; i < nData; i++) {
        if (Error(model, i) <= precision) {
            inliers.push_back(i);
        }
    }
}

/// Compute residuals for all data points.
void ModelEstimator::computeResiduals(const Model &model, std::vector<double> &residuals) const {
    for (int i = 0; i < NbData(); i++) {
        residuals.push_back(Error(model, i));
    }
}

/// RMSE/max error of inliers of model \a M.
std::pair<double, double>
ModelEstimator::ErrorStats(const std::vector<int> &in, const Model &M) const {
    std::vector<int>::const_iterator it = in.begin();
    double l2 = 0, linf = 0;
    for (; it != in.end(); ++it) {
        double e = Error(M, *it);
        l2 += e;
        if (linf < e)
            linf = e;
    }
    return std::pair<double, double>(sqrt(l2 / in.size()), sqrt(linf));
}

} // namespace orsa
