/**
 * @file model_estimator.hpp
 * @brief Model regression from sample matches.
 * @author Pascal Monasse, Pierre Moulon
 *
 * Copyright (c) 2011 Pierre Moulon
 * Copyright (c) 2011,2020-2021 Pascal Monasse
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

#ifndef MODEL_ESTIMATOR_H
#define MODEL_ESTIMATOR_H

#include <vector>
#include <utility>
#include "libNumerics/matrix.h"

namespace orsa {

/// Generic class for parametric model estimation.
///
/// Subclasses must follow this interface:
///   1. SizeSample() Number correspondences necessary to compute a model.
///   2. NbModels() Number models from sample of SizeSample() correspondences.
///   3. DistToPoint() Residual is distance to a point or to a line?
///   4. Fit(const vector<size_t> &indices, vector<Model> *models)
///        Compute model(s) compatible with indexed correspondences.
///   5. Error(const Model &model, size_t index)
///        Reprojection square error for indexed correspondence.
///   6. pSigma(double sigma)
///        Compute the area of the inlier region (for LRTSAC)
///   7. toPixelSpace(const model &E)
///        Returns the associated fundamental matrix for the essential model.
///        Otherwise it just passes the input.
/// Points 2 and 3 are used only by ORSA, but they still need to be defined.
class ModelEstimator {
public:
    typedef libNumerics::matrix<double> Mat;
    typedef Mat Model;

    /// With true \c symError, the error is the maximum of dist(M(x1),x2) (side 1)
    /// and dist(x1,M'(x2)) (side 0). M' is the dual model of M (M'=M^-1 for
    /// homography, M'=M^T for fundamental).
    /// Return in \a side of method Error (if non-null pointer) the side this
    /// maximum is reached.
    bool symError;

    /// Constructor
    ModelEstimator(const Mat &data, bool symmetricError = false);
    virtual ~ModelEstimator() {}

    /// Number of data matches.
    int NbData() const { return data_.ncol(); }

    /// Compute model from points.
    bool ComputeModel(const std::vector<int> &indices, Model *model) const;

    /// Minimum number of points required to compute a model.
    /// - homography -> 4
    /// - fundamental 7 pts -> 7
    /// - fundamental 8 pts -> 8
    virtual int SizeSample() const = 0;

    /// Maximum number of models possibly computed from a sample.
    /// - homography -> 1
    /// - fundamental 7 pts -> 3
    /// - fundamental 8 pts -> 1
    virtual int NbModels() const = 0;

    /// Degree of freedom of the model
    virtual int nDegreeOfFreedom() const = 0;
    /// Degree of freedom of the data (for magsac++)
    virtual int dataDegreeOfFreedom() const = 0;
    /// C(n) = 1/(n/2 * Gamma(n/2)) where n is DoF (for magsac++)
    virtual double C() const = 0;
    /// UG((DoF - 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
    virtual double upperGamma() const = 0;
    /// LG((DoF + 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
    virtual double lowerGamma() const =0;

    /// Indicate if distance used to distinguish inlier/outlier is to a point
    /// (true) or a line (false). Most regression routines, such as RANSAC, do not
    /// care about that, bur ORSA needs it.
    /// - homography -> true
    /// - fundamental -> false
    virtual bool DistToPoint() const = 0;

    /// Computes the square error of a correspondence wrt \a model.
    /// \param model The model to evaluate.
    /// \param index The point index stored in \c data.
    /// \param[out] side In which image is the error measured?
    virtual double Error(const Model &model, int index, int *side = 0) const=0;
    /// Computes the square error of a correspondence wrt \a model.
    /// \param model The model to evaluate.
    /// \param testData a point to compute the error on.
    /// \param[out] side In which image is the error measured?
    virtual double Error(const Model &model, Mat testData, int *side=0) const=0;

    /// Computes the models associated to indexed sample.
    /// \param indices Indices of points to consider for model estimation.
    /// \param models  Estimated model(s) from sampled point.
    virtual void Fit(const std::vector<int> &indices,
                     std::vector<Model> *models) const = 0;

    /// Computes the models associated to indexed sample.
    /// \param indices Indices of points to consider for model estimation.
    /// \param models  Estimated model(s) from sampled point.
    /// \param models  Estimated model(s) from sampled point. //TODO
    virtual void Fit(const std::vector<int> &indices,
                     std::vector<Model> *models, const double *weights_) const = 0;

    /// Computes the area of the inlier region according to the AC-RANSAC
    /// paper for a given error margin. The inlier region is over-estimated to
    /// ease computation.
    /// \param sigma The error margin.
    virtual double pSigma(double sigma, bool leftSide = true) const = 0;

    /// Returns the data stored in the estimator. Useful to access
    /// artificially created data.
    Mat data() const;

    /// Returns the fundamental model from an essential matrix.
    /// For other models it just returns the input.
    virtual Model toPixelSpace(const Model &M) const = 0;

    /// Find inliers of \a model within \a precision.
    void FindInliers(const Model &model, double precision,
                      std::vector<int> &inliers) const;

    /// Compute residuals for all data points.
    void computeResiduals(const Model &model,
                          std::vector<double> &residuals) const;
    /// RMSE/Max error for inliers \a in and model \a M.
    std::pair<double, double>
    ErrorStats(const std::vector<int> &in, const Model &M) const;
protected:
    Mat data_; ///< Data to estimate a model
};

}  // namespace orsa

#endif
