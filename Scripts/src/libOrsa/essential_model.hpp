/**
* @file essential_model.hpp
* @brief Compute essential matrix
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

#ifndef MMM_ORSA_ESSENTIAL_MODEL_HPP
#define MMM_ORSA_ESSENTIAL_MODEL_HPP

#include "eigen/Eigen/Eigen"
#include "fundamental_model.hpp"
#include "match.hpp"
#include "model_estimator.hpp"


namespace orsa {

/// Convert an eigen matrix to our own matrix defined in libNumerics/matrix.h
libNumerics::matrix<double> eigenToLibMatrix(const Eigen::MatrixXd &matrixToChange);

/// Apply calibration matrix normalisation to the points to go from pixel points to camera points.
void normalisePoints(const std::vector<Match> &points_,
                     const libNumerics::matrix<double> &intrinsicsSrc_,
                     const libNumerics::matrix<double> &intrinsicsDst_,
                     std::vector<Match> &normalisedPoints_);

/// Transforms an essential matrix to a fundamental matrix.
libNumerics::matrix<double>
essentialToFundamental(const libNumerics::matrix<double> &E, const libNumerics::matrix<double> &intrinsicsSrc,
                       const libNumerics::matrix<double> &intrinsicsDst);

/// Transforms a fundamental matrix to an essential matrix.
libNumerics::matrix<double>
fundamentalToEssential(const libNumerics::matrix<double> &F, const libNumerics::matrix<double> &intrinsicsSrc,
                       const libNumerics::matrix<double> &intrinsicsDst);


/// Script extracted from MAGSAC implementation.
// Multiply two degree one polynomials of variables x, y, z.
// E.g. p1 = a[0]x + a[1]y + a[2]z + a[3]
// Output order: x^2 xy y^2 xz yz z^2 x y z 1 (GrevLex)
inline Eigen::Matrix<double, 1, 10> multiplyDegOnePoly(
        const Eigen::RowVector4d &a,
        const Eigen::RowVector4d &b);

/// Script extracted from MAGSAC implementation.
// Multiply a 2 deg poly (in x, y, z) and a one deg poly in GrevLex order.
// x^3 x^2y xy^2 y^3 x^2z xyz y^2z xz^2 yz^2 z^3 x^2 xy y^2 xz yz z^2 x y z 1
inline Eigen::Matrix<double, 1, 20> multiplyDegTwoDegOnePoly(
        const Eigen::Matrix<double, 1, 10> &a,
        const Eigen::RowVector4d &b);

/// Script extracted from MAGSAC implementation.
// Shorthand for multiplying the Essential matrix with its transpose.
inline Eigen::Matrix<double, 1, 10> computeEETranspose(
        const Eigen::Matrix<double, 1, 4> nullSpace[3][3],
        int i,
        int j);

/// Essential matrix model.
class EssentialModel : public FundamentalModel {
public:
    /// Both normalised and un-normalised matches are necessary for computation.
    /// Calibration matrixes are also requiered.
    EssentialModel(const std::vector<Match> &m, const std::vector<Match> &mNormalised,
                   int width1, int height1, int width2, int height2,
                   const libNumerics::matrix<double> &intrinsicsSrc_, const libNumerics::matrix<double> &intrinsicsDst_,
                   bool symError = false);

    /// 5 point correspondences required to compute an essential matrix.
    int SizeSample() const { return 5; }

    /// 10 Essential matrix can be estimated from a sample of 5 points.
    int NbModels() const { return 10; }

    /// Degree of freedom of the model
    int nDegreeOfFreedom() const { return 5;};
    /// Degree of freedom of the data (for magsac++)
    int dataDegreeOfFreedom() const {return 4;};
    /// C(n) = 1/(n/2 * Gamma(n/2)) where n is DoF (for magsac++)
    double C() const {return 0.25;};
    /// UG((DoF - 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
    double upperGamma() const {return 0.0036572608340910764;};
    /// LG((DoF + 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
    double lowerGamma() const {return 1.3012265540498875;};

    /// Distance used to distinguish inlier/outlier is to a line
    virtual bool DistToPoint() const { return false; }

    /// Computes a essential matrix given a set of indeces in the matches.
    void Fit(const std::vector<int> &indices, std::vector<Model> *Es) const;

    void Fit(const std::vector<int> &indices, std::vector<Model> *Es, const double *weights_) const;

    /// Sampson error in pixel space for a given point through E.
    double Error(const Model &E, int index, int *side = 0) const;

    double Error(const Model &E, Mat testMat, int *side = 0) const;

    /// Computes the fundamental matrix associated with E.
    inline Model toPixelSpace(const Model &E) const;

private:
    /// Essential matrix specific elements are added to the fundamental model.
    const Mat _dataNormalised;
    libNumerics::matrix<double> _intrinsicSrc, _invIntrinsicSrc;
    libNumerics::matrix<double> _intrinsicDst, _invIntrinsicDst;
    mutable Model _E;
    mutable Model _F;

    void EpipolarEquation(const std::vector<int> &indices, Mat *A) const;
    void algo8pt(const Mat &A, std::vector<Mat> *Es) const;

    void algo5pt(const std::vector<int> &indices, std::vector<Model> *Es) const;

    void algo5pt(const std::vector<int> &indices, std::vector<Model> *Es, const double *weights_) const;
};

} // Namespace Orsa

#endif //MMM_ORSA_ESSENTIAL_MODEL_HPP
