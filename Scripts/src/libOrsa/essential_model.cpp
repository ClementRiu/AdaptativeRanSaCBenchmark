/** Copyright (C) 2019 Czech Technical University.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following
*       disclaimer in the documentation and/or other materials provided
*       with the distribution.
*
*     * Neither the name of Czech Technical University nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
* Please contact the author of this library if you have any questions.
* Author: Daniel Barath (barath.daniel@sztaki.mta.hu)
*/
// Adaptation by Clément Riu (2021)

#include "essential_model.hpp"
#include "e5p.h"
#include "libNumerics/numerics.h"

namespace orsa {

/// Convert an eigen matrix to our own matrix defined in libNumerics/matrix.h
libNumerics::matrix<double> eigenToLibMatrix(const Eigen::MatrixXd &matrixToChange) {
    int nRow = matrixToChange.rows();
    int nCol = matrixToChange.cols();
    libNumerics::matrix<double> matrixToReturn(nRow, nCol);
    for (int iRow = 0; iRow < nRow; iRow++) {
        for (int iCol = 0; iCol < nCol; iCol++) {
            matrixToReturn(iRow, iCol) = matrixToChange(iRow, iCol);
        }
    }
    return matrixToReturn;
}

/// Apply calibration matrix normalisation to the points to go from pixel points to camera points.
void normalisePoints(const std::vector<Match> &points_,
                     const libNumerics::matrix<double> &intrinsicsSrc_,
                     const libNumerics::matrix<double> &intrinsicsDst_,
                     std::vector<Match> &normalisedPoints_) {

    const libNumerics::matrix<double> inverse_intrinsics_src = intrinsicsSrc_.inv(),
        inverse_intrinsics_dst = intrinsicsDst_.inv();

    // Most likely, this is not the fastest solution, but it does
    // not affect the speed of Graph-cut RANSAC, so not a crucial part of
    // this example.
    double x0, y0, x1, y1;
    std::vector<Match>::const_iterator itPoints = points_.begin();
    for (; itPoints != points_.end(); ++itPoints) {
        libNumerics::vector<double> point_src(3),
            point_dst(3),
            normalized_point_src(3),
            normalized_point_dst(3);

        x0 = (*itPoints).x1;
        y0 = (*itPoints).y1;
        x1 = (*itPoints).x2;
        y1 = (*itPoints).y2;

        point_src(0) = x0;
        point_src(1) = y0;
        point_src(2) = 1.0; // Homogeneous point in the first image
        point_dst(0) = x1;
        point_dst(1) = y1;
        point_dst(2) = 1.0; // Homogeneous point in the second image

        // Normalized homogeneous point in the first image
        normalized_point_src =
            inverse_intrinsics_src * point_src;
        // Normalized homogeneous point in the second image
        normalized_point_dst =
            inverse_intrinsics_dst * point_dst;

        Match normalisedMatch(normalized_point_src(0), normalized_point_src(1), normalized_point_dst(0),
                              normalized_point_dst(1));
        // The second four columns contain the normalized coordinates.
        normalisedPoints_.push_back(normalisedMatch);
    }
}

/// Transforms an essential matrix to a fundamental matrix.
libNumerics::matrix<double>
essentialToFundamental(const libNumerics::matrix<double> &E, const libNumerics::matrix<double> &intrinsicsSrc,
                       const libNumerics::matrix<double> &intrinsicsDst) {
    return intrinsicsSrc.t().inv() * E * intrinsicsDst.inv();
}

/// Transforms a fundamental matrix to an essential matrix.
libNumerics::matrix<double>
fundamentalToEssential(const libNumerics::matrix<double> &F, const libNumerics::matrix<double> &intrinsicsSrc,
                       const libNumerics::matrix<double> &intrinsicsDst) {
    return intrinsicsSrc.t() * F * intrinsicsDst;
}

/// Constructor
EssentialModel::EssentialModel(const std::vector<Match> &m, const std::vector<Match> &mNormalised,
                               int width1, int height1, int width2, int height2,
                               const libNumerics::matrix<double> &intrinsicsSrc_,
                               const libNumerics::matrix<double> &intrinsicsDst_,
                               bool symError) : FundamentalModel(m, width1, height1, width2, height2, symError),
                                                _dataNormalised(Match::toMat(mNormalised)),
                                                _intrinsicSrc(intrinsicsSrc_), _intrinsicDst(intrinsicsDst_) {
    _invIntrinsicSrc = intrinsicsSrc_.inv();
    _invIntrinsicDst = intrinsicsDst_.inv();
}

/// Script extracted from MAGSAC implementation.
/// Computes aN essential matrix given a set of indices in matches.
void EssentialModel::Fit(const std::vector<int> &indices, std::vector<Model> *Es) const {
    if (indices.size() >= 8) {
        std::vector<Model> Fs;
        FundamentalModel::Fit(indices, &Fs);
        std::vector<Model>::const_iterator itF = Fs.begin();
        for (; itF != Fs.end(); itF++){
            Model E = fundamentalToEssential(*itF, _intrinsicSrc, _intrinsicDst);
            typedef Eigen::Matrix<double,3,3,Eigen::RowMajor> Mat;
            Eigen::Map<Mat> E2eigen(&E(0));
            Eigen::JacobiSVD<Mat> svd(E2eigen,
                                      Eigen::ComputeFullU|Eigen::ComputeFullV);
            Eigen::Vector<double,3> D = svd.singularValues();
            D(0)=D(1)=1; D(2)=0;
            E2eigen = svd.matrixU()*D.asDiagonal()*svd.matrixV().transpose();
            Es->push_back(eigenToLibMatrix(E2eigen));

        }
    } else {
        algo5pt(indices, Es);
    }
}

// TODO
void EssentialModel::Fit(const std::vector<int> &indices, std::vector<Model> *Es, const double *weights_) const {
    if (indices.size() >= 8) {
        std::vector<Model> Fs;
        FundamentalModel::Fit(indices, &Fs, weights_);
        std::vector<Model>::const_iterator itF = Fs.begin();
        for (; itF != Fs.end(); itF++){
            Model E = fundamentalToEssential(*itF, _intrinsicSrc, _intrinsicDst);
            typedef Eigen::Matrix<double,3,3,Eigen::RowMajor> Mat;
            Eigen::Map<Mat> E2eigen(&E(0));
            Eigen::JacobiSVD<Mat> svd(E2eigen,
                                      Eigen::ComputeFullU|Eigen::ComputeFullV);
            Eigen::Vector<double,3> D = svd.singularValues();
            D(0)=D(1)=1; D(2)=0;
            E2eigen = svd.matrixU()*D.asDiagonal()*svd.matrixV().transpose();
            Es->push_back(eigenToLibMatrix(E2eigen));

        }
    } else {
        algo5pt(indices, Es, weights_);
    }
}

void EssentialModel::EpipolarEquation(const std::vector<int> &indices,
                                      ModelEstimator::Mat *A) const {
    for (size_t i = 0; i < indices.size(); ++i) {
        int j = indices[i];
        libNumerics::vector<double> x1(3), x2(3);
        x1(0) = _dataNormalised(0, j);
        x1(1) = _dataNormalised(1, j);
        x1(2) = 1;
        x2(0) = _dataNormalised(2, j);
        x2(1) = _dataNormalised(3, j);
        x2(2) = 1;
        (*A)(i, 0) = x1(0) * x2(0);  // 0 represents x coords,
        (*A)(i, 1) = x1(0) * x2(1);  // 1 represents y coords.
        (*A)(i, 2) = x1(0);
        (*A)(i, 3) = x1(1) * x2(0);
        (*A)(i, 4) = x1(1) * x2(1);
        (*A)(i, 5) = x1(1);
        (*A)(i, 6) = x2(0);
        (*A)(i, 7) = x2(1);
        (*A)(i, 8) = 1.0;
    }
}

void EssentialModel::algo8pt(const Mat &A, std::vector<Mat> *Es) const {
    // Without rank constraint
    libNumerics::vector<double> vecNullspace(9);
    libNumerics::SVD::Nullspace(A, &vecNullspace, 1, 1);
    libNumerics::matrix<double> E(3, 3);
    E.read(vecNullspace);

    // Force the rank 2 constraint
    typedef Eigen::Matrix<double,3,3,Eigen::RowMajor> Mat;
    Eigen::Map<Mat> E2eigen(&E(0));
    Eigen::JacobiSVD<Mat> svd(E2eigen,
                              Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::Vector<double,3> D = svd.singularValues();
    D(0)=D(1)=1; D(2)=0;
    E2eigen = svd.matrixU()*D.asDiagonal()*svd.matrixV().transpose();
    Es->push_back(eigenToLibMatrix(E2eigen));
}

void EssentialModel::algo5pt(const std::vector<int> &indices, std::vector<Model> *Es) const {
    Eigen::MatrixXd coefficients(indices.size(), 9);

    // Step 1. Create the nx9 matrix containing epipolar constraints.
    //   Essential matrix is a linear combination of the 4 vectors spanning the null space of this
    //   matrix.
    double x0, y0, x1, y1, weight = 1.0;
    for (size_t i = 0; i < indices.size(); i++) {
        int index = indices[i];
        x0 = _dataNormalised(0, index);
        y0 = _dataNormalised(1, index);
        x1 = _dataNormalised(2, index);
        y1 = _dataNormalised(3, index);

        // Precalculate these values to avoid calculating them multiple times
        const double
                weight_times_x0 = weight * x0,
                weight_times_x1 = weight * x1,
                weight_times_y0 = weight * y0,
                weight_times_y1 = weight * y1;

        coefficients.row(i) <<
                weight_times_x0 * x1,
                weight_times_x0 * y1,
                weight_times_x0,
                weight_times_y0 * x1,
                weight_times_y0 * y1,
                weight_times_y0,
                weight_times_x1,
                weight_times_y1,
                weight;
    }

    // Extract the null space from a minimal sampling (using LU) or non-minimal sampling (using SVD).
    Eigen::Matrix<double, 9, 4> nullSpace;

    if (indices.size() == 5) {
        const Eigen::FullPivLU<Eigen::MatrixXd> lu(coefficients);
        if (lu.dimensionOfKernel() != 4) {
            return;
        }
        nullSpace = lu.kernel();
    } else {
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(
                coefficients.transpose() * coefficients);
        const Eigen::MatrixXd &Q = qr.matrixQ();
        nullSpace = Q.rightCols<4>();
    }

    const Eigen::Matrix<double, 1, 4> nullSpaceMatrix[3][3] = {
            {nullSpace.row(0), nullSpace.row(3), nullSpace.row(6)},
            {nullSpace.row(1), nullSpace.row(4), nullSpace.row(7)},
            {nullSpace.row(2), nullSpace.row(5), nullSpace.row(8)}};

    // Step 2. Expansion of the epipolar constraints on the determinant and trace.
    const Eigen::Matrix<double, 10, 20> constraintMatrix = buildConstraintMatrix(nullSpaceMatrix);

    // Step 3. Eliminate part of the matrix to isolate polynomials in z.
    Eigen::FullPivLU<Eigen::Matrix<double, 10, 10>> c_lu(constraintMatrix.block<10, 10>(0, 0));
    const Eigen::Matrix<double, 10, 10> eliminatedMatrix = c_lu.solve(constraintMatrix.block<10, 10>(0, 10));


    Eigen::Matrix<double, 10, 10> actionMatrix = Eigen::Matrix<double, 10, 10>::Zero();
    actionMatrix.block<3, 10>(0, 0) = eliminatedMatrix.block<3, 10>(0, 0);
    actionMatrix.row(3) = eliminatedMatrix.row(4);
    actionMatrix.row(4) = eliminatedMatrix.row(5);
    actionMatrix.row(5) = eliminatedMatrix.row(7);
    actionMatrix(6, 0) = -1.0;
    actionMatrix(7, 1) = -1.0;
    actionMatrix(8, 3) = -1.0;
    actionMatrix(9, 6) = -1.0;
    Eigen::EigenSolver<Eigen::Matrix<double, 10, 10>> eigensolver(actionMatrix);
    const Eigen::VectorXcd &eigenvalues = eigensolver.eigenvalues();

    // Now that we have x, y, and z we need to substitute them back into the null space to get a valid
    // essential matrix solution.
    for (size_t i = 0; i < 10; i++) {
        // Only consider real solutions.
        if (eigenvalues(i).imag() != 0) {
            continue;
        }

        Eigen::Matrix3d E_dst_src;
        Eigen::Map<Eigen::Matrix<double, 9, 1>>(E_dst_src.data()) =
                nullSpace * eigensolver.eigenvectors().col(i).tail<4>().real();

        Eigen::MatrixXd model;
        model = E_dst_src;
        Es->push_back(eigenToLibMatrix(model).t());
    }

    return;
}

    void EssentialModel::algo5pt(const std::vector<int> &indices,
                                 std::vector<Model> *Es,
                                 const double *weights_) const {
        Eigen::MatrixXd coefficients(indices.size(), 9);

        // Step 1. Create the nx9 matrix containing epipolar constraints.
        //   Essential matrix is a linear combination of the 4 vectors spanning the null space of this
        //   matrix.
        double x0, y0, x1, y1, weight = 1.0;
        for (size_t i = 0; i < indices.size(); i++) {
            int index = indices[i];

            if (weights_ != nullptr)
                weight = weights_[index];

            x0 = _dataNormalised(0, index);
            y0 = _dataNormalised(1, index);
            x1 = _dataNormalised(2, index);
            y1 = _dataNormalised(3, index);

            // Precalculate these values to avoid calculating them multiple times
            const double
                    weight_times_x0 = weight * x0,
                    weight_times_x1 = weight * x1,
                    weight_times_y0 = weight * y0,
                    weight_times_y1 = weight * y1;

            coefficients.row(i) <<
                    weight_times_x0 * x1,
                    weight_times_x0 * y1,
                    weight_times_x0,
                    weight_times_y0 * x1,
                    weight_times_y0 * y1,
                    weight_times_y0,
                    weight_times_x1,
                    weight_times_y1,
                    weight;
        }

        // Extract the null space from a minimal sampling (using LU) or non-minimal sampling (using SVD).
        Eigen::Matrix<double, 9, 4> nullSpace;

        if (indices.size() == 5) {
            const Eigen::FullPivLU<Eigen::MatrixXd> lu(coefficients);
            if (lu.dimensionOfKernel() != 4) {
                return;
            }
            nullSpace = lu.kernel();
        } else {
            const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(
                    coefficients.transpose() * coefficients);
            const Eigen::MatrixXd &Q = qr.matrixQ();
            nullSpace = Q.rightCols<4>();
        }

        const Eigen::Matrix<double, 1, 4> nullSpaceMatrix[3][3] = {
                {nullSpace.row(0), nullSpace.row(3), nullSpace.row(6)},
                {nullSpace.row(1), nullSpace.row(4), nullSpace.row(7)},
                {nullSpace.row(2), nullSpace.row(5), nullSpace.row(8)}};

        // Step 2. Expansion of the epipolar constraints on the determinant and trace.
        const Eigen::Matrix<double, 10, 20> constraintMatrix = buildConstraintMatrix(nullSpaceMatrix);

        // Step 3. Eliminate part of the matrix to isolate polynomials in z.
        Eigen::FullPivLU<Eigen::Matrix<double, 10, 10>> c_lu(constraintMatrix.block<10, 10>(0, 0));
        const Eigen::Matrix<double, 10, 10> eliminatedMatrix = c_lu.solve(constraintMatrix.block<10, 10>(0, 10));


        Eigen::Matrix<double, 10, 10> actionMatrix = Eigen::Matrix<double, 10, 10>::Zero();
        actionMatrix.block<3, 10>(0, 0) = eliminatedMatrix.block<3, 10>(0, 0);
        actionMatrix.row(3) = eliminatedMatrix.row(4);
        actionMatrix.row(4) = eliminatedMatrix.row(5);
        actionMatrix.row(5) = eliminatedMatrix.row(7);
        actionMatrix(6, 0) = -1.0;
        actionMatrix(7, 1) = -1.0;
        actionMatrix(8, 3) = -1.0;
        actionMatrix(9, 6) = -1.0;
        Eigen::EigenSolver<Eigen::Matrix<double, 10, 10>> eigensolver(actionMatrix);
        const Eigen::VectorXcd &eigenvalues = eigensolver.eigenvalues();

        // Now that we have x, y, and z we need to substitute them back into the null space to get a valid
        // essential matrix solution.
        for (size_t i = 0; i < 10; i++) {
            // Only consider real solutions.
            if (eigenvalues(i).imag() != 0) {
                continue;
            }

            Eigen::Matrix3d E_dst_src;
            Eigen::Map<Eigen::Matrix<double, 9, 1>>(E_dst_src.data()) =
                    nullSpace * eigensolver.eigenvectors().col(i).tail<4>().real();

            Eigen::MatrixXd model;
            model = E_dst_src;
            Es->push_back(eigenToLibMatrix(model).t());
        }

        return;
    }


/// Computes the fundamental matrix associated with E.
inline EssentialModel::Model EssentialModel::toPixelSpace(const Model &E) const {
    if (E != _E) {
        _E = E;
        _F = essentialToFundamental(_E, _intrinsicSrc, _intrinsicDst);
    }
    return _F;
}

/// Sampson error in pixel space for a given point through E.
/// BEWARE : this takes the fundamental matrix as input.
/// \param E The essential matrix.
/// \param index The point correspondence.
/// \param side In which image is the error measured?
/// \return The square reprojection error.
double EssentialModel::Error(const Model &E, int index, int *side) const {
    double xa = _dataNormalised(0, index), ya = _dataNormalised(1, index);
    double xb = data_(2, index), yb = data_(3, index);

    double a, b, c, a2, b2, c2, d;
    // Transfer error in image 2
    if (side) *side = 1;
    a = E(0, 0) * xa + E(1, 0) * ya + E(2, 0);
    b = E(0, 1) * xa + E(1, 1) * ya + E(2, 1);
    c = E(0, 2) * xa + E(1, 2) * ya + E(2, 2);
    a2 = _invIntrinsicDst(0,0)*a + _invIntrinsicDst(1,0)*b + _invIntrinsicDst(2,0)*c;
    b2 = _invIntrinsicDst(0,1)*a + _invIntrinsicDst(1,1)*b + _invIntrinsicDst(2,1)*c;
    c2 = _invIntrinsicDst(0,2)*a + _invIntrinsicDst(1,2)*b + _invIntrinsicDst(2,2)*c;

    d = a2 * xb + b2 * yb + c2;
    double err = (d * d) / (a2 * a2 + b2 * b2);
    // Transfer error in image 1
    if (symError) { // ... but only if requested
        xa = data_(0, index), ya = data_(1, index);
        xb = _dataNormalised(2, index), yb = _dataNormalised(3, index);
        a = E(0, 0) * xb + E(0, 1) * yb + E(0, 2);
        b = E(1, 0) * xb + E(1, 1) * yb + E(1, 2);
        c = E(2, 0) * xb + E(2, 1) * yb + E(2, 2);
        a2 = _invIntrinsicSrc(0,0)*a + _invIntrinsicSrc(1,0)*b + _invIntrinsicSrc(2,0)*c;
        b2 = _invIntrinsicSrc(0,1)*a + _invIntrinsicSrc(1,1)*b + _invIntrinsicSrc(2,1)*c;
        c2 = _invIntrinsicSrc(0,2)*a + _invIntrinsicSrc(1,2)*b + _invIntrinsicSrc(2,2)*c;
        d = a2 * xa + b2 * ya + c2;
        double err1 = (d * d) / (a2 * a2 + b2 * b2);
        if (err1 > err) {
            err = err1;
            if (side) *side = 0;
        }
    }
    return err;
}

/// Sampson error in pixel space for a given point through E.
/// BEWARE : this takes the essential matrix as input.
/// \param E The essential matrix.
/// \param testMat The point correspondence.
/// \param side In which image is the error measured?
/// \return The square reprojection error.
double EssentialModel::Error(const Model &E, Mat testMat, int *side) const {
    const Model F = toPixelSpace(E);
    double xa = testMat(0, 0), ya = testMat(1, 0);
    double xb = testMat(2, 0), yb = testMat(3, 0);

    return this->FundamentalModel::Error(F, xa, xb, ya, yb, side);
}

}  // namespace orsa
