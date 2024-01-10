// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

//
// Adaptation by Clément Riu (2021)
//

#include "pnp_model.hpp"
#include "pnp_poly.hpp"
#include "essential_model.hpp"

namespace orsa {
/// Apply calibration matrix normalisation to the points to go from pixel points to camera points.
void normalisePoints(const std::vector<Match2D3D> &points_,
                     const libNumerics::matrix<double> &calibration,
                     std::vector<Match2D3D> &normalisedPoints_) {
    const libNumerics::matrix<double> inverse_calibration = calibration.inv();

    // Most likely, this is not the fastest solution, but it does
    // not affect the speed of Graph-cut RANSAC, so not a crucial part of
    // this example.
    double x0, y0;
    std::vector<Match2D3D>::const_iterator itPoints = points_.begin();
    for (; itPoints != points_.end(); ++itPoints) {
        libNumerics::vector<double> point_src(3),
                normalized_point_src(3);

        x0 = (*itPoints).p2D(0);
        y0 = (*itPoints).p2D(1);

        point_src(0) = x0;
        point_src(1) = y0;
        point_src(2) = 1.0; // Homogeneous point in the first image

        // Normalized homogeneous point in the first image
        normalized_point_src =
                inverse_calibration * point_src;

        Match2D3D normalisedMatch(normalized_point_src(0), normalized_point_src(1),
                                  (*itPoints).p3D(0), (*itPoints).p3D(1), (*itPoints).p3D(2));
        // The second four columns contain the normalized coordinates.
        normalisedPoints_.push_back(normalisedMatch);
    }
}

void convert3Dto2D(const ModelEstimator::Model &proj_matrix,
                   const libNumerics::matrix<double> &intrinsic,
                   Match2D3D &m) {
    double r11 = proj_matrix(0, 0), r12 = proj_matrix(0, 1), r13 = proj_matrix(0, 2),
            r21 = proj_matrix(1, 0), r22 = proj_matrix(1, 1), r23 = proj_matrix(1, 2),
            r31 = proj_matrix(2, 0), r32 = proj_matrix(2, 1), r33 = proj_matrix(2, 2),
            t1 = proj_matrix(0, 3), t2 = proj_matrix(1, 3), t3 = proj_matrix(2, 3);
    double x = r11 * m.p3D(0) + r12 * m.p3D(1) + r13 * m.p3D(2) + t1;
    double y = r21 * m.p3D(0) + r22 * m.p3D(1) + r23 * m.p3D(2) + t2;
    double z = r31 * m.p3D(0) + r32 * m.p3D(1) + r33 * m.p3D(2) + t3;
    double u = intrinsic(0, 0) * x + intrinsic(0, 1) * y + intrinsic(0, 2) * z;
    double v = intrinsic(1, 1) * y + intrinsic(1, 2) * z;
    if (z > 0){
        m.p2D(0) = u / z;
        m.p2D(1) = v / z;
    } else {
        m.p2D(0) = - 1;
        m.p2D(1) = - 1;
    }
}

Eigen::Vector3d LiftImagePoint(const Eigen::Vector2d &point) {
    return point.homogeneous() / std::sqrt(point.squaredNorm() + 1);
}

ModelEstimator::Mat PnPModel::RTtoMat(const Eigen::Matrix3d &R, const Eigen::Vector3d &T) {
    Mat matrixToReturn(3, 4);
    for (int iRow = 0; iRow < 3; iRow++) {
        for (int iCol = 0; iCol < 3; iCol++) {
            matrixToReturn(iRow, iCol) = R(iRow, iCol);
        }
        matrixToReturn(iRow, 3) = T(iRow);
    }
    return matrixToReturn;
}

PnPModel::PnPModel(const std::vector<Match2D3D> &m, const std::vector<Match2D3D> &mNormalised,
                   const libNumerics::matrix<double> &intrinsics,
                   const int width, const int height) :
        ModelEstimator(Match2D3D::toMat(m), false),
        _area(width * height), _matches(m), _matchesNormalised(mNormalised),
        _intrinsic(intrinsics) {
    _invIntrinsic = intrinsics.inv();
}

void PnPModel::Fit(const std::vector<int> &indices, std::vector<Model> *RTs) const {
    if (indices.size() > 3) {
        EPNP(indices, RTs);
    } else {
        P3P(indices, RTs);
    }
}

void PnPModel::Fit(const std::vector<int> &indices, std::vector<Model> *RTs, const double *weights_) const {
    if (indices.size() > 3) {
        EPNP(indices, RTs, weights_);
    } else {
        P3P(indices, RTs);
    }
}

void PnPModel::EPNP(const std::vector<int> &indices, std::vector<Model> *RTs, const double *weights_) const {
    Model proj_matrix;
    if (!ComputePose(&proj_matrix, indices, weights_)) {
        return;
    }

    RTs->push_back(proj_matrix);
}

void PnPModel::P3P(const std::vector<int> &indices, std::vector<Model> *RTs) const {
    Eigen::Matrix3d points3D_world;
    points3D_world.col(0) = _matchesNormalised[indices[0]].p3D;
    points3D_world.col(1) = _matchesNormalised[indices[1]].p3D;
    points3D_world.col(2) = _matchesNormalised[indices[2]].p3D;

    const Eigen::Vector3d u = LiftImagePoint(_matchesNormalised[indices[0]].p2D);
    const Eigen::Vector3d v = LiftImagePoint(_matchesNormalised[indices[1]].p2D);
    const Eigen::Vector3d w = LiftImagePoint(_matchesNormalised[indices[2]].p2D);

    // Angles between 2D points.
    const double cos_uv = u.transpose() * v;
    const double cos_uw = u.transpose() * w;
    const double cos_vw = v.transpose() * w;

    // Distances between 2D points.
    const double dist_AB_2 = (_matchesNormalised[indices[0]].p3D - _matchesNormalised[indices[1]].p3D).squaredNorm();
    const double dist_AC_2 = (_matchesNormalised[indices[0]].p3D - _matchesNormalised[indices[2]].p3D).squaredNorm();
    const double dist_BC_2 = (_matchesNormalised[indices[1]].p3D - _matchesNormalised[indices[2]].p3D).squaredNorm();

    const double dist_AB = std::sqrt(dist_AB_2);

    const double a = dist_BC_2 / dist_AB_2;
    const double b = dist_AC_2 / dist_AB_2;

    // Helper variables for calculation of coefficients.
    const double a2 = a * a;
    const double b2 = b * b;
    const double p = 2 * cos_vw;
    const double q = 2 * cos_uw;
    const double r = 2 * cos_uv;
    const double p2 = p * p;
    const double p3 = p2 * p;
    const double q2 = q * q;
    const double r2 = r * r;
    const double r3 = r2 * r;
    const double r4 = r3 * r;
    const double r5 = r4 * r;

    // Build polynomial coefficients: a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0.
    Eigen::Matrix<double, 5, 1> coeffs;
    coeffs(0) = -2 * b + b2 + a2 + 1 + a * b * (2 - r2) - 2 * a;
    coeffs(1) = -2 * q * a2 - r * p * b2 + 4 * q * a + (2 * q + p * r) * b +
                (r2 * q - 2 * q + r * p) * a * b - 2 * q;
    coeffs(2) = (2 + q2) * a2 + (p2 + r2 - 2) * b2 - (4 + 2 * q2) * a -
                (p * q * r + p2) * b - (p * q * r + r2) * a * b + q2 + 2;
    coeffs(3) = -2 * q * a2 - r * p * b2 + 4 * q * a +
                (p * r + q * p2 - 2 * q) * b + (r * p + 2 * q) * a * b - 2 * q;
    coeffs(4) = a2 + b2 - 2 * a + (2 - p2) * b - 2 * a * b + 1;

    Eigen::VectorXd roots_real;
    Eigen::VectorXd roots_imag;
    if (!FindPolynomialRootsCompanionMatrix(coeffs, &roots_real, &roots_imag)) {
        return;
    }

    RTs->reserve(roots_real.size());

    for (Eigen::VectorXd::Index i = 0; i < roots_real.size(); ++i) {
        const double kMaxRootImag = 1e-10;
        if (std::abs(roots_imag(i)) > kMaxRootImag) {
            continue;
        }

        const double x = roots_real(i);
        if (x < 0) {
            continue;
        }

        const double x2 = x * x;
        const double x3 = x2 * x;

        // Build polynomial coefficients: b1*y + b0 = 0.
        const double bb1 =
                (p2 - p * q * r + r2) * a + (p2 - r2) * b - p2 + p * q * r - r2;
        const double b1 = b * bb1 * bb1;
        const double b0 =
                ((1 - a - b) * x2 + (a - 1) * q * x - a + b + 1) *
                (r3 * (a2 + b2 - 2 * a - 2 * b + (2 - r2) * a * b + 1) * x3 +
                 r2 *
                 (p + p * a2 - 2 * r * q * a * b + 2 * r * q * b - 2 * r * q -
                  2 * p * a - 2 * p * b + p * r2 * b + 4 * r * q * a +
                  q * r3 * a * b - 2 * r * q * a2 + 2 * p * a * b + p * b2 -
                  r2 * p * b2) *
                 x2 +
                 (r5 * (b2 - a * b) - r4 * p * q * b +
                  r3 * (q2 - 4 * a - 2 * q2 * a + q2 * a2 + 2 * a2 - 2 * b2 + 2) +
                  r2 * (4 * p * q * a - 2 * p * q * a * b + 2 * p * q * b - 2 * p * q -
                        2 * p * q * a2) +
                  r * (p2 * b2 - 2 * p2 * b + 2 * p2 * a * b - 2 * p2 * a + p2 +
                       p2 * a2)) *
                 x +
                 (2 * p * r2 - 2 * r3 * q + p3 - 2 * p2 * q * r + p * q2 * r2) * a2 +
                 (p3 - 2 * p * r2) * b2 +
                 (4 * q * r3 - 4 * p * r2 - 2 * p3 + 4 * p2 * q * r - 2 * p * q2 * r2) *
                 a +
                 (-2 * q * r3 + p * r4 + 2 * p2 * q * r - 2 * p3) * b +
                 (2 * p3 + 2 * q * r3 - 2 * p2 * q * r) * a * b + p * q2 * r2 -
                 2 * p2 * q * r + 2 * p * r2 + p3 - 2 * r3 * q);

        // Solve for y.
        const double y = b0 / b1;
        const double y2 = y * y;

        const double nu = x2 + y2 - 2 * x * y * cos_uv;

        const double dist_PC = dist_AB / std::sqrt(nu);
        const double dist_PB = y * dist_PC;
        const double dist_PA = x * dist_PC;

        Eigen::Matrix3d points3D_camera;
        points3D_camera.col(0) = u * dist_PA;  // A'
        points3D_camera.col(1) = v * dist_PB;  // B'
        points3D_camera.col(2) = w * dist_PC;  // C'

        // Find transformation from the world to the camera system.
        const Eigen::Matrix4d transform =
                Eigen::umeyama(points3D_world, points3D_camera, false);
        RTs->push_back(eigenToLibMatrix(transform.topLeftCorner<3, 4>()));
    }
}

double PnPModel::Error(const Model &RT, int index, int *side) const {
    double xa = _matches[index].p2D(0), ya = _matches[index].p2D(1);
    double xb = _matches[index].p3D(0), yb = _matches[index].p3D(1), zb = _matches[index].p3D(2);
    return Error(RT, xa, xb, ya, yb, zb, side);
}

double PnPModel::Error(const Model &RT, Mat testData, int *side) const {
    double xa = testData(0, 0), ya = testData(1, 0);
    double xb = testData(2, 0), yb = testData(3, 0), zb = testData(4, 0);
    return Error(RT, xa, xb, ya, yb, zb, side);
}

/// Square reprojection error for a given point given R and T.
double PnPModel::Error(const Model &RT,
                       const double xa, const double xb,
                       const double ya, const double yb,
                       const double zb, int *side) const {
    // Note that this code might not be as nice as Eigen expressions,
    // but it is significantly faster in various tests.
    if (side) *side = 0;
    const double P_00 = RT(0, 0);
    const double P_01 = RT(0, 1);
    const double P_02 = RT(0, 2);
    const double P_03 = RT(0, 3);
    const double P_10 = RT(1, 0);
    const double P_11 = RT(1, 1);
    const double P_12 = RT(1, 2);
    const double P_13 = RT(1, 3);
    const double P_20 = RT(2, 0);
    const double P_21 = RT(2, 1);
    const double P_22 = RT(2, 2);
    const double P_23 = RT(2, 3);

    // Project 3D point from world to camera.
    const double px_2 = P_20 * xb + P_21 * yb + P_22 * zb + P_23;

    // Check if 3D point is in front of camera.
    if (px_2 > std::numeric_limits<double>::epsilon()) {
        const double px_0 = P_00 * xb + P_01 * yb + P_02 * zb + P_03;
        const double px_1 = P_10 * xb + P_11 * yb + P_12 * zb + P_13;

        const double xc = px_0 / px_2;
        const double yc = px_1 / px_2;

        const double u = _intrinsic(0, 0) * xc + _intrinsic(0, 1) * yc + _intrinsic(0, 2);
        const double v = _intrinsic(1, 1) * yc + _intrinsic(1, 2);

        const double dx_0 = xa - u;
        const double dx_1 = ya - v;

        return dx_0 * dx_0 + dx_1 * dx_1;
    } else {
        return std::numeric_limits<double>::max();
    }

}

void PnPModel::ComputeSquaredReprojectionError(
        const Model &proj_matrix, std::vector<double> *residuals) const {

    residuals->resize(NbData());

    for (int i = 0; i < NbData(); ++i) {
        (*residuals)[i] = Error(proj_matrix, i);
    }
}

void PnPModel::Residuals(const Model &proj_matrix,
                         std::vector<double> *residuals) const {
    ComputeSquaredReprojectionError(proj_matrix, residuals);
}

bool PnPModel::ComputePose(Model *proj_matrix, const std::vector<int> &indices, const double *weights_) const {
    ChooseControlPoints(indices);

    if (!ComputeBarycentricCoordinates(indices)) {
        return false;
    }

    const Eigen::Matrix<double, Eigen::Dynamic, 12> M = ComputeM(indices, weights_);
    const Eigen::Matrix<double, 12, 12> MtM = M.transpose() * M;

    Eigen::JacobiSVD<Eigen::Matrix<double, 12, 12>> svd(
            MtM, Eigen::ComputeFullV | Eigen::ComputeFullU);
    const Eigen::Matrix<double, 12, 12> Ut = svd.matrixU().transpose();

    const Eigen::Matrix<double, 6, 10> L6x10 = ComputeL6x10(Ut);
    const Eigen::Matrix<double, 6, 1> rho = ComputeRho();

    Eigen::Vector4d betas[4];
    std::array<double, 4> reproj_errors;
    std::array<Eigen::Matrix3d, 4> Rs;
    std::array<Eigen::Vector3d, 4> ts;

    FindBetasApprox1(L6x10, rho, &betas[1]);
    RunGaussNewton(L6x10, rho, &betas[1]);
    reproj_errors[1] = ComputeRT(Ut, betas[1], indices, &Rs[1], &ts[1]);

    FindBetasApprox2(L6x10, rho, &betas[2]);
    RunGaussNewton(L6x10, rho, &betas[2]);
    reproj_errors[2] = ComputeRT(Ut, betas[2], indices, &Rs[2], &ts[2]);

    FindBetasApprox3(L6x10, rho, &betas[3]);
    RunGaussNewton(L6x10, rho, &betas[3]);
    reproj_errors[3] = ComputeRT(Ut, betas[3], indices, &Rs[3], &ts[3]);

    int best_idx = 1;
    if (reproj_errors[2] < reproj_errors[1]) {
        best_idx = 2;
    }
    if (reproj_errors[3] < reproj_errors[best_idx]) {
        best_idx = 3;
    }

    *proj_matrix = RTtoMat(Rs[best_idx], ts[best_idx]);
    return true;
}

void PnPModel::ChooseControlPoints(const std::vector<int> &indices) const {
    // Take C0 as the reference points centroid:
    cws_[0].setZero();
    std::vector<int>::const_iterator itIndices = indices.begin();
    for (; itIndices != indices.end(); itIndices++) {
        cws_[0] += _matchesNormalised[*itIndices].p3D;
    }
    cws_[0] /= indices.size();

    Eigen::Matrix<double, Eigen::Dynamic, 3> PW0(indices.size(), 3);
    for (size_t i = 0; i < indices.size(); ++i) {
        PW0.row(i) = _matchesNormalised[indices[i]].p3D - cws_[0];
    }

    const Eigen::Matrix3d PW0tPW0 = PW0.transpose() * PW0;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
            PW0tPW0, Eigen::ComputeFullV | Eigen::ComputeFullU);
    const Eigen::Vector3d D = svd.singularValues();
    const Eigen::Matrix3d Ut = svd.matrixU().transpose();

    for (int i = 1; i < 4; ++i) {
        const double k = std::sqrt(D(i - 1) / indices.size());
        cws_[i] = cws_[0] + k * Ut.row(i - 1).transpose();
    }
}

bool PnPModel::ComputeBarycentricCoordinates(const std::vector<int> &indices) const {
    Eigen::Matrix3d CC;
    for (int i = 0; i < 3; ++i) {
        for (int j = 1; j < 4; ++j) {
            CC(i, j - 1) = cws_[j][i] - cws_[0][i];
        }
    }

    if (CC.colPivHouseholderQr().rank() < 3) {
        return false;
    }

    const Eigen::Matrix3d CC_inv = CC.inverse();

    alphas_.resize(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            alphas_[i][1 + j] = CC_inv(j, 0) * (_matchesNormalised[indices[i]].p3D(0) - cws_[0][0]) +
                                CC_inv(j, 1) * (_matchesNormalised[indices[i]].p3D(1) - cws_[0][1]) +
                                CC_inv(j, 2) * (_matchesNormalised[indices[i]].p3D(2) - cws_[0][2]);
        }
        alphas_[i][0] = 1.0 - alphas_[i][1] - alphas_[i][2] - alphas_[i][3];
    }

    return true;
}

Eigen::Matrix<double, Eigen::Dynamic, 12> PnPModel::ComputeM(const std::vector<int> &indices, const double *weights_) const {
    Eigen::Matrix<double, Eigen::Dynamic, 12> M(2 * indices.size(), 12);
    for (size_t i = 0; i < indices.size(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
            double weight = 1.0;
            if (weights_ != nullptr) {
                weight = weights_[indices[i]];
            }
            M(2 * i, 3 * j) = alphas_[i][j] * weight;
            M(2 * i, 3 * j + 1) = 0.0;
            M(2 * i, 3 * j + 2) = -alphas_[i][j] * weight * _matchesNormalised[indices[i]].p2D.x() * weight;

            M(2 * i + 1, 3 * j) = 0.0;
            M(2 * i + 1, 3 * j + 1) = alphas_[i][j] * weight;
            M(2 * i + 1, 3 * j + 2) = -alphas_[i][j] * _matchesNormalised[indices[i]].p2D.y() * weight;
        }
    }
    return M;
}

Eigen::Matrix<double, 6, 10> PnPModel::ComputeL6x10(
        const Eigen::Matrix<double, 12, 12> &Ut) const {
    Eigen::Matrix<double, 6, 10> L6x10;

    std::array<std::array<Eigen::Vector3d, 6>, 4> dv;
    for (int i = 0; i < 4; ++i) {
        int a = 0, b = 1;
        for (int j = 0; j < 6; ++j) {
            dv[i][j][0] = Ut(11 - i, 3 * a) - Ut(11 - i, 3 * b);
            dv[i][j][1] = Ut(11 - i, 3 * a + 1) - Ut(11 - i, 3 * b + 1);
            dv[i][j][2] = Ut(11 - i, 3 * a + 2) - Ut(11 - i, 3 * b + 2);

            b += 1;
            if (b > 3) {
                a += 1;
                b = a + 1;
            }
        }
    }

    for (int i = 0; i < 6; ++i) {
        L6x10(i, 0) = dv[0][i].transpose() * dv[0][i];
        L6x10(i, 1) = 2.0 * dv[0][i].transpose() * dv[1][i];
        L6x10(i, 2) = dv[1][i].transpose() * dv[1][i];
        L6x10(i, 3) = 2.0 * dv[0][i].transpose() * dv[2][i];
        L6x10(i, 4) = 2.0 * dv[1][i].transpose() * dv[2][i];
        L6x10(i, 5) = dv[2][i].transpose() * dv[2][i];
        L6x10(i, 6) = 2.0 * dv[0][i].transpose() * dv[3][i];
        L6x10(i, 7) = 2.0 * dv[1][i].transpose() * dv[3][i];
        L6x10(i, 8) = 2.0 * dv[2][i].transpose() * dv[3][i];
        L6x10(i, 9) = dv[3][i].transpose() * dv[3][i];
    }

    return L6x10;
}

Eigen::Matrix<double, 6, 1> PnPModel::ComputeRho() const {
    Eigen::Matrix<double, 6, 1> rho;
    rho[0] = (cws_[0] - cws_[1]).squaredNorm();
    rho[1] = (cws_[0] - cws_[2]).squaredNorm();
    rho[2] = (cws_[0] - cws_[3]).squaredNorm();
    rho[3] = (cws_[1] - cws_[2]).squaredNorm();
    rho[4] = (cws_[1] - cws_[3]).squaredNorm();
    rho[5] = (cws_[2] - cws_[3]).squaredNorm();
    return rho;
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_1 = [B11 B12     B13         B14]

void PnPModel::FindBetasApprox1(const Eigen::Matrix<double, 6, 10> &L6x10,
                                const Eigen::Matrix<double, 6, 1> &rho,
                                Eigen::Vector4d *betas) const {
    Eigen::Matrix<double, 6, 4> L_6x4;
    for (int i = 0; i < 6; ++i) {
        L_6x4(i, 0) = L6x10(i, 0);
        L_6x4(i, 1) = L6x10(i, 1);
        L_6x4(i, 2) = L6x10(i, 3);
        L_6x4(i, 3) = L6x10(i, 6);
    }

    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 4>> svd(
            L_6x4, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::Matrix<double, 6, 1> Rho_temp = rho;
    const Eigen::Matrix<double, 4, 1> b4 = svd.solve(Rho_temp);

    if (b4[0] < 0) {
        (*betas)[0] = std::sqrt(-b4[0]);
        (*betas)[1] = -b4[1] / (*betas)[0];
        (*betas)[2] = -b4[2] / (*betas)[0];
        (*betas)[3] = -b4[3] / (*betas)[0];
    } else {
        (*betas)[0] = std::sqrt(b4[0]);
        (*betas)[1] = b4[1] / (*betas)[0];
        (*betas)[2] = b4[2] / (*betas)[0];
        (*betas)[3] = b4[3] / (*betas)[0];
    }
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_2 = [B11 B12 B22                            ]

void PnPModel::FindBetasApprox2(const Eigen::Matrix<double, 6, 10> &L6x10,
                                const Eigen::Matrix<double, 6, 1> &rho,
                                Eigen::Vector4d *betas) const {
    Eigen::Matrix<double, 6, 3> L_6x3(6, 3);

    for (int i = 0; i < 6; ++i) {
        L_6x3(i, 0) = L6x10(i, 0);
        L_6x3(i, 1) = L6x10(i, 1);
        L_6x3(i, 2) = L6x10(i, 2);
    }

    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 3>> svd(
            L_6x3, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::Matrix<double, 6, 1> Rho_temp = rho;
    const Eigen::Matrix<double, 3, 1> b3 = svd.solve(Rho_temp);

    if (b3[0] < 0) {
        (*betas)[0] = std::sqrt(-b3[0]);
        (*betas)[1] = (b3[2] < 0) ? std::sqrt(-b3[2]) : 0.0;
    } else {
        (*betas)[0] = std::sqrt(b3[0]);
        (*betas)[1] = (b3[2] > 0) ? std::sqrt(b3[2]) : 0.0;
    }

    if (b3[1] < 0) {
        (*betas)[0] = -(*betas)[0];
    }

    (*betas)[2] = 0.0;
    (*betas)[3] = 0.0;
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_3 = [B11 B12 B22 B13 B23                    ]

void PnPModel::FindBetasApprox3(const Eigen::Matrix<double, 6, 10> &L6x10,
                                const Eigen::Matrix<double, 6, 1> &rho,
                                Eigen::Vector4d *betas) const {
    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 5>> svd(
            L6x10.leftCols<5>(), Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::Matrix<double, 6, 1> Rho_temp = rho;
    const Eigen::Matrix<double, 5, 1> b5 = svd.solve(Rho_temp);

    if (b5[0] < 0) {
        (*betas)[0] = std::sqrt(-b5[0]);
        (*betas)[1] = (b5[2] < 0) ? std::sqrt(-b5[2]) : 0.0;
    } else {
        (*betas)[0] = std::sqrt(b5[0]);
        (*betas)[1] = (b5[2] > 0) ? std::sqrt(b5[2]) : 0.0;
    }
    if (b5[1] < 0) {
        (*betas)[0] = -(*betas)[0];
    }
    (*betas)[2] = b5[3] / (*betas)[0];
    (*betas)[3] = 0.0;
}

void PnPModel::RunGaussNewton(const Eigen::Matrix<double, 6, 10> &L6x10,
                              const Eigen::Matrix<double, 6, 1> &rho,
                              Eigen::Vector4d *betas) const {
    Eigen::Matrix<double, 6, 4> A;
    Eigen::Matrix<double, 6, 1> b;

    const int kNumIterations = 5;
    for (int k = 0; k < kNumIterations; ++k) {
        for (int i = 0; i < 6; ++i) {
            A(i, 0) = 2 * L6x10(i, 0) * (*betas)[0] + L6x10(i, 1) * (*betas)[1] +
                      L6x10(i, 3) * (*betas)[2] + L6x10(i, 6) * (*betas)[3];
            A(i, 1) = L6x10(i, 1) * (*betas)[0] + 2 * L6x10(i, 2) * (*betas)[1] +
                      L6x10(i, 4) * (*betas)[2] + L6x10(i, 7) * (*betas)[3];
            A(i, 2) = L6x10(i, 3) * (*betas)[0] + L6x10(i, 4) * (*betas)[1] +
                      2 * L6x10(i, 5) * (*betas)[2] + L6x10(i, 8) * (*betas)[3];
            A(i, 3) = L6x10(i, 6) * (*betas)[0] + L6x10(i, 7) * (*betas)[1] +
                      L6x10(i, 8) * (*betas)[2] + 2 * L6x10(i, 9) * (*betas)[3];

            b(i) = rho[i] - (L6x10(i, 0) * (*betas)[0] * (*betas)[0] +
                             L6x10(i, 1) * (*betas)[0] * (*betas)[1] +
                             L6x10(i, 2) * (*betas)[1] * (*betas)[1] +
                             L6x10(i, 3) * (*betas)[0] * (*betas)[2] +
                             L6x10(i, 4) * (*betas)[1] * (*betas)[2] +
                             L6x10(i, 5) * (*betas)[2] * (*betas)[2] +
                             L6x10(i, 6) * (*betas)[0] * (*betas)[3] +
                             L6x10(i, 7) * (*betas)[1] * (*betas)[3] +
                             L6x10(i, 8) * (*betas)[2] * (*betas)[3] +
                             L6x10(i, 9) * (*betas)[3] * (*betas)[3]);
        }

        const Eigen::Vector4d x = A.colPivHouseholderQr().solve(b);

        (*betas) += x;
    }
}

double PnPModel::ComputeRT(const Eigen::Matrix<double, 12, 12> &Ut,
                           const Eigen::Vector4d &betas, const std::vector<int> &indices,
                           Eigen::Matrix3d *R, Eigen::Vector3d *t) const {
    ComputeCcs(betas, Ut);
    ComputePcs(indices.size());

    SolveForSign(indices.size());

    EstimateRT(R, t, indices);

    return ComputeTotalReprojectionError(*R, *t);
}

void PnPModel::ComputeCcs(const Eigen::Vector4d &betas,
                          const Eigen::Matrix<double, 12, 12> &Ut) const {
    for (int i = 0; i < 4; ++i) {
        ccs_[i][0] = ccs_[i][1] = ccs_[i][2] = 0.0;
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 3; ++k) {
                ccs_[j][k] += betas[i] * Ut(11 - i, 3 * j + k);
            }
        }
    }
}

void PnPModel::ComputePcs(const size_t sampleSize) const {
    pcs_.resize(sampleSize);
    for (size_t i = 0; i < sampleSize; ++i) {
        for (int j = 0; j < 3; ++j) {
            pcs_[i][j] = alphas_[i][0] * ccs_[0][j] + alphas_[i][1] * ccs_[1][j] +
                         alphas_[i][2] * ccs_[2][j] + alphas_[i][3] * ccs_[3][j];
        }
    }
}

void PnPModel::SolveForSign(const size_t sampleSize) const {
    if (pcs_[0][2] < 0.0) {
        for (int i = 0; i < 4; ++i) {
            ccs_[i] = -ccs_[i];
        }
        for (size_t i = 0; i < sampleSize; ++i) {
            pcs_[i] = -pcs_[i];
        }
    }
}

void PnPModel::EstimateRT(Eigen::Matrix3d *R, Eigen::Vector3d *t, const std::vector<int> &indices) const {
    Eigen::Vector3d pc0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d pw0 = Eigen::Vector3d::Zero();

    for (size_t i = 0; i < indices.size(); ++i) {
        pc0 += pcs_[i];
        double weight = 1.0;
        pw0 += _matchesNormalised[indices[i]].p3D;
    }
    pc0 /= indices.size();
    pw0 /= indices.size();

    Eigen::Matrix3d abt = Eigen::Matrix3d::Zero();
    for (size_t i = 0; i < indices.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            abt(j, 0) += (pcs_[i][j] - pc0[j]) * (_matchesNormalised[indices[i]].p3D(0) - pw0[0]);
            abt(j, 1) += (pcs_[i][j] - pc0[j]) * (_matchesNormalised[indices[i]].p3D(1) - pw0[1]);
            abt(j, 2) += (pcs_[i][j] - pc0[j]) * (_matchesNormalised[indices[i]].p3D(2) - pw0[2]);
        }
    }

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
            abt, Eigen::ComputeFullV | Eigen::ComputeFullU);
    const Eigen::Matrix3d abt_U = svd.matrixU();
    const Eigen::Matrix3d abt_V = svd.matrixV();

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            (*R)(i, j) = abt_U.row(i) * abt_V.row(j).transpose();
        }
    }

    if (R->determinant() < 0) {
        Eigen::Matrix3d Abt_v_prime = abt_V;
        Abt_v_prime.col(2) = -abt_V.col(2);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                (*R)(i, j) = abt_U.row(i) * Abt_v_prime.row(j).transpose();
            }
        }
    }

    *t = pc0 - *R * pw0;
}

double PnPModel::ComputeTotalReprojectionError(const Eigen::Matrix3d &R,
                                               const Eigen::Vector3d &t) const {
    Model proj_matrix(3, 4);
    proj_matrix = RTtoMat(R, t);

    std::vector<double> residuals;
    ComputeSquaredReprojectionError(proj_matrix, &residuals);

    double reproj_error = 0.0;
    for (const double residual : residuals) {
        reproj_error += std::sqrt(residual);
    }

    return reproj_error;
}

/// Computes the area of the inlier region according to the AC-RANSAC
/// paper for a given error margin. The inlier region is over-estimated to
/// ease computation. This is a point-to-point error. Used for LRT.
/// \param sigma The error margin.
double PnPModel::pSigma(const double sigma, const bool leftSide) const {
    return M_PI * sigma * sigma / (double) _area;
}

} // namespace orsa
