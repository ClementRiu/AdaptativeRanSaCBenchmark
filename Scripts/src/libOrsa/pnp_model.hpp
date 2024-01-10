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

#ifndef MMM_ORSA_PNP_MODEL_HPP
#define MMM_ORSA_PNP_MODEL_HPP

#include "eigen/Eigen/Eigen"
#include <vector>

#include "model_estimator.hpp"
#include "match2d3d.hpp"

namespace orsa {
/// Apply calibration matrix normalisation to the points to go from pixel points to camera points.
void normalisePoints(const std::vector<Match2D3D> &points_,
                     const libNumerics::matrix<double> &calibration,
                     std::vector<Match2D3D> &normalisedPoints_);

void convert3Dto2D(const ModelEstimator::Model &proj_matrix,
                   const libNumerics::matrix<double> &intrinsic,
                   Match2D3D &m);

class PnPModel : public ModelEstimator {
public:
    PnPModel(const std::vector<Match2D3D> &m, const std::vector<Match2D3D> &mNormalised,
             const libNumerics::matrix<double> &intrinsics_,
             int width, int height);
    virtual ~PnPModel() {} // Necessary for safe subclass

    /// 3 points are required to compute the model.
    int SizeSample() const { return 3; }

    /// 1 model is computed from a sample of kMinNumSamples points.
    int NbModels() const { return 4; }

    /// Degree of freedom of the model.
    int nDegreeOfFreedom() const { return 6; };
    /// Degree of freedom of the data (for magsac++)
    int dataDegreeOfFreedom() const {return 5;};
    /// C(n) = 1/(n/2 * Gamma(n/2)) where n is DoF (for magsac++)
    double C() const {return 0.13298076;};
    /// UG((DoF - 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
    double upperGamma() const {return 0.0101184589199229014;};
    /// LG((DoF + 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
    double lowerGamma() const {return 1.9215217286137212;};

    /// Distance used to distinguish inlier/outlier is to a line.
    virtual bool DistToPoint() const { return true; }

    // Estimate the most probable solution of the P3P problem from a set of
    // kMinNumSamples 2D-3D point correspondences.
    void Fit(const std::vector<int> &indices, std::vector<Model> *RTs) const;

    void Fit(const std::vector<int> &indices,
             std::vector<Model> *RTs,
             const double *weights_) const;

    double Error(const Model &RT, int index, int *side = 0) const;

    double Error(const Model &RT, Mat testData, int *side = 0) const;

    /// Computes the area of the inlier region according to the AC-RANSAC
    /// paper for a given error margin. The inlier region is over-estimated to
    /// ease computation. This is a point-to-line error. Used for LRT.
    /// \param sigma The error margin.
    double pSigma(double sigma, bool leftSide = true) const;

    /// Returns the fundamental matrix. For Fundamental matrix the output is the input.
    /// Necessary for essential matrix.
    inline Model toPixelSpace(const Model &RT) const { return RT; }

private:
    void EPNP(const std::vector<int> &indices, std::vector<Model> *RTs, const double *weights_ = nullptr) const;
    void P3P(const std::vector<int> &indices, std::vector<Model> *RTs) const;

    /// Square reprojection error for a given point given R and T.
    double Error(const Model &RT,
                 double xa, double xb,
                 double ya, double yb,
                 double yc, int *side = 0) const;

    // Calculate the squared reprojection error given a set of 2D-3D point
    // correspondences and a projection matrix.
    //
    // @param points2D     Normalized 2D image points as Nx2 matrix.
    // @param points3D     3D world points as Nx3 matrix.
    // @param proj_matrix  3x4 projection matrix.
    // @param residuals    Output vector of residuals.
    void Residuals(const Model &proj_matrix, std::vector<double> *residuals) const;

    static Mat RTtoMat(const Eigen::Matrix3d &R, const Eigen::Vector3d &T);

    void ComputeSquaredReprojectionError(
            const Model &proj_matrix, std::vector<double> *residuals) const;

    bool ComputePose(Model *proj_matrix, const std::vector<int> &indices, const double *weights_=nullptr) const;

    void ChooseControlPoints(const std::vector<int> &indices) const;

    bool ComputeBarycentricCoordinates(const std::vector<int> &indices) const;

    Eigen::Matrix<double, Eigen::Dynamic, 12> ComputeM(const std::vector<int> &indices, const double *weights_ = nullptr) const;

    Eigen::Matrix<double, 6, 10> ComputeL6x10(
            const Eigen::Matrix<double, 12, 12> &Ut) const;

    Eigen::Matrix<double, 6, 1> ComputeRho() const;

    void FindBetasApprox1(const Eigen::Matrix<double, 6, 10> &L_6x10,
                          const Eigen::Matrix<double, 6, 1> &rho,
                          Eigen::Vector4d *betas) const;

    void FindBetasApprox2(const Eigen::Matrix<double, 6, 10> &L_6x10,
                          const Eigen::Matrix<double, 6, 1> &rho,
                          Eigen::Vector4d *betas) const;

    void FindBetasApprox3(const Eigen::Matrix<double, 6, 10> &L_6x10,
                          const Eigen::Matrix<double, 6, 1> &rho,
                          Eigen::Vector4d *betas) const;

    void RunGaussNewton(const Eigen::Matrix<double, 6, 10> &L_6x10,
                        const Eigen::Matrix<double, 6, 1> &rho,
                        Eigen::Vector4d *betas) const;

    double ComputeRT(const Eigen::Matrix<double, 12, 12> &Ut,
                     const Eigen::Vector4d &betas, const std::vector<int> &indices,
                     Eigen::Matrix3d *R, Eigen::Vector3d *t) const;

    void ComputeCcs(const Eigen::Vector4d &betas,
                    const Eigen::Matrix<double, 12, 12> &Ut) const;

    void ComputePcs(size_t sampleSize) const;

    void SolveForSign(const size_t sampleSize) const;

    void EstimateRT(Eigen::Matrix3d *R, Eigen::Vector3d *t, const std::vector<int> &indices) const;

    double ComputeTotalReprojectionError(const Eigen::Matrix3d &R,
                                         const Eigen::Vector3d &t) const;

    const int _area;
    const std::vector<Match2D3D> _matches;
    const std::vector<Match2D3D> _matchesNormalised;
    libNumerics::matrix<double> _intrinsic, _invIntrinsic;
    mutable std::vector<Eigen::Vector3d> pcs_;
    mutable std::vector<Eigen::Vector4d> alphas_;
    mutable std::array<Eigen::Vector3d, 4> cws_;
    mutable std::array<Eigen::Vector3d, 4> ccs_;
};
} // namespace orsa

#endif //MMM_ORSA_PNP_MODEL_HPP
