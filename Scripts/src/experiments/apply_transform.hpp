//
// Created by clementriu on 8/17/20.
//

#ifndef MMM_ORSA_APPLY_TRANSFORM_HPP
#define MMM_ORSA_APPLY_TRANSFORM_HPP

#include <random>

#include "libImage/image_io.hpp"

#include "libOrsa/match.hpp"
#include "libOrsa/match2d3d.hpp"
#include "libOrsa/model_estimator.hpp"
#include "libOrsa/libNumerics/homography.h"

#include "utilities/Rect.hpp"


/// Apply an homography to an image and to the set of points the homography was computed with.
/// \param[in] image1: The image on which to apply the homography.
/// \param[in] image2: The true second image from which the homography was computed.
/// \param[in] H: The homography.
/// \param[in] rect: The intersection of the two images.
/// \param[in] vec_inliers: Vector of the indexes of the matches from which the homography was computed.
/// \param[in] matchings: Vector of all the matches.
/// \param[out] newMatchings: Vector of the new matches.
/// \param[in] noise: Type of noise: -1 noise is applied before the homography, 1 noise is applied after the homography,
/// otherwise no noise is applied. (Default: 0)
/// \param[in] std: Standard deviation of the noise if some is applied. (Default: 0)
/// \return The warped image.
Image<RGBColor> apply_homography(const Image<RGBColor> &image1,
                                 const Image<RGBColor> &image2,
                                 const libNumerics::matrix<double> &H,
                                 const utility::Rect &rect,
                                 const std::vector<int> &vec_inliers,
                                 const std::vector<Match> &matchings,
                                 std::vector<Match> &newMatchings,
                                 std::default_random_engine &generator,
                                 int noise = 0,
                                 double std = 0);

/// Transform a set of matches to perfectly match a given fundamenal matrix.
/// Works by projecting a inlier-match's second-image-point on its associated epipolar line.
/// \param[in] image2: second image, on which width and height will be extracted.
/// \param[in] F: Fundamental matrix to fit.
/// \param[in] vec_inliers: Vectors of index of matches to transform in matchings.
/// \param[in] matchings: Vectors of matches to transform.
/// \param[out] newMatchings: Vectors of transformed matches.
void apply_fundamental(const Image<RGBColor> &image2,
                       const libNumerics::matrix<double> &F,
                       const std::vector<int> &vec_inliers,
                       const std::vector<Match> &matchings,
                       std::vector<Match> &newMatchings);

void applyPnP(const Image<RGBColor> &image,
              const libNumerics::matrix<double> &RT,
              const std::vector<int> &vec_inliers,
              const std::vector<Match2D3D> &matchingsNormalised,
              const libNumerics::matrix<double> &calib,
              std::vector<Match2D3D> &newMatchings);

/// Generate random matches in two images. Those matches can be chosen to be outliers of a given model.
/// \param[in] w1: Width of the first image.
/// \param[in] h1: Height of the first image.
/// \param[in] w2: Width of the second image.
/// \param[in] h2: Height of the second image.
/// \param[in] outlierRatio: Ratio of outliers to create.
/// If it is a value between 0 and 1 excluded it will be a ratio of the inlierCount parameter.
/// If a value greater than 1 included it will be the number of created outliers.
/// \param[in] inlierCount: Number of inliers. Used when the outlierRatio is below 1.
/// \param[in] addOutlier: Type of outliers: if addOutlier = 2 the matches will not correspond to inliers of the given model.
/// \param[in] model: Used model when addOutlier = 2.
/// \param[in] params: Parameters of the model when addOutlier = 2.
/// \param[in] outlierThresholdSq: Squared inlier/outlier threshold when addOutlier = 2.
/// \param[out] artificialOutliers: Vector to store the generated matches.
/// \param[in] maxIterOutlier: Number of rejects authorised when addOutlier = 2. (Default: 1000)
/// \return True if as many ouliers are created as wanted, false otherwise.
bool generateOutliers(int w1, int h1, int w2, int h2, double outlierRatio,
                      int inlierCount, int addOutlier, const orsa::ModelEstimator *model,
                      const libNumerics::matrix<double> &H, double outlierThresholdSq,
                      std::vector<Match> &artificialOutliers, std::default_random_engine &generator, int maxIterOutlier = 1000);

/// Generate random matches in two images. Those matches can be chosen to be outliers of a given model. Those matches are generated to ensure a more realistic distribution of outliers in the error space. Only works for point error space.
/// \param[in] w1: Width of the first image.
/// \param[in] h1: Height of the first image.
/// \param[in] w2: Width of the second image.
/// \param[in] h2: Height of the second image.
/// \param[in] outlierRatio: Ratio of outliers to create.
/// If it is a value between 0 and 1 excluded it will be a ratio of the inlierCount parameter.
/// If a value greater than 1 included it will be the number of created outliers.
/// \param[in] inlierCount: Number of inliers. Used when the outlierRatio is below 1.
/// \param[in] addOutlier: Type of outliers: if addOutlier = 2 the matches will not correspond to inliers of the given model.
/// \param[in] model: Used model when addOutlier = 2.
/// \param[in] params: Parameters of the model when addOutlier = 2.
/// \param[in] outlierThresholdSq: Squared inlier/outlier threshold when addOutlier = 2.
/// \param[out] artificialOutliers: Vector to store the generated matches.
/// \param[in] maxIterOutlier: Number of rejects authorised when addOutlier = 2. (Default: 1000)
/// \return True if as many ouliers are created as wanted, false otherwise.
bool generateHomOutliersUniformError(int w1, int h1, int w2, int h2,
                                  double outlierRatio, int inlierCount,
                                  int addOutlier,
                                  const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                                  double outlierThresholdSq,
                                  std::vector<Match> &artificialOutliers,
                                  std::default_random_engine &generator,
                                  int maxIterOutlier);

/// Generate random matches in two images. Those matches can be chosen to be outliers of a given model. Those matches are generated to ensure a more realistic distribution of outliers in the error space. Only works for line error space.
/// \param[in] w1: Width of the first image.
/// \param[in] h1: Height of the first image.
/// \param[in] w2: Width of the second image.
/// \param[in] h2: Height of the second image.
/// \param[in] outlierRatio: Ratio of outliers to create.
/// If it is a value between 0 and 1 excluded it will be a ratio of the inlierCount parameter.
/// If a value greater than 1 included it will be the number of created outliers.
/// \param[in] inlierCount: Number of inliers. Used when the outlierRatio is below 1.
/// \param[in] addOutlier: Type of outliers: if addOutlier = 2 the matches will not correspond to inliers of the given model.
/// \param[in] model: Used model when addOutlier = 2.
/// \param[in] params: Parameters of the model when addOutlier = 2.
/// \param[in] outlierThresholdSq: Squared inlier/outlier threshold when addOutlier = 2.
/// \param[out] artificialOutliers: Vector to store the generated matches.
/// \param[in] maxIterOutlier: Number of rejects authorised when addOutlier = 2. (Default: 1000)
/// \return True if as many ouliers are created as wanted, false otherwise.
bool generateFundOutlierUniform(int w1, int h1, int w2, int h2,
                                double outlierRatio, int inlierCount,
                                int addOutlier,
                                const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                                double outlierThresholdSq,
                                std::vector<Match> &artificialOutliers,
                                std::default_random_engine &generator,
                                int maxIterOutlier);

bool generatePNPOutlierUniform(double xmin, double xmax,
                               double ymin, double ymax,
                               double zmin, double zmax,
                               int w, int h,
                               double outlierRatio, int inlierCount,
                               int addOutlier,
                               const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                               const libNumerics::matrix<double> &intrinsics,
                               double outlierThresholdSq,
                               std::vector<Match2D3D> &artificialOutliers,
                               std::default_random_engine &generator,
                               int maxIterOutlier);

/// Add noise to matches of a vector while staying in a specified range.
/// \param[in] vectToModify: Vector on which to apply the noising.
/// \param[in] width: Maximum value of x.
/// \param[in] height: Maximum value of y.
/// \param[out] vectModified: Vector of modified matches.
/// \param[in] std: Standard deviation of the gaussian noise or range of the uniform noise. (Default: 0)
/// \param[in] gaussian: If true gaussian noise is used, uniform noise othewise. (Default: false)
void safeAddNoise(const std::vector<Match> &vectToModify,
                  int width, int height,
                  std::vector<Match> &vectModified,
                  std::default_random_engine &generator,
                  double std = 0, bool gaussian = false,
                  int numInliers = 0);

/// Add noise to matches of a vector while staying in a specified range.
/// \param[in] vectToModify: Vector on which to apply the noising.
/// \param[in] width: Maximum value of x.
/// \param[in] height: Maximum value of y.
/// \param[out] vectModified: Vector of modified matches.
/// \param[in] std: Standard deviation of the gaussian noise or range of the uniform noise. (Default: 0)
/// \param[in] gaussian: If true gaussian noise is used, uniform noise otherwise. (Default: false)
void safeAddNoise(const std::vector<Match2D3D> &vectToModify,
                  int width, int height,
                  std::vector<Match2D3D> &vectModified,
                  std::default_random_engine &generator,
                  double std, bool gaussian,
                  const int numInliers);

///// Find the maximum error of a model.
///// \param[in] vectToEvaluate: The vector of points to evaluate.
///// \param[in] model: The model the data fits.
///// \param[in] modelParams: The model parameters to use for computation.
///// \return The squared maximum error found.
//double findMaxErrorSq(const std::vector<Match> &vectToEvaluate,
//                      const orsa::ModelEstimator *model, const libNumerics::matrix<double> &modelParams);

/// Find the maximum error of a model.
/// \param[in] vectToEvaluate: The vector of points to evaluate.
/// \param[in] model: The model the data fits.
/// \param[in] modelParams: The model parameters to use for computation.
/// \return The squared maximum error found.
template <typename M>
double findMaxErrorSq(const std::vector<M> &vectToEvaluate,
                      const orsa::ModelEstimator *model, const libNumerics::matrix<double> &modelParams) {
    double errorMax = 0;
    typename std::vector<M>::const_iterator itVectToEvaluate = vectToEvaluate.begin();
    for (; itVectToEvaluate != vectToEvaluate.end(); ++itVectToEvaluate) {
        typename std::vector<M> outlierVect{*itVectToEvaluate};
        double error = model->Error(modelParams, M::toMat(outlierVect));
        if (error > errorMax) {
            errorMax = error;
        }
    }
    return errorMax;
}

/// Find the maximum error of a model.
/// \param[in] vectIdxToEvaluate: The vector of index of points to evaluate in model.
/// \param[in] model: The model the data fits.
/// \param[in] modelParams: The model parameters to use for computation.
/// \return The squared maximum error found.
double findMaxErrorSq(const std::vector<int> &vectIdxToEvaluate,
                      const orsa::ModelEstimator *model, const libNumerics::matrix<double> &modelParams);

void findOutlier2D3DRange(const std::vector<Match2D3D> &matches,
                          double &xMin, double &xMax,
                          double &yMin, double &yMax,
                          double &zMin, double &zMax,
                          double increaseRatio);

#endif //MMM_ORSA_APPLY_TRANSFORM_HPP
