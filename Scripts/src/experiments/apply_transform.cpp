//
// Created by clementriu on 8/17/20.
//

#include "apply_transform.hpp"

#include <cmath>

#include "libOrsa/pnp_model.hpp"
#include "libOrsa/sampling.hpp"

#include "utilities/fundamental_graphical_output.hpp"
#include "utilities/homography_graphical_output.hpp"
#include "utilities/warping.hpp"


/// Truncate a value to a range.
/// \param[in] min: Lower bound of the range.
/// \param[in] max: Upper bound of range.
/// \param[in, out] value: Value that will be updated to be included in the desired range.
void truncate(const double min, const double max, double &value) {
    assert(min <= max);
    if (value < min) {
        value = min;
    }
    if (value > max) {
        value = max;
    }
}

/// Add gaussian noise of mean 0 and specified std to a point.
/// \param[in] std: Standard deviation of the gaussian noise.
/// \param[in, out] x: First coordinate of the point. Will be updated with the noise.
/// \param[in, out] y: Second coordinate of the point. Will be updated with the noise.
/// \param[in] generator Random value generator.
void addNoise(const double std, double &x, double &y, std::default_random_engine &generator) {
    std::normal_distribution<double> distributionGaussian(0.0, std);

    double xNoise, yNoise;
    xNoise = distributionGaussian(generator);
    yNoise = distributionGaussian(generator);

    x += xNoise;
    y += yNoise;
}


/// Add gaussian noise of mean 0 and specified std to a point bounded to NOISEBOUND * std.
/// Or add uniform noise in [-std, +std].
/// Resulting values are cropped to a range [0, width] x [0, height].
/// \param[in] std: Standard deviation of the gaussian noise or range of the uniform noise.
/// \param[in, out] x: First coordinate of the point. Will be updated with the noise.
/// \param[in, out] y: Second coordinate of the point. Will be updated with the noise.
/// \param[in] generator: Random value generator.
/// \param[in] width: Maximum value of x.
/// \param[in] height: Maximum value of y.
/// \param[in] gaussian: If true gaussian noise is used, uniform noise othewise. (Default: false)
void addNoiseBounded(const double std,
                     double &x, double &y,
                     std::default_random_engine &generator,
                     const double width, const double height,
                     const bool gaussian = false) {
    double xNoise, yNoise;

    if (gaussian) {
        const double NOISEBOUND = 3.0;
        std::normal_distribution<double> distributionGaussian(0.0, std);
        xNoise = distributionGaussian(generator);
        yNoise = distributionGaussian(generator);
        truncate(-NOISEBOUND * std, NOISEBOUND * std, xNoise);
        truncate(-NOISEBOUND * std, NOISEBOUND * std, yNoise);
    } else {
        std::uniform_real_distribution<double> distributionUniform(-std, std);
        xNoise = distributionUniform(generator);
        yNoise = distributionUniform(generator);
    }

    x += xNoise;
    y += yNoise;

    truncate(0, width, x);
    truncate(0, height, y);
}

/// Apply an homography to a list of inlier points and leaves possibility to add noise to those points.
/// Will create a list of matches corresponding to the desired homography.
/// \param[in] vec_inliers: Vector of indexes of the matches on which to apply the homography.
/// \param[in] matchings: Vector of matches from which to get the points.
/// \param[in] H: Homography to apply.
/// \param[out] newMatchings: Vector of new matches.
/// \param[in] w: Width of the image.
/// \param[in] h: Height of the image.
/// \param[in] noise: Type of noise: -1 noise is applied before the homography, 1 noise is applied after the homography,
/// otherwise no noise is applied. (Default: 0)
/// \param[in] std: Standard deviation of the noise if some is applied. (Default: 0)
void transformPointsWithNoise(const std::vector<int> &vec_inliers,
                              const std::vector<Match> &matchings,
                              const libNumerics::matrix<double> &H,
                              std::vector<Match> &newMatchings,
                              const int w,
                              const int h,
                              std::default_random_engine &generator,
                              const int noise = 0,
                              const double std = 0) {


    std::vector<int>::const_iterator it = vec_inliers.begin();
    for (; it != vec_inliers.end(); ++it) {
        double x2 = static_cast<double>(matchings[*it].x1);
        double y2 = static_cast<double>(matchings[*it].y1);
        if (noise == -1) {
            addNoise(std, x2, y2, generator);
        }
        if (TransformH(H, x2, y2)) { // Transformation of the points.
            if (noise == 1) {
                addNoise(std, x2, y2, generator);
            }
            truncate(0, w, x2);
            truncate(0, h, y2);
            Match newMatch = matchings[*it];
            newMatch.x2 = static_cast<float>(x2);
            newMatch.y2 = static_cast<float>(y2);
            newMatchings.push_back(newMatch);
        }
    }
}

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
                                 const int noise,
                                 const double std) {
    int xc = (rect.left + rect.right) / 2;
    int yc = (rect.top + rect.bottom) / 2;
    size_t wM = std::max(image1.Width(), image2.Width());
    size_t hM = std::max(image1.Height(), image2.Height());
    int xo = static_cast<int>(wM / 2);
    int yo = static_cast<int>(hM / 2);
    libNumerics::matrix<double> T = utility::translation(xo - xc, yo - yc);

    // Modification of the image:
    Image<RGBColor> warpedImage(wM, hM, WHITE);
    utility::Warp(image1, T * H, warpedImage);

    // Mofidication of the points:
    transformPointsWithNoise(vec_inliers, matchings, H, newMatchings, wM, hM, generator, noise, std);
    return warpedImage;
}

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
                       std::vector<Match> &newMatchings) {

    std::vector<int>::const_iterator itInliers = vec_inliers.begin();

    for (; itInliers != vec_inliers.end(); ++itInliers) {
        Match matchToModify = matchings[*itInliers];
        libNumerics::vector<double> epi = utility::epipolar_line(F.t(), matchToModify);
        if (epi.qnorm() == 0) continue;
        double xp = matchToModify.x2, yp = matchToModify.y2;
        utility::project_on_epi(xp, yp, epi);

        double width = image2.Width();
        double height = image2.Height();
        if (xp >= 0 && xp <= width && yp >= 0 && yp <= height) {
            Match newMatch = matchToModify;
            newMatch.x2 = static_cast<float>(xp);
            newMatch.y2 = static_cast<float>(yp);
            newMatchings.push_back(newMatch);
        }
    }
}

/// Transform a set of matches to perfectly match a given fundamenal matrix. TODO
/// Works by projecting a inlier-match's second-image-point on its associated epipolar line.
/// \param[in] image2: second image, on which width and height will be extracted.
/// \param[in] F: Fundamental matrix to fit.
/// \param[in] vec_inliers: Vectors of index of matches to transform in matchings.
/// \param[in] matchings: Vectors of matches to transform.
/// \param[out] newMatchings: Vectors of transformed matches.
void applyPnP(const Image<RGBColor> &image,
              const libNumerics::matrix<double> &RT,
              const std::vector<int> &vec_inliers,
              const std::vector<Match2D3D> &matchingsNormalised,
              const libNumerics::matrix<double> &calib,
              std::vector<Match2D3D> &newMatchings) {

    std::vector<int>::const_iterator itInliers = vec_inliers.begin();

    for (; itInliers != vec_inliers.end(); ++itInliers) {
        Match2D3D matchToModify = matchingsNormalised[*itInliers];
        double x = matchToModify.p3D(0), y = matchToModify.p3D(1), z = matchToModify.p3D(2);
        Match2D3D tempMatch(0, 0, static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
        orsa::convert3Dto2D(RT, calib, tempMatch);

        double width = image.Width();
        double height = image.Height();
        if (tempMatch.p2D(0) >= 0 && tempMatch.p2D(0) <= width
            && tempMatch.p2D(1) >= 0 && tempMatch.p2D(1) <= height) {
            newMatchings.push_back(tempMatch);
        }
    }
}

/// Generate random uniform matches in two images. Those matches can be chosen to be outliers of a given model.
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
bool generateOutliers(const int w1, const int h1, const int w2, const int h2,
                      const double outlierRatio, const int inlierCount,
                      const int addOutlier,
                      const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                      const double outlierThresholdSq,
                      std::vector<Match> &artificialOutliers,
                      std::default_random_engine &generator,
                      const int maxIterOutlier) {
    int iter = 0;
    int numOutlierCreated = 0;
    int numOutlierWanted;

    std::uniform_real_distribution<double> distributionUniformWidth1(0, w1);
    std::uniform_real_distribution<double> distributionUniformHeight1(0, h1);
    std::uniform_real_distribution<double> distributionUniformWidth2(0, w2);
    std::uniform_real_distribution<double> distributionUniformHeight2(0, h2);

    // Selection of the number of outlier to generate:
    if (outlierRatio < 1.0) {
        numOutlierWanted = static_cast<int>(std::floor(outlierRatio * inlierCount / (1 - outlierRatio)));
    } else {
        numOutlierWanted = static_cast<int>(outlierRatio);
    }

    while ((numOutlierCreated < numOutlierWanted) && (iter < numOutlierWanted + maxIterOutlier)) {
        double x1, y1, x2, y2;
        x1 = distributionUniformWidth1(generator);
        y1 = distributionUniformHeight1(generator);
        x2 = distributionUniformWidth2(generator);
        y2 = distributionUniformHeight2(generator);

        Match outlier(static_cast<float>(x1),
                      static_cast<float>(y1),
                      static_cast<float>(x2),
                      static_cast<float>(y2));

        std::vector<Match> outlierVect{outlier};

        // The point is kept if addOutlier != 2 or if addOutlier = 2 and the point is not an inlier:
        if (!((addOutlier == 2) &&
              (model->Error(params, Match::toMat(outlierVect)) <= outlierThresholdSq))) {
            numOutlierCreated += 1;

            artificialOutliers.push_back(outlier);
        }
        iter += 1;
    }
    return (numOutlierCreated == numOutlierWanted);
}

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
bool generateHomOutliersUniformError(const int w1, const int h1, const int w2, const int h2,
                                     const double outlierRatio, const int inlierCount,
                                     const int addOutlier,
                                     const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                                     const double outlierThresholdSq,
                                     std::vector<Match> &artificialOutliers,
                                     std::default_random_engine &generator,
                                     const int maxIterOutlier) {
    int iter = 0;
    int numOutlierCreated = 0;
    int numOutlierWanted;

    // Definition of the basic generators: for the first point in the left image and the orientation.
    std::uniform_real_distribution<double> distributionUniformWidth1(0, w1);
    std::uniform_real_distribution<double> distributionUniformHeight1(0, h1);
    std::uniform_real_distribution<long double> distributionUniformOrientation(0, 2 * M_PI);

    // Selection of the number of outlier to generate:
    if (outlierRatio < 1.0) {
        numOutlierWanted = static_cast<int>(std::floor(outlierRatio * inlierCount / (1 - outlierRatio)));
    } else {
        numOutlierWanted = static_cast<int>(outlierRatio);
    }

    while ((numOutlierCreated < numOutlierWanted) && (iter < numOutlierWanted + maxIterOutlier)) {

        // Generation of a point in the left image and its perfect match.
        double x1, y1, x2, y2;
        x1 = distributionUniformWidth1(generator);
        y1 = distributionUniformHeight1(generator);
        x2 = x1;
        y2 = y1;
        TransformH(params, x2, y2);
        if (x2 >= 0 && x2 <= w2 && y2 >= 0 && y2 <= h2) {

            double minOffSet, maxOffSet;

            // If addOutlier == 2 : we define a minimum distance between the perfect match and the outlier match.
            if (addOutlier == 2) {
                minOffSet = std::sqrt(outlierThresholdSq);
            } else {
                minOffSet = 0;
            }

            // To generate a second point in the image, the space is cut in 4 parts, to determine the available
            // distances between the match and the image border.
            double theta1, theta2, theta3, theta4;
            theta1 = std::atan((h2 - y2) / (double) (w2 - x2));
            theta2 = std::atan(x2 / (double) (h2 - y2)) + M_PI / (double) 2;
            theta3 = std::atan((y2) / (double) (x2)) + M_PI;
            theta4 = std::atan((w2 - x2) / (double) (y2)) + M_PI * 3 / (double) 2;

            assert(theta1 > 0 && theta1 < M_PI / (double) 2);
            assert(theta2 > M_PI / (double) 2 && theta2 < M_PI);
            assert(theta3 > M_PI && theta3 < M_PI * 3 / (double) 2);
            assert(theta4 > M_PI * 3 / (double) 2 && theta4 < 2 * M_PI);

            // A direction is sampled and the distance between the match and the image border is computed.
            double theta = distributionUniformOrientation(generator);

            if (theta >= theta4 || theta < theta1) {
                maxOffSet = (double) (w2 - x2) / std::cos(theta);
            }
            if (theta >= theta1 && theta < theta2) {
                maxOffSet = (double) (h2 - y2) / std::cos(theta - (M_PI / (double) 2));
            }
            if (theta >= theta2 && theta < theta3) {
                maxOffSet = (double) (x2) / std::cos(theta - M_PI);
            }
            if (theta >= theta3 && theta < theta4) {
                maxOffSet = (double) (y2) / std::cos(theta - (M_PI * 3 / (double) 2));
            }

            if (maxOffSet >= minOffSet) {
                // A point is drawn along the drawn direction.
                std::uniform_real_distribution<double> distributionUniformOutlierError(minOffSet, maxOffSet);
                double outlierError = distributionUniformOutlierError(generator);

                x2 += std::cos(theta) * outlierError;
                y2 += std::sin(theta) * outlierError;

                assert(x2 >= 0 && x2 <= w2);
                assert(y2 >= 0 && y2 <= h2);

                Match outlier(static_cast<float>(x1),
                              static_cast<float>(y1),
                              static_cast<float>(x2),
                              static_cast<float>(y2));

                std::vector<Match> outlierVect{outlier};

                // The point is kept if addOutlier != 2 or if addOutlier = 2 and the point is not an inlier:
                assert(!((addOutlier == 2) &&
                         (model->Error(params, Match::toMat(outlierVect)) <= outlierThresholdSq)));
                numOutlierCreated += 1;

                artificialOutliers.push_back(outlier);
            }
        }
        iter += 1;
    }

    return (numOutlierCreated == numOutlierWanted);
}

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
bool generateFundOutlierUniform(const int w1, const int h1, const int w2, const int h2,
                                const double outlierRatio, const int inlierCount,
                                const int addOutlier,
                                const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                                const double outlierThresholdSq,
                                std::vector<Match> &artificialOutliers,
                                std::default_random_engine &generator,
                                const int maxIterOutlier) {
    int iter = 0;
    int numOutlierCreated = 0;
    int numOutlierWanted;

    // Definition of the basic generators: for the first point in the left image and the orientation (above or under the
    // epipolar line).
    std::uniform_real_distribution<double> distributionUniformWidth1(0, w1);
    std::uniform_real_distribution<double> distributionUniformHeight1(0, h1);
    std::uniform_int_distribution<int> distributionUniform01(0, 1);

    // Selection of the number of outlier to generate:
    if (outlierRatio < 1.0) {
        numOutlierWanted = static_cast<int>(std::floor(outlierRatio * inlierCount / (1 - outlierRatio)));
    } else {
        numOutlierWanted = static_cast<int>(outlierRatio);
    }

    while ((numOutlierCreated < numOutlierWanted) && (iter < numOutlierWanted + maxIterOutlier)) {

        // Generation of a point in the left image and its epipolar line.
        double x1, y1, x2, y2;
        x1 = distributionUniformWidth1(generator);
        y1 = distributionUniformHeight1(generator);
        Match tempMatch(static_cast<float>(x1), static_cast<float>(y1), 0, 0);
        libNumerics::vector<double> epi = utility::epipolar_line(params.t(), tempMatch);

        double a = epi(0), b = epi(1), c = epi(2);

        // Computation of the possible values of x to find a point on the epipolar line.
        iter += 1;
        if (std::abs(a) > 1.0e-10 &&
            std::abs(b) > 1.0e-10) { // TODO: use global parameter for definition of "too small".
            // Definition of the acceptable range of x to chose a point on the epipolar line.
            double xMin = 0.0;
            double xMax = (double) w2;
            double slope = -a / b;
            if (slope > 1.0e-10) { // TODO: use global parameter for definition of "too small".
                xMin = std::ceil(std::max(-c / a, (double) 0));
                xMax = std::floor(std::min(h2 / slope - (c / a), (double) w2));
            }
            if (slope < -1.0e-10) { // TODO: use global parameter for definition of "too small".
                xMin = std::ceil(std::max(h2 / slope - (c / a), (double) 0));
                xMax = std::floor(std::min(-c / a, (double) w2));
            }
            // Drawing a match on the epipolar line.
            std::uniform_real_distribution<double> distributionUniformX2(xMin, xMax);
            x2 = distributionUniformX2(generator);
            y2 = (a * x2 + c) / (-b);

            // Drawing a direction to perturbate the perfect match.
            int direction = distributionUniform01(generator);
            double theta = std::atan(slope);
            if (direction == 0) {
                theta -= M_PI / (double) 2;
            } else {
                theta += M_PI / (double) 2;
            }
            while (theta < 0) {
                theta += 2 * M_PI;
            }
            while (theta > 2 * M_PI) {
                theta -= 2 * M_PI;
            }
            double minOffSet, maxOffSet;

            if (addOutlier == 2) {
                minOffSet = std::sqrt(outlierThresholdSq);
            } else {
                minOffSet = 0;
            }

            // To generate a second point in the image, the space is cut in 4 parts, to determine the available
            // distances between the match and the image border.
            double theta1, theta2, theta3, theta4;
            theta1 = std::atan((h2 - y2) / (double) (w2 - x2));
            theta2 = std::atan(x2 / (double) (h2 - y2)) + M_PI / (double) 2;
            theta3 = std::atan((y2) / (double) (x2)) + M_PI;
            theta4 = std::atan((w2 - x2) / (double) (y2)) + M_PI * 3 / (double) 2;

            assert(theta1 > 0 && theta1 < M_PI / (double) 2);
            assert(theta2 > M_PI / (double) 2 && theta2 < M_PI);
            assert(theta3 > M_PI && theta3 < M_PI * 3 / (double) 2);
            assert(theta4 > M_PI * 3 / (double) 2 && theta4 < 2 * M_PI);

            // A direction is sampled and the distance between the match and the image border is computed.
            if (theta >= theta4 || theta < theta1) {
                maxOffSet = (double) (w2 - x2) / std::cos(theta);
            }
            if (theta >= theta1 && theta < theta2) {
                maxOffSet = (double) (h2 - y2) / std::cos(theta - (M_PI / (double) 2));
            }
            if (theta >= theta2 && theta < theta3) {
                maxOffSet = (double) (x2) / std::cos(theta - M_PI);
            }
            if (theta >= theta3 && theta < theta4) {
                maxOffSet = (double) (y2) / std::cos(theta - (M_PI * 3 / (double) 2));
            }

            if (maxOffSet >= minOffSet) {

                // Drawing an error between the outlier and the epipolar line.
                std::uniform_real_distribution<double> distributionUniformOutlierError(minOffSet, maxOffSet);
                double outlierError = distributionUniformOutlierError(generator);
                x2 += outlierError * std::cos(theta);
                y2 += outlierError * std::sin(theta);

                assert(x2 >= 0 && x2 <= w2);
                assert(y2 >= 0 && y2 <= h2);

                Match outlier(static_cast<float>(x1),
                              static_cast<float>(y1),
                              static_cast<float>(x2),
                              static_cast<float>(y2));

                std::vector<Match> outlierVect{outlier};

                // The point is kept if addOutlier != 2 or if addOutlier = 2 and the point is not an inlier:
                assert(!((addOutlier == 2) &&
                         (model->Error(params, Match::toMat(outlierVect)) <= outlierThresholdSq)));

                numOutlierCreated += 1;

                artificialOutliers.push_back(outlier);

            }
        }
        // Special cases when the epipolar line is horizontal.
        if (std::abs(b) < 1.0e-10 && std::abs(a) > 1.0e-10) {
            double xTheorique = -c / a;
            std::uniform_real_distribution<double> distributionUniformX2(0, w2 - 2 * std::sqrt(outlierThresholdSq));
            x2 = distributionUniformX2(generator);
            if (x2 > xTheorique - std::sqrt(outlierThresholdSq)) {
                x2 += 2 * std::sqrt(outlierThresholdSq);
            }
            std::uniform_real_distribution<double> distributionUniformY2(0, h2);
            y2 = distributionUniformY2(generator);

            assert(x2 >= 0 && x2 <= w2);
            assert(y2 >= 0 && y2 <= h2);

            Match outlier(static_cast<float>(x1),
                          static_cast<float>(y1),
                          static_cast<float>(x2),
                          static_cast<float>(y2));

            std::vector<Match> outlierVect{outlier};

            // The point is kept if addOutlier != 2 or if addOutlier = 2 and the point is not an inlier:
            assert(!((addOutlier == 2) &&
                     (model->Error(params, Match::toMat(outlierVect)) <= outlierThresholdSq)));
            numOutlierCreated += 1;

            artificialOutliers.push_back(outlier);
        }
        // Special cases when the epipolar line is vertical.
        if (std::abs(a) < 1.0e-10 && std::abs(b) > 1.0e-10) {
            std::uniform_real_distribution<double> distributionUniformX2(0, h2);
            x2 = distributionUniformX2(generator);

            double yTheorique = -c / b;
            std::uniform_real_distribution<double> distributionUniformy2(0, h2 - 2 * std::sqrt(outlierThresholdSq));
            y2 = distributionUniformy2(generator);
            if (y2 > yTheorique - std::sqrt(outlierThresholdSq)) {
                y2 += 2 * std::sqrt(outlierThresholdSq);
            }

            assert(x2 >= 0 && x2 <= w2);
            assert(y2 >= 0 && y2 <= h2);

            Match outlier(static_cast<float>(x1),
                          static_cast<float>(y1),
                          static_cast<float>(x2),
                          static_cast<float>(y2));

            std::vector<Match> outlierVect{outlier};

            // The point is kept if addOutlier != 2 or if addOutlier = 2 and the point is not an inlier:
            assert(!((addOutlier == 2) &&
                     (model->Error(params, Match::toMat(outlierVect)) <= outlierThresholdSq)));
            numOutlierCreated += 1;

            artificialOutliers.push_back(outlier);
        }
    }
    return (numOutlierCreated == numOutlierWanted);
}

/// Generate random matches in two images. Those matches can be chosen to be outliers of a given model. Those matches are generated to ensure a more realistic distribution of outliers in the error space. Only works for line error space. TODO
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
bool generatePNPOutlierUniform(const double xmin, const double xmax,
                               const double ymin, const double ymax,
                               const double zmin, const double zmax,
                               const int w, const int h,
                               const double outlierRatio, const int inlierCount,
                               const int addOutlier,
                               const orsa::ModelEstimator *model, const libNumerics::matrix<double> &params,
                               const libNumerics::matrix<double> &intrinsics,
                               const double outlierThresholdSq,
                               std::vector<Match2D3D> &artificialOutliers,
                               std::default_random_engine &generator,
                               const int maxIterOutlier) {
    int iter = 0;
    int numOutlierCreated = 0;
    int numOutlierWanted;

    // Definition of the basic generators: for the first point in the left image and the orientation (above or under the
    // epipolar line).
    std::uniform_real_distribution<double> distributionUniformX(xmin, xmax);
    std::uniform_real_distribution<double> distributionUniformY(ymin, ymax);
    std::uniform_real_distribution<double> distributionUniformZ(zmin, zmax);
    std::uniform_real_distribution<long double> distributionUniformOrientation(0, 2 * M_PI);


    // Selection of the number of outlier to generate:
    if (outlierRatio < 1.0) {
        numOutlierWanted = static_cast<int>(std::floor(outlierRatio * inlierCount / (1 - outlierRatio)));
    } else {
        numOutlierWanted = static_cast<int>(outlierRatio);
    }

    while ((numOutlierCreated < numOutlierWanted) && (iter < numOutlierWanted + maxIterOutlier)) {

        // Generation of a point in the left image and its epipolar line.
        double x, y, z;
        x = distributionUniformX(generator);
        y = distributionUniformY(generator);
        z = distributionUniformZ(generator);
        Match2D3D tempMatch(0, 0, static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
        orsa::convert3Dto2D(params, intrinsics, tempMatch);
        double u = tempMatch.p2D(0), v = tempMatch.p2D(1);


        if (u >= 0 && u < w && v >= 0 && v < h) {

            double minOffSet, maxOffSet;

            // If addOutlier == 2 : we define a minimum distance between the perfect match and the outlier match.
            if (addOutlier == 2) {
                minOffSet = std::sqrt(outlierThresholdSq);
            } else {
                minOffSet = 0;
            }

            // To generate a second point in the image, the space is cut in 4 parts, to determine the available
            // distances between the match and the image border.
            double theta1, theta2, theta3, theta4;
            theta1 = std::atan((h - v) / (double) (w - u));
            theta2 = std::atan(u / (double) (h - v)) + M_PI / (double) 2;
            theta3 = std::atan((v) / (double) (u)) + M_PI;
            theta4 = std::atan((w - u) / (double) (v)) + M_PI * 3 / (double) 2;

            assert(theta1 > 0 && theta1 < M_PI / (double) 2);
            assert(theta2 > M_PI / (double) 2 && theta2 < M_PI);
            assert(theta3 > M_PI && theta3 < M_PI * 3 / (double) 2);
            assert(theta4 > M_PI * 3 / (double) 2 && theta4 < 2 * M_PI);

            // A direction is sampled and the distance between the match and the image border is computed.
            double theta = distributionUniformOrientation(generator);

            if (theta >= theta4 || theta < theta1) {
                maxOffSet = (double) (w - u) / std::cos(theta);
            }
            if (theta >= theta1 && theta < theta2) {
                maxOffSet = (double) (h - v) / std::cos(theta - (M_PI / (double) 2));
            }
            if (theta >= theta2 && theta < theta3) {
                maxOffSet = (double) (u) / std::cos(theta - M_PI);
            }
            if (theta >= theta3 && theta < theta4) {
                maxOffSet = (double) (v) / std::cos(theta - (M_PI * 3 / (double) 2));
            }

            if (maxOffSet >= minOffSet) {
                // A point is drawn along the drawn direction.
                std::uniform_real_distribution<double> distributionUniformOutlierError(minOffSet, maxOffSet);
                double outlierError = distributionUniformOutlierError(generator);

                u += std::cos(theta) * outlierError;
                v += std::sin(theta) * outlierError;

                assert(u >= 0 && u <= w);
                assert(v >= 0 && v <= h);

                Match2D3D outlier(static_cast<float>(u), static_cast<float>(v),
                                  tempMatch.p3D(0), tempMatch.p3D(1), tempMatch.p3D(2));

                std::vector<Match2D3D> outlierVect{outlier};

                // The point is kept if addOutlier != 2 or if addOutlier = 2 and the point is not an inlier:
                assert(!((addOutlier == 2) &&
                         (model->Error(params, Match2D3D::toMat(outlierVect)) <= outlierThresholdSq)));
                numOutlierCreated += 1;

                artificialOutliers.push_back(outlier);
            }
        }
        iter += 1;
    }
    return (numOutlierCreated == numOutlierWanted);
}

/// Add noise to matches of a vector while staying in a specified range.
/// \param[in] vectToModify: Vector on which to apply the noising.
/// \param[in] width: Maximum value of x.
/// \param[in] height: Maximum value of y.
/// \param[out] vectModified: Vector of modified matches.
/// \param[in] std: Standard deviation of the gaussian noise or range of the uniform noise. (Default: 0)
/// \param[in] gaussian: If true gaussian noise is used, uniform noise otherwise. (Default: false)
void safeAddNoise(const std::vector<Match> &vectToModify,
                  const int width, const int height,
                  std::vector<Match> &vectModified,
                  std::default_random_engine &generator,
                  const double std, const bool gaussian,
                  const int numInliers) {

    std::vector<int> inlierSample;
    if (numInliers > 0) {
        orsa::UniformSample(numInliers, vectToModify.size(), &inlierSample);
    } else {
        for (int i = 0; i < (int) vectToModify.size(); i++) {
            inlierSample.push_back(i);
        }
    }

    if (std > 0) {
        std::vector<int>::iterator itInlierSample = inlierSample.begin();
        for (; itInlierSample != inlierSample.end(); ++itInlierSample) {
            Match matchToModify = vectToModify[*itInlierSample];
            double x2 = static_cast<double>(matchToModify.x2);
            double y2 = static_cast<double>(matchToModify.y2);
            addNoiseBounded(std, x2, y2, generator, width, height, gaussian);

            Match newMatch = matchToModify;
            newMatch.x2 = static_cast<float>(x2);
            newMatch.y2 = static_cast<float>(y2);
            vectModified.push_back(newMatch);
        }
    } else {
        std::vector<int>::iterator itInlierSample = inlierSample.begin();
        for (; itInlierSample != inlierSample.end(); ++itInlierSample) {
            vectModified.push_back(vectToModify[*itInlierSample]);
        }
    }
}

/// Add noise to matches of a vector while staying in a specified range.
/// \param[in] vectToModify: Vector on which to apply the noising.
/// \param[in] width: Maximum value of x.
/// \param[in] height: Maximum value of y.
/// \param[out] vectModified: Vector of modified matches.
/// \param[in] std: Standard deviation of the gaussian noise or range of the uniform noise. (Default: 0)
/// \param[in] gaussian: If true gaussian noise is used, uniform noise otherwise. (Default: false)
void safeAddNoise(const std::vector<Match2D3D> &vectToModify,
                  const int width, const int height,
                  std::vector<Match2D3D> &vectModified,
                  std::default_random_engine &generator,
                  const double std, const bool gaussian,
                  const int numInliers) {

    std::vector<int> inlierSample;
    if (numInliers > 0) {
        orsa::UniformSample(numInliers, vectToModify.size(), &inlierSample);
    } else {
        for (int i = 0; i < (int) vectToModify.size(); i++) {
            inlierSample.push_back(i);
        }
    }

    if (std > 0) {
        std::vector<int>::iterator itInlierSample = inlierSample.begin();
        for (; itInlierSample != inlierSample.end(); ++itInlierSample) {
            Match2D3D matchToModify = vectToModify[*itInlierSample];
            double u = static_cast<double>(matchToModify.p2D(0));
            double v = static_cast<double>(matchToModify.p2D(1));
            addNoiseBounded(std, u, v, generator, width, height, gaussian);

            Match2D3D newMatch = matchToModify;
            newMatch.p2D(0) = static_cast<float>(u);
            newMatch.p2D(1) = static_cast<float>(v);
            vectModified.push_back(newMatch);
        }
    } else {
        std::vector<int>::iterator itInlierSample = inlierSample.begin();
        for (; itInlierSample != inlierSample.end(); ++itInlierSample) {
            vectModified.push_back(vectToModify[*itInlierSample]);
        }
    }
}

/// Find the maximum error of a model.
/// \param[in] vectIdxToEvaluate: The vector of index of points to evaluate in model.
/// \param[in] model: The model the data fits.
/// \param[in] modelParams: The model parameters to use for computation.
/// \return The squared maximum error found.
double findMaxErrorSq(const std::vector<int> &vectIdxToEvaluate,
                      const orsa::ModelEstimator *model, const libNumerics::matrix<double> &modelParams) {
    double errorMax = 0;
    std::vector<int>::const_iterator itIdxVectToEvaluate = vectIdxToEvaluate.begin();
    for (; itIdxVectToEvaluate != vectIdxToEvaluate.end(); ++itIdxVectToEvaluate) {
        double error = model->Error(modelParams, *itIdxVectToEvaluate);
        if (error > errorMax) {
            errorMax = error;
        }
    }
    return errorMax;
}

void findOutlier2D3DRange(const std::vector<Match2D3D> &matches,
                          double &xMin, double &xMax,
                          double &yMin, double &yMax,
                          double &zMin, double &zMax,
                          const double increaseRatio) {
    xMin = std::numeric_limits<double>::infinity();
    xMax = -1 * std::numeric_limits<double>::infinity();
    yMin = std::numeric_limits<double>::infinity();
    yMax = -1 * std::numeric_limits<double>::infinity();
    zMin = std::numeric_limits<double>::infinity();
    zMax = -1 * std::numeric_limits<double>::infinity();
    std::vector<Match2D3D>::const_iterator itMatches = matches.begin();
    for (; itMatches != matches.end(); ++itMatches){
        double x, y, z;
        x = itMatches->p3D(0);
        y = itMatches->p3D(1);
        z = itMatches->p3D(2);
        if (x < xMin) {
            xMin = x;
        }
        if (x > xMax) {
            xMax = x;
        }
        if (y < yMin) {
            yMin = y;
        }
        if (y > yMax) {
            yMax = y;
        }
        if (z < zMin) {
            zMin = z;
        }
        if (z > zMax) {
            zMax = z;
        }
    }
    xMin *= increaseRatio;
    xMax *= increaseRatio;
    yMin *= increaseRatio;
    yMax *= increaseRatio;
    zMin *= increaseRatio;
    zMax *= increaseRatio;
}