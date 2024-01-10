//
// Created by clementriu on 8/31/20.
//

#ifndef MMM_ORSA_READ_WRITE_EXPERIMENT_HPP
#define MMM_ORSA_READ_WRITE_EXPERIMENT_HPP

#include "libOrsa/libNumerics/matrix.h"
#include "libOrsa/match.hpp"

//// Save all parameters of an artificial match generation to a file.
//// \param[in] fileName: Emplacement of the file.
//// \param[in] seed: Value of the random seed.
//// \param[in] iterMax: Maximum number of iteration allowed during computation.
//// \param[in] noiseType: When inlier noise is applied : 0 for no noise, 1 after model, -1 before model.
//// \param[in] stdNoise: Std of the inlier noise.
//// \param[in] outlierType: Outlier type : 0 for no outliers, 1 for uniform outliers, 2 for uniform outliers that can't be inliers.
//// \param[in] outlierRatio: Ratio of outliers.
//// \param[in] outlierThreshold: Inlier/Outlier threshold.
//// \param[in] maxOutlier: Maximum number of iterations when adding outliers.
//// \param[in] w1: h1: w2: h2: Dimensions of images 1 and 2.
//// \param[in] modelParams: Used model matrix.
//// \param[in] modelError: Error of the used model.
//// \param[in] inlierCount: Number of inliers of the model.
//// \param[in] computedSigma: Estimated inlier/outlier threshold.
//// \return True if file was successfully opened, False, otherwise.
bool saveArtificialGenerationParams(const char *fileName, unsigned int seed,
                                    double precision, int iterMax,
                                    int noiseType, double stdNoise,
                                    int outlierType, double outlierRatio, double outlierThreshold,
                                    int maxOutlier,
                                    int w1, int h1, int w2, int h2,
                                    const libNumerics::matrix<double> &modelParams, double modelError,
                                    int inlierCount, double computedSigma);

//// Save all parameters of an artificial match generation to a file.
//// \param[in] fileName: Emplacement of the file.
//// \param[out] w1: h1: w2: h2: Dimensions of images 1 and 2.
//// \param[out] modelParams: Used model matrix.
//// \param[out] inlierCount: Number of inliers of the model.
//// \return True if file was successfully opened, False, otherwise.
bool readArtificialGenerationParams(const char *fileName,
                                    int &w1, int &h1, int &w2, int &h2,
                                    libNumerics::matrix<double> &modelParams, int &inlierCount);

//// Read a set of points and mixes them.
//// \param[in] fileInliers: Path to the inlier matches.
//// \param[in] fileOutliers: Path to the oulier matches.
//// \param[in] nGen: Number of dataset to read.
//// \param[out] pointsAll: All points read.
//// \param[out] groundTruthLabelsAll: All labels: 0 if from outlier, 1 if from inlier.
//// \param[out] inliersAll: Only inliers.
//// \param[out] outliersAll: Only outliers.
//// \param[in] readOutliers: If False, no outliers will be read.
bool ReadPoints(const char *fileInliers, const char *fileOutliers, int nGen,
                std::vector<std::vector<Match>> &pointsAll, std::vector<std::vector<int>> &groundTruthLabelsAll,
                std::vector<std::vector<Match>> &inliersAll, std::vector<std::vector<Match>> &outliersAll,
                bool readOutliers);

#endif //MMM_ORSA_READ_WRITE_EXPERIMENT_HPP
