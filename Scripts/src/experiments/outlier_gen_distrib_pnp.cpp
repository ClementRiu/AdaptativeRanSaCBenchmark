//
// Created by riuclement on 9/9/21.
//

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>

#include "libOrsa/pnp_model.hpp"

#include "libOrsa/ransac.hpp"
#include "libOrsa/orsa.hpp"
#include "libOrsa/lrtsac.hpp"
#include "libOrsa/muse.hpp"
#include "libOrsa/orsa_fast.hpp"
#include "libOrsa/starsac.hpp"
#include "libOrsa/magsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/data_handler.hpp"
#include "utilities/metrics.hpp"
#include "utilities/usac_file_handler.hpp"

#include "apply_transform.hpp"
#include "read_write_experiment.hpp"


template<typename T>
bool saveVectOfVect(const char *nameFile, const std::vector<std::vector<T>> &vectOfVect) {
    std::ofstream f(nameFile);
    if (f.is_open()) {
        typename std::vector<std::vector<T>>::const_iterator itVectOfVect = vectOfVect.begin();
        for (; itVectOfVect != vectOfVect.end(); ++itVectOfVect) {
            typename std::vector<T>::const_iterator itVect = (*itVectOfVect).begin();
            for (; itVect != (*itVectOfVect).end(); ++itVect) {
                f << *itVect << ", ";
            }
            f << "\n";
        }
    }
    return f.is_open();
}

const char *concatenateTextValue(char *text, const double value) {
    std::string valueStr = std::to_string(value);
    return strcat(text, valueStr.c_str());
}

int main(int argc, char **argv) {
    int nGen = 1; // Number of different datasets to generate.
    int nRun = 1; // Number of run for each different datasets.

    int maxMatchNumber = 0; // Maximum number of matches in the dataset. If 0, no limit is set.
    // Inlier noise parameters:
    double stdNoise = 1.0;
    double increaseRatio = 1.1; // The increase of the 3D region compared to the inlier bounding box.
    // Random outlier parameters:
    int outlierType = 2; // If 0, no outliers,
    // if 1 add outliers that can be in the inlier region,
    // if 2 add outliers that are true outliers.
    double outlierRatio = 1000; // Outlier/Inlier ratio.
    // If integer, the number of outlier to generate,
    // if decimal between 0 and 1 (excluded) ratio.
    int maxIterOutlier = 100000; // Maximum number of iterations to generate the outliers.
    // Must be set as outliers can be refused when in inlier region.

    int modelToAnalyse = 0; // 0 for PnP.

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    utility::CmdLine cmd;
    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for PnP."));

    cmd.add(utility::make_option('g', nGen, "num-gen-exp")
                    .doc("Number of noisy datasets generated for the experiment."));
    cmd.add(utility::make_option('n', nRun, "num-run")
                    .doc("Number of run of the algorithm."));
    cmd.add(utility::make_option(0, maxMatchNumber, "max-match")
                    .doc("Maximum number of matches in the dataset."));
    cmd.add(utility::make_option('s', stdNoise, "noise-level")
                    .doc("Value of the noise std. If > 0 noise will be added."));
    cmd.add(utility::make_option(0, increaseRatio, "increase-ratio")
                    .doc("The increase of the 3D region compared to the inlier bounding box."));
    cmd.add(utility::make_switch('u', "uniform-noise")
                    .doc("Use uniform noise or gaussian noise."));
    cmd.add(utility::make_option('o', outlierType, "outlier-type")
                    .doc("Add outliers: 0 for no, 1 for uniform outliers, 2 for uniform outliers that can't be inliers."));
    cmd.add(utility::make_option('r', outlierRatio, "outlier-ratio")
                    .doc("Ratio of outlier or number: if in [0, 1[ then it is the ratio, if integer equal or greater than 1 it is the number of outlier."));

    cmd.add(utility::make_switch('f', "force-compute")
                    .doc("Run the computation even if a file already exists."));
    cmd.add(utility::make_switch('v', "verbose")
                    .doc("Print info during the run."));
    cmd.add(utility::make_option('t', seed, "time-seed")
                    .doc("Use value instead of time for random seed (for debug)."));

    try {
        cmd.process(argc, argv);
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return 1;
    }

    int minArgc = 7;

    if (argc != minArgc) {
        std::cerr << "Usage: " << argv[0]
                  << " [options] good_matches.txt info.txt calib.txt outputInfo.txt [noisyInliers.txt [artificialOutliers.txt]]\n"
                  << "- good_matches.txt: path to the good matches to read.\n"
                  << "- info.txt: path to the run info (like dimensions of images).\n"
                  << "- calib.txt: path to the calib matrix of both images.\n"
                  << "- outputInfo.txt: path to the output file with all the informations on the experiment.\n"
                  << "- noisyInliers.txt: path to the output file with the inliers after noise addition.\n"
                  << "- artificialOutliers.txt: path to the output file with the artificial outliers.\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }

    bool forceCompute = cmd.used('f');
    bool verbose = cmd.used('v');
    bool gaussian = !cmd.used('u');
    // Init random seed
    std::srand(seed);

    const char *pathInGoodMatches = argv[1]; // Output match file of a generate_artificial_dataset run.
    const char *pathInInfo = argv[2]; // Output info file of a generate_artificial_dataset run.
    const char *pathInCalib = argv[3]; // Essential only: calibration matrix of both images
    const char *pathOutput = argv[4]; // Output file with all metrics and info.
    const char *pathOutOutlierSmartError = argv[5];
    const char *pathOutOutlierUniError = argv[6];

    std::ifstream test_f(pathOutput);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    std::cout << "Random seed: " << seed << "\n" << std::endl;

    std::cout << "\nReading " << pathInGoodMatches << std::endl;

    std::vector<Match2D3D> inlierMatches;
    if (Match2D3D::loadMatch2D3D(pathInGoodMatches, inlierMatches)) {
        std::cout << "Read " << inlierMatches.size() << " matches" << std::endl;
    } else {
        std::cerr << "Failed reading matches from " << pathInGoodMatches << std::endl;
        return 1;
    }

    std::cout << "\nReading " << pathInInfo << std::endl;
    int w, h, zero1, zero2;
    libNumerics::matrix<double> trueModelParams(3, 4);
    int inlierCount;
    if (readArtificialGenerationParams(pathInInfo, w, h, zero1, zero2, trueModelParams, inlierCount)) {
        std::cout << "Read experiment information.\n" << std::endl;
    } else {
        std::cerr << "Failed reading experiment information from " << pathInInfo << std::endl;
        return 1;
    }

    libNumerics::matrix<double> intrinsics(3, 3); // The intrinsic parameters of the source camera
    if (modelToAnalyse == 0) {
        if (utility::loadCalibration(pathInCalib, intrinsics))
            std::cout << "Read " << pathInCalib << " calibration file." << std::endl;
        else {
            std::cerr << "Failed reading calibration matrix from " << pathInCalib << std::endl;
            return 1;
        }
    }
    // Generic algorithm pointer that will be used to call each algorithm.
    orsa::RansacAlgorithm *algorithm;

    // Preparation of the random generator:
    std::default_random_engine generator;
    generator.seed(seed);
    // Generate multiple semi-artificial dataset.
    for (int gen = 0; gen < nGen; gen++) {
        std::vector<Match2D3D> noisyInlierMatches;
        int numInliers = inlierCount;
        if (outlierRatio < 1) {
            numInliers = std::min(inlierCount, (int) ((double) maxMatchNumber * (1 - outlierRatio)));
        }

        // Adding noise to the perfect inlier matches and removing matches if necessary.
        safeAddNoise(inlierMatches, w, h, noisyInlierMatches, generator, stdNoise, gaussian, numInliers);

        if (verbose) {
            std::cout << "Number of noisy inliers: " << noisyInlierMatches.size() << std::endl;
        }

        std::vector<Match2D3D> normalisedNoisyInlierMatches;
        orsa::ModelEstimator *tempModel;

        // Temporary model to create the artificial outliers.
        switch (modelToAnalyse) {
            case 0:
                orsa::normalisePoints(noisyInlierMatches, intrinsics, normalisedNoisyInlierMatches);

                tempModel = new orsa::PnPModel(noisyInlierMatches, normalisedNoisyInlierMatches, intrinsics, w, h);
                break;
            default:
                std::cerr << "Model number not supported. Choose in {0}" << std::endl;
                return 1;
                break;
        }

        double xMin, xMax, yMin, yMax, zMin, zMax;
        findOutlier2D3DRange(noisyInlierMatches, xMin, xMax, yMin, yMax, zMin, zMax, increaseRatio);
        // Find the true inlier/outlier threshold.
        double outlierThresholdSq = findMaxErrorSq(noisyInlierMatches, tempModel, trueModelParams);

        std::vector<Match2D3D> matchesToEvaluate;
        std::vector<int> labels;

        // Create the artificial outliers, depending on the model and merge inliers and outliers
        // in a matchesToEvaluate vector with labels.
        std::vector<Match2D3D> artificialOutlierSmarts;
        std::vector<Match2D3D> artificialOutlierUnis;
        if (outlierType > 0 && outlierRatio > 0.0) {
            switch (modelToAnalyse) {
                case 0:
                    if (!generatePNPOutlierUniform(xMin, xMax, yMin, yMax, zMin, zMax, w, h,
                                                   outlierRatio, noisyInlierMatches.size(), outlierType,
                                                   tempModel, trueModelParams, intrinsics,
                                                   outlierThresholdSq, artificialOutlierSmarts, generator, maxIterOutlier)) {
                        std::cerr << "Not enough outliers created: " << artificialOutlierSmarts.size() << std::endl;
                    } else {
                        if (verbose) {
                            std::cout << "Number of outliers: " << artificialOutlierSmarts.size() << std::endl;
                        }

                        utility::randomDataFusionWithMemory(noisyInlierMatches, artificialOutlierSmarts, matchesToEvaluate,
                                                            1, 0,
                                                            labels);
                    }
                    break;
                default:
                    std::cerr << "Model number not supported. Choose in {0}" << std::endl;
                    return 1;
                    break;

            }
        } else {
            std::vector<Match2D3D>::const_iterator itNoisyInlierMatches = noisyInlierMatches.begin();
            for (; itNoisyInlierMatches != noisyInlierMatches.end(); itNoisyInlierMatches++) {
                matchesToEvaluate.push_back(*itNoisyInlierMatches);
                labels.push_back(1);
            }
        }


        std::vector<double> errorsOutlierUnis;

        std::uniform_real_distribution<double> distributionUniformWidth1(xMin, xMax);
        std::uniform_real_distribution<double> distributionUniformHeight1(yMin, yMax);
        std::uniform_real_distribution<double> distributionUniformDepth1(zMin, zMax);
        std::uniform_real_distribution<double> distributionUniformWidth2(0, w);
        std::uniform_real_distribution<double> distributionUniformHeight2(0, h);
        int iter = 0;
        int numOutlierCreated = 0;
        int numOutlierWanted;
        if (outlierRatio < 1.0) {
            numOutlierWanted = static_cast<int>(std::floor(outlierRatio * inlierCount / (1 - outlierRatio)));
        } else {
            numOutlierWanted = static_cast<int>(outlierRatio);
        }

        while ((numOutlierCreated < numOutlierWanted) && (iter < numOutlierWanted + maxIterOutlier)) {
            double x1 = distributionUniformWidth1(generator);
            double y1 = distributionUniformHeight1(generator);
            double z1 = distributionUniformDepth1(generator);
            double x2 = distributionUniformWidth2(generator);
            double y2 = distributionUniformHeight2(generator);
            Match2D3D testOutlier(static_cast<float>(x1),
                                  static_cast<float>(y1),
                                  static_cast<float>(z1),
                              static_cast<float>(x2),
                              static_cast<float>(y2));
            std::vector<Match2D3D> outlierVect{testOutlier};
            double testError = tempModel->Error(trueModelParams, Match2D3D::toMat(outlierVect));
            if (testError > outlierThresholdSq){
                errorsOutlierUnis.push_back(std::sqrt(testError));
                artificialOutlierUnis.push_back(testOutlier);
                numOutlierCreated ++;
            }
            iter ++;
        }


        std::vector<double> errorsOutlierSmarts;
        std::vector<Match2D3D>::const_iterator itOutlierSmart = artificialOutlierSmarts.begin();
        for (; itOutlierSmart != artificialOutlierSmarts.end(); itOutlierSmart++) {
            typename std::vector<Match2D3D> outlierVect{*itOutlierSmart};
            errorsOutlierSmarts.push_back(std::sqrt(tempModel->Error(trueModelParams, Match2D3D::toMat(outlierVect))));
        }
        std::ofstream foutSmart(pathOutOutlierSmartError);
        if (foutSmart.is_open()) {
            std::vector<double>::const_iterator VectIt = errorsOutlierSmarts.begin();
            for (; VectIt != errorsOutlierSmarts.end(); ++VectIt) {
                foutSmart << *VectIt << " ";
            }
        }
        foutSmart.close();
        std::ofstream foutUni(pathOutOutlierUniError);
        if (foutUni.is_open()) {
            std::vector<double>::const_iterator VectIt = errorsOutlierUnis.begin();
            for (; VectIt != errorsOutlierUnis.end(); ++VectIt) {
                foutUni << *VectIt << " ";
            }
        }
        foutUni.close();

        delete tempModel;

    } // End of the 'gen' loop.

    return 0;
}
