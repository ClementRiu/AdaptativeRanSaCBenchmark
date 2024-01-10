//
// Created by riuclement on 3/22/21.
//

#include <cstdlib>
#include <ctime>
#include <iostream>

#include "libOrsa/homography_model.hpp"
#include "libOrsa/fundamental_model.hpp"
#include "libOrsa/essential_model.hpp"

#include "libOrsa/lrtsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/data_handler.hpp"
#include "utilities/metrics.hpp"
#include "utilities/usac_file_handler.hpp"

#include "apply_transform.hpp"
#include "read_write_experiment.hpp"


/// File to run an experiment on an semi-artificial dataset:
/// - Add noise to the inliers.
/// - Creates outliers.
/// - Run Ransac with two different thresholds, AC-RANSAC and LRTSac.
/// - Generates different datasets and run multiple time each algorithms on each dataset.

//// Function to estimate and evaluate a model parameters with a given algorithm.
//// \param model: A ModelEstimator pointer that contains the matches.
//// \param algorithm: An RansacAlgorithm pointer that is linked to the model.
//// \param labels: The vector of match labels: 1 for inliers, 0 for outliers.
//// \param noisyInlierMatches: The vector of semi-artificial inlier matches.
//// \param algoNumber: The index of the algorithm: 0 and 1 should be RANSAC, 1 AC-RANSAC, 2 LRTSac.
//// \param sigmaEstimatedVectVect: For each algorithm the estimated sigmas are saved for each run.
//// First dimension is the algorithm, second dimension is the run.
//// \param runtimeVectVect: For each algorithm the runtime are saved for each run.
////// First dimension is the algorithm, second dimension is the run.
//// \param precisionVectVect: For each algorithm the precision are saved for each run.
////// First dimension is the algorithm, second dimension is the run.
//// \param recallVectVect: For each algorithm the recall are saved for each run.
////// First dimension is the algorithm, second dimension is the run.
//// \param algorithmNames: Name of the different algorithms for printing.
//// \param iterMax: Maximum number of iterations for the run.
//// \param verbose: Verbose setting.
bool runExperiment(const orsa::ModelEstimator *model, const orsa::RansacAlgorithm *algorithm,
                   const std::vector<int> &labels, const std::vector<Match> &noisyInlierMatches,
                   const int algoNumber,
                   std::vector<std::vector<double>> &sigmaEstimatedVectVect,
                   std::vector<std::vector<double>> &runtimeVectVect,
                   std::vector<std::vector<double>> &precisionVectVect,
                   std::vector<std::vector<double>> &recallVectVect,
                   const std::vector<const char *> &algorithmNames,
                   const int iterMax, const bool verbose) {
    // Output values:
    bool ok = false;
    double runtime;
    orsa::RansacAlgorithm::RunResult res;

    // Run the given model with given algorithm:
    ok = algorithm->evalModel(res, runtime, iterMax, verbose);

    double error=0;
    if (ok)
        error = model->ErrorStats(res.vInliers, res.model).first;

    int numTruePositives = 0;
    double precision = 0;
    double recall = 0;
    int inlierCount = res.vInliers.size();

    // Evaluate the estimated model:
    if (ok) {
        sigmaEstimatedVectVect[algoNumber].push_back(res.sigma);
        runtimeVectVect[algoNumber].push_back(runtime);

        numTruePositives = utility::computeTruePositive(res.vInliers, labels);
        precision = utility::computePrecision(numTruePositives, res.vInliers.size());
        recall = utility::computeRecall(numTruePositives, noisyInlierMatches.size());
    }

    precisionVectVect[algoNumber].push_back(precision);
    recallVectVect[algoNumber].push_back(recall);

    if (ok && verbose) {
        std::cout << algorithmNames[algoNumber] << std::endl;
        std::cout << "\tModelParam=" << res.model << std::endl;
        std::cout << "\tInlier error: " << error << "\n\tComputed sigma: " << res.sigma
                  << "\n\tNumber of inliers: " << inlierCount
                  << "\n\tNumber of iterations: " << res.T << "\n\tVPM: " << res.vpm
                  << "\n\tRuntime: " << runtime
                  << std::endl;
        std::cout << "\tNumber true positives: " << numTruePositives
                  << "\n\tPrecision: " << precision
                  << "\n\tRecall: " << recall
                  << std::endl;
    }

    if (!ok && verbose) {
        std::cerr << "Failed to estimate a model with " << algorithmNames[algoNumber] << std::endl;
    }

    return ok;
}

const char *concatenateTextValue(char *text, const double value) {
    std::string valueStr = std::to_string(value);
    return strcat(text, valueStr.c_str());
}

int main(int argc, char **argv) {
    const int nAlgo = 2; // Number of algorithms in this experiment.

    int iterMax = 10000; // Maximum number of iterations of the algorithm.
    double maxPrecision = 16; // Default pixel precision of the algorithm.
    double cpIIT = 0.99; // Default confidence of the algorithm.

    double maxSigmaLRT = maxPrecision; // LRT specific precision /// Redundant initialisation for clarity.
    double cpI = 0.99;
    double cpIIB = 0.95;

    int nGen = 1; // Number of different datasets to generate.
    int nRun = 1; // Number of run for each different datasets.

    int maxMatchNumber = 0; // Maximum number of matches in the dataset. If 0, no limit is set.
    // Inlier noise parameters:
    double stdNoise = 1.0;
    // Random outlier parameters:
    int outlierType = 0; // If 0, no outliers,
    // if 1 add outliers that can be in the inlier region,
    // if 2 add outliers that are true outliers.
    double outlierRatio = 0.0; // Outlier/Inlier ratio.
    // If integer, the number of outlier to generate,
    // if decimal between 0 and 1 (excluded) ratio.
    int maxIterOutlier = 100000; // Maximum number of iterations to generate the outliers.
    // Must be set as outliers can be refused when in inlier region.

    int modelToAnalyse = 0; // 0 for Homography, 1 for Fundamental, 2 for Essential.

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    utility::CmdLine cmd;
    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of algorithms."));

    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential."));

    cmd.add(utility::make_option('p', maxPrecision, "maximum-precision")
                    .doc("Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of AC-RANSAC and LRT."));
    cmd.add(utility::make_option(0, cpIIT, "cpIIT")
                    .doc("Value of confidence (Ransac/LRT)"));

    cmd.add(utility::make_option(0, cpI, "cpI")
                    .doc("Confidence proba wrt type I error (LRT)"));
    cmd.add(utility::make_option(0, cpIIB, "cpIIB")
                    .doc("Confidence proba wrt to bailout (LRT)"));

    cmd.add(utility::make_option('g', nGen, "num-gen-exp")
                    .doc("Number of noisy datasets generated for the ipol experience."));
    cmd.add(utility::make_option('n', nRun, "num-run")
                    .doc("Number of run of the algorithm."));
    cmd.add(utility::make_option(0, maxMatchNumber, "max-match")
                    .doc("Maximum number of matches in the dataset."));
    cmd.add(utility::make_option('s', stdNoise, "noise-level")
                    .doc("Value of the noise std. If > 0 noise will be added."));
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

    int minArgcHF = 6;
    int minArgcE = 7;

    if ((modelToAnalyse < 2 && argc != minArgcHF) || (modelToAnalyse == 2 && argc != minArgcE)) {
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
    const char *pathInCalib; // Essential only: calibration matrix of both images
    const char *pathOutput; // Output file with all metrics and info.
    const char *pathNoisyInliers; // Output file with semi-artificial inlier matches.
    const char *pathArtificialOutliers; // Output file with artificial outlier matches.

    if (modelToAnalyse < 2) {
        pathOutput = argv[3];
        pathNoisyInliers = argv[4];
        pathArtificialOutliers = argv[5];
    }
    if (modelToAnalyse == 2) {
        pathInCalib = argv[3];
        pathOutput = argv[4];
        pathNoisyInliers = argv[5];
        pathArtificialOutliers = argv[6];
    }

    std::ifstream test_f(pathOutput);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    maxSigmaLRT = maxPrecision; // LRT specific precision

    char LRTSacBailoutText[100] = "LRTSac with bailout and threshold ";
    char LRTSacText[100] = "LRTSac without bailout and threshold ";
    std::vector<const char *> algorithmNames = {
            concatenateTextValue(LRTSacBailoutText, maxSigmaLRT),
            concatenateTextValue(LRTSacText, maxSigmaLRT),
    };

    std::cout << "\nReading " << pathInGoodMatches << std::endl;

    std::vector<Match> inlierMatches;
    if (Match::loadMatch(pathInGoodMatches, inlierMatches)) {
        std::cout << "Read " << inlierMatches.size() << " matches" << std::endl;
    } else {
        std::cerr << "Failed reading matches from " << pathInGoodMatches << std::endl;
        return 1;
    }

    std::cout << "\nReading " << pathInInfo << std::endl;
    int w1, h1, w2, h2;
    libNumerics::matrix<double> trueModelParams(3, 3);
    int inlierCount;
    if (readArtificialGenerationParams(pathInInfo, w1, h1, w2, h2, trueModelParams, inlierCount)) {
        std::cout << "Read experiment information.\n" << std::endl;
    } else {
        std::cerr << "Failed reading experiment information from " << pathInInfo << std::endl;
        return 1;
    }

    libNumerics::matrix<double> intrinsics_source(3, 3), // The intrinsic parameters of the source camera
    intrinsics_destination(3, 3); // The intrinsic parameters of the destination camera
    if (modelToAnalyse == 2) {
        if (utility::loadCalibration(pathInCalib, intrinsics_source, intrinsics_destination))
            std::cout << "Read " << pathInCalib << " calibration file." << std::endl;
        else {
            std::cerr << "Failed reading matches from " << pathInCalib << std::endl;
            return 1;
        }
    }
    // Generic algorithm pointer that will be used to call each algorithm.
    orsa::RansacAlgorithm *algorithm;

    // All metrics are stored across each algorithm across each run.
    std::vector<std::vector<double>> precisionVectVect(nAlgo);
    std::vector<std::vector<double>> recallVectVect(nAlgo);
    std::vector<std::vector<double>> runtimeVectVect(nAlgo);
    std::vector<std::vector<double>> sigmaEstimatedVectVect(nAlgo);

    std::vector<std::vector<Match>> noisyInliersVect;
    std::vector<std::vector<Match>> artificialOutliersVect;

    // Preparation of the random generator:
    std::default_random_engine generator;
    generator.seed(seed);
    // Generate multiple semi-artificial dataset.
    for (int gen = 0; gen < nGen; gen++) {
        std::vector<Match> noisyInlierMatches;
        int numInliers = inlierCount;
        if (outlierRatio < 1) {
            numInliers = std::min(inlierCount, (int) ((double) maxMatchNumber * (1 - outlierRatio)));
        }

        // Adding noise to the perfect inlier matches and removing matches if necessary.
        safeAddNoise(inlierMatches, w2, h2, noisyInlierMatches, generator, stdNoise, gaussian, numInliers);

        if (verbose) {
            std::cout << "Number of noisy inliers: " << noisyInlierMatches.size() << std::endl;
        }

        std::vector<Match> normalisedNoisyInlierMatches;
        orsa::ModelEstimator *tempModel;

        // Temporary model to create the artificial outliers.
        switch (modelToAnalyse) {
            case 0:
                std::cout << "Running a Homography estimation..." << std::endl;
                tempModel = new orsa::HomographyModel(noisyInlierMatches, w1, h1, w2, h2, true);
                break;
            case 1:
                std::cout << "Running a Fundamental matrix estimation..." << std::endl;
                tempModel = new orsa::FundamentalModel(noisyInlierMatches, w1, h1, w2, h2, true);
                break;
            case 2:

                orsa::normalisePoints(noisyInlierMatches, intrinsics_source, intrinsics_destination,
                                      normalisedNoisyInlierMatches);

                tempModel = new orsa::EssentialModel(noisyInlierMatches, normalisedNoisyInlierMatches, w1, h1, w2, h2,
                                                     intrinsics_source,
                                                     intrinsics_destination, true);
                break;
            default:
                std::cerr << "Model number not supported. Choose in {0, 1, 2}" << std::endl;
                return 1;
                break;
        }

        // Find the true inlier/outlier threshold.
        double outlierThresholdSq = findMaxErrorSq(noisyInlierMatches, tempModel, trueModelParams);

        std::vector<Match> matchesToEvaluate;
        std::vector<int> labels;

        // Create the artificial outliers, depending on the model and merge inliers and outliers
        // in a matchesToEvaluate vector with labels.
        std::vector<Match> artificialOutliers;
        if (outlierType > 0 && outlierRatio > 0.0) {
            switch (modelToAnalyse) {
                case 0:
                    if (!generateHomOutliersUniformError(w1, h1, w2, h2,
                                                         outlierRatio, noisyInlierMatches.size(), outlierType,
                                                         tempModel, trueModelParams,
                                                         outlierThresholdSq, artificialOutliers, generator,
                                                         maxIterOutlier)) {
                        std::cerr << "Not enough outliers created: " << artificialOutliers.size() << std::endl;
                    } else {
                        if (verbose) {
                            std::cout << "Number of outliers: " << artificialOutliers.size() << std::endl;
                        }

                        utility::randomDataFusionWithMemory(noisyInlierMatches, artificialOutliers, matchesToEvaluate,
                                                            1, 0,
                                                            labels);
                    }
                    break;
                case 1:
                case 2:
                    if (!generateFundOutlierUniform(w1, h1, w2, h2,
                                                    outlierRatio, noisyInlierMatches.size(), outlierType,
                                                    tempModel, trueModelParams,
                                                    outlierThresholdSq, artificialOutliers, generator, maxIterOutlier)) {
                        std::cerr << "Not enough outliers created: " << artificialOutliers.size() << std::endl;
                    } else {
                        if (verbose) {
                            std::cout << "Number of outliers: " << artificialOutliers.size() << std::endl;
                        }

                        utility::randomDataFusionWithMemory(noisyInlierMatches, artificialOutliers, matchesToEvaluate,
                                                            1, 0,
                                                            labels);
                    }
                    break;
                default:
                    std::cerr << "Model number not supported. Choose in {0, 1, 2}" << std::endl;
                    return 1;
                    break;

            }
        } else {
            std::vector<Match>::const_iterator itNoisyInlierMatches = noisyInlierMatches.begin();
            for (; itNoisyInlierMatches != noisyInlierMatches.end(); itNoisyInlierMatches++) {
                matchesToEvaluate.push_back(*itNoisyInlierMatches);
                labels.push_back(1);
            }
        }

        delete tempModel;

        noisyInliersVect.push_back(noisyInlierMatches);
        artificialOutliersVect.push_back(artificialOutliers);

        std::vector<Match> normalisedPoints;
        orsa::normalisePoints(matchesToEvaluate, intrinsics_source, intrinsics_destination, normalisedPoints);

        if (verbose) {
            std::cout << "\nReady to run experiment: " << gen + 1 << " out of " << nGen << "\n\t"
                      << matchesToEvaluate.size() << " matches.\n\t"
                      << nRun << " Runs.\n"
                      << std::endl;
        }

        // Create the ModelEstimator:
        orsa::ModelEstimator *model;

        switch (modelToAnalyse) {
            case 0:
                std::cout << "Running a Homography estimation..." << std::endl;
                model = new orsa::HomographyModel(matchesToEvaluate, w1, h1, w2, h2, true);
                break;
            case 1:
                std::cout << "Running a Fundamental matrix estimation..." << std::endl;
                model = new orsa::FundamentalModel(matchesToEvaluate, w1, h1, w2, h2, true);
                break;
            case 2:
                model = new orsa::EssentialModel(matchesToEvaluate, normalisedPoints, w1, h1, w2, h2,
                                                 intrinsics_source,
                                                 intrinsics_destination, true);
                break;
            default:
                std::cerr << "Model number not supported. Choose in {0, 1, 2}" << std::endl;
                return 1;
                break;
        }

        // Prepare all algorithms for this artificial dataset:
        orsa::LRTSac lrtAlgorithm = orsa::LRTSac(model);
        lrtAlgorithm.setHyperParameters(cpI, cpIIB, cpIIT, maxSigmaLRT);
        algorithm = &lrtAlgorithm;

        // For each dataset, multiple run are possible:
        for (int run = 0; run < nRun; run++) {
            if (!verbose) {
                std::cout << "\rDataset " << gen + 1 << " out of " << nGen << " - Experiment " << run + 1 << " out of "
                          << nRun << std::flush;
            }

            // Running all algorithms:
            lrtAlgorithm.setCpIIB(cpIIB);
            runExperiment(model, algorithm, labels, noisyInlierMatches, 0,
                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect, algorithmNames,
                          iterMax, verbose);

            lrtAlgorithm.setCpIIB(1.0);
            runExperiment(model, algorithm, labels, noisyInlierMatches, 1,
                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect, algorithmNames,
                          iterMax, verbose);

            if (verbose && nRun > 1) {
                std::cout << std::endl;
            }
        } // End of the 'run' loop.

        delete model;
    } // End of the 'gen' loop.
    if (!verbose) {
        std::cout << std::endl;
    }

    std::cout << std::endl;

    // Computing the precision, recall and runtime mean and std over generated dataset x run
    // per dataset for each algorithm.
    std::vector<double> precisionMeanVect;
    utility::meanOfVect(precisionVectVect, precisionMeanVect);
    std::vector<double> recallMeanVect;
    utility::meanOfVect(recallVectVect, recallMeanVect);
    std::vector<double> runtimeMeanVect;
    utility::meanOfVect(runtimeVectVect, runtimeMeanVect);

    std::vector<double> precisionStdVect;
    utility::standardDeviation(precisionMeanVect, precisionVectVect, precisionStdVect);
    std::vector<double> recallStdVect;
    utility::standardDeviation(recallMeanVect, recallVectVect, recallStdVect);
    std::vector<double> runtimeStdVect;
    utility::standardDeviation(runtimeMeanVect, runtimeVectVect, runtimeStdVect);

    for (size_t i = 0; i < precisionVectVect.size(); i++) {
        switch (i) {
            case 0:
                std::cout << "For LRT with threshold: " << maxSigmaLRT << std::endl;
                break;
            case 1:
                std::cout << "For LRT without early bailout with threshold: " << maxSigmaLRT << std::endl;
                break;
            default:
                break;
        }
        std::cout << "Over " << nGen << " datasets with " << nRun << " runs:\n"
                  << "\tPrecision: " << precisionMeanVect[i] << " +/- " << precisionStdVect[i] << "\n"
                  << "\tRecall: " << recallMeanVect[i] << " +/- " << recallStdVect[i] << "\n"
                  << "\tRuntime: " << runtimeMeanVect[i] << " +/- " << runtimeStdVect[i] << "\n"
                  << std::endl;
    }

    // Saving the experiment info.
    utility::saveExpInfo(pathOutput, seed, nGen, nRun, pathInGoodMatches, algorithmNames,
                         precisionMeanVect, precisionStdVect, precisionVectVect,
                         recallMeanVect, recallStdVect, recallVectVect,
                         runtimeMeanVect, runtimeStdVect, runtimeVectVect,
                         sigmaEstimatedVectVect,
                         true);

    // Saving the artificial dataset.
    utility::saveMatches(pathNoisyInliers, noisyInliersVect);
    utility::saveMatches(pathArtificialOutliers, artificialOutliersVect);

    return 0;
}
