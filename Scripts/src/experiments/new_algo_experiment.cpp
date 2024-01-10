//
// Created by clementriu on 2/24/21.
//

#include <cstdlib>
#include <ctime>
#include <iostream>

#include "libOrsa/homography_model.hpp"
#include "libOrsa/fundamental_model.hpp"
#include "libOrsa/essential_model.hpp"

#include "libOrsa/muse.hpp"
#include "libOrsa/orsa_fast.hpp"
#include "libOrsa/magsac.hpp"
#include "libOrsa/starsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/data_handler.hpp"
#include "utilities/metrics.hpp"
#include "utilities/usac_file_handler.hpp"

#include "read_write_experiment.hpp"


bool ReadPoints(const char *fileInliers, const char *fileOutliers, const int nGen,
                std::vector<std::vector<Match>> &pointsAll, std::vector<std::vector<int>> &groundTruthLabelsAll,
                std::vector<std::vector<Match>> &inliersAll, const bool readOutliers) {
    std::vector<Match> ptsInliersAll;
    std::vector<Match> ptsOutliersAll;

    if (!Match::loadMatch(fileInliers, ptsInliersAll)) {
        std::cerr << "Problem loading inliers: " << fileInliers << std::endl;
        return false;
    }
    if (readOutliers) {
        if (!Match::loadMatch(fileOutliers, ptsOutliersAll)) {
            std::cerr << "Problem loading outliers: " << fileOutliers << std::endl;
            return false;
        }
    }
    const int numPointsIn = ptsInliersAll.size();
    const int numPointsOut = ptsOutliersAll.size();

    assert(numPointsIn % nGen == 0);
    assert(numPointsOut % nGen == 0);
    const int numPointsInPerGen = numPointsIn / nGen;
    const int numPointsOutPerGen = numPointsOut / nGen;

    for (int gen = 0; gen < nGen; gen++) {
        std::vector<Match> ptsInliers;
        std::vector<Match> ptsOutliers;

        for (int i = gen * numPointsInPerGen; i < (gen + 1) * numPointsInPerGen; i++) {
            ptsInliers.push_back(ptsInliersAll[i]);
        }
        if (readOutliers) {
            for (int i = gen * numPointsOutPerGen; i < (gen + 1) * numPointsOutPerGen; i++) {
                ptsOutliers.push_back(ptsOutliersAll[i]);
            }
        }

        std::vector<Match> ptsMixed;
        std::vector<int> groundTruthLabels;

        utility::randomDataFusionWithMemory(ptsInliers, ptsOutliers, ptsMixed, 1, 0, groundTruthLabels);

        pointsAll.push_back(ptsMixed);
        groundTruthLabelsAll.push_back(groundTruthLabels);
        inliersAll.push_back(ptsInliers);
    }
    return true;
}



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
bool runExperiment(
        const orsa::ModelEstimator *model, const orsa::RansacAlgorithm *algorithm,
        const std::vector<int> &labels, const std::vector<Match> &noisyInlierMatches,
        const int algoNumber,
        std::vector<std::vector<double>> &sigmaEstimatedVectVect,
        std::vector<std::vector<double>> &runtimeVectVect,
        std::vector<std::vector<double>> &precisionVectVect,
        std::vector<std::vector<double>> &recallVectVect,
        const std::vector<const char *> &algorithmNames,
        const int iterMax, const bool verbose) {
    // Output values:
    double error = 0;
    orsa::RansacAlgorithm::RunResult res;
    double runtime;
    // Run the given model with given algorithm:
    bool ok = algorithm->evalModel(res, runtime, iterMax, verbose);

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
bool runExperiment(
        const orsa::ModelEstimator *model,
        const orsa::RansacAlgorithm *algorithm, const orsa::MAGSAC &magsacAlgorithm,
        const std::vector<int> &labels, const std::vector<Match> &noisyInlierMatches,
        const int algoNumber,
        std::vector<std::vector<double>> &sigmaEstimatedVectVect,
        std::vector<std::vector<double>> &runtimeVectVect,
        std::vector<std::vector<double>> &precisionVectVect,
        std::vector<std::vector<double>> &recallVectVect,
        std::vector<std::vector<int>> &possibleInliersVect,
        std::vector<std::vector<double>> &weightsVect,
        std::vector<std::vector<double>> &errorsVect,
        std::vector<std::vector<double>> &errorsAllVect,
        const std::vector<const char *> &algorithmNames,
        const int iterMax, const bool verbose) {
    // Output values:
    double error = 0;
    orsa::RansacAlgorithm::RunResult res;
    double runtime;
    // Run the given model with given algorithm:
    bool ok = algorithm->evalModel(res, runtime, iterMax, verbose);

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

    if (ok) {
        possibleInliersVect.push_back(res.vInliers);
        weightsVect.push_back(magsacAlgorithm.getWeights());
        std::vector<double> error;
        utility::computeModelError(res.vInliers, model, res.model, error);
        errorsVect.push_back(error);
        error.clear();
        utility::computeModelError(model, res.model, error);
        errorsAllVect.push_back(error);
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
    int iterMax = 10000; // Maximum number of iterations of the algorithm.
    int iterMaxStarSac = iterMax; // Maximum number of iterations of the algorithm.
    double maxPrecision = 16; // Default pixel precision of the algorithm.
    double cpII = 0.99; // Default confidence of the algorithm.

    double maxSigmaOrsa = maxPrecision; // ORSA specific precision /// Redundant initialisation for clarity.

    double maxSigmaStarSac = maxPrecision; // LRT specific precision /// Redundant initialisation for clarity.

    int nGen = 1; // Number of different datasets to generate.
    int nRun = 1; // Number of run for each different datasets.

    int modelToAnalyse = 0; // 0 for Homography, 1 for Fundamental, 2 for Essential.

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    utility::CmdLine cmd;
    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of algorithms."));
    cmd.add(utility::make_option('s', iterMaxStarSac, "iterMaxStarSac")
                    .doc("Number of iterations of algorithm starsac."));

    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential."));

    cmd.add(utility::make_option(0, maxPrecision, "maximum-precision")
                    .doc("Max precision (in pixels) of registration of Fast-AC-RANSAC and StarSac."));
    cmd.add(utility::make_option(0, cpII, "cpII")
                    .doc("Value of confidence (Fast-AC-RANSAC and StarSac)"));

    cmd.add(utility::make_option('g', nGen, "num-gen-exp")
                    .doc("Number of noisy datasets generated for the experiment."));
    cmd.add(utility::make_option('n', nRun, "num-run")
                    .doc("Number of run of the algorithm."));
    cmd.add(utility::make_switch('o', "ignore-outliers")
                    .doc("This dataset does not contain outliers"));

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

    int minArgcHF = 10;
    int minArgcE = 11;

    if ((modelToAnalyse < 2 && argc != minArgcHF) || (modelToAnalyse == 2 && argc != minArgcE)) {
        std::cerr << "Usage: " << argv[0]
                  << " [options] good_matches.txt info.txt calib.txt outputInfo.txt\n"
                  << "- good_matches.txt: path to the good matches to read.\n"
                  << "- info.txt: path to the run info (like dimensions of images).\n"
                  << "- calib.txt: path to the calib matrix of both images.\n"
                  << "- outputInfo.txt: path to the output file with all the informations on the experiment.\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }

    bool forceCompute = cmd.used('f');
    bool verbose = cmd.used('v');
    bool readOutliers = !cmd.used('o');
    // Init random seed
    std::srand(seed);

    int indexArg = 1;
    const char *pathInInfo = argv[indexArg++]; // Output info file of a generate_artificial_dataset run.
    const char *pathInCalib; // Essential only: calibration matrix of both images.
    const char *pathInGoodMatches; // Path to all matches, either as input or as output.
    const char *pathInBadMatches; // Path to all matches, either as input or as output.
    const char *pathOutput; // Path to only good matches.

    // MAGSAC specific paths:
    const char *pathToOutLabels;
    const char *pathToOutWeights;
    const char *pathToOutPossibleInliers;
    const char *pathToOutErrors;
    const char *pathToOutErrorsAll;

    pathInGoodMatches = argv[indexArg++];
    pathInBadMatches = argv[indexArg++];
    pathOutput = argv[indexArg++];
    if (modelToAnalyse == 2) {
        pathInCalib = argv[indexArg++];
    }

    pathToOutLabels = argv[indexArg++];
    pathToOutWeights = argv[indexArg++];
    pathToOutPossibleInliers = argv[indexArg++];
    pathToOutErrors = argv[indexArg++];
    pathToOutErrorsAll = argv[indexArg++];


    std::ifstream test_f(pathOutput);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    maxSigmaOrsa = maxPrecision; // ORSA specific precision
    maxSigmaStarSac = maxPrecision; // LRT specific precision
    if (!cmd.used('s')) {
        iterMaxStarSac = iterMax;
    }

    char FastACRansacText[100] = "Fast-AC-Ransac with threshold ";
    char StarSacText[100] = "StarSac with threshold ";
    char MuseText[100] = "MUSE";
    char MagsacText[100] = "MAGSAC";
    std::vector<const char *> algorithmNames = {
            concatenateTextValue(FastACRansacText, maxSigmaOrsa),
//            concatenateTextValue(StarSacText, maxSigmaStarSac),
            MuseText,
            MagsacText
    };

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

    // Find matches with SIFT or read correspondence file
    std::vector<std::vector<Match>> pointsAll; // The point correspondences, each is of format "x1 y1 x2 y2"
    std::vector<std::vector<int>> groundTruthLabelsAll; // The ground truth labeling provided in the dataset
    std::vector<std::vector<Match>> inliersAll; // The point correspondences, each is of format "x1 y1 x2 y2"

    std::cout << "\nReading " << pathInGoodMatches << " and " << pathInBadMatches << std::endl;
    if (!ReadPoints(pathInGoodMatches, pathInBadMatches, nGen, pointsAll, groundTruthLabelsAll, inliersAll,
                    readOutliers)) {
        std::cerr << "Problem loading points !" << std::endl;
        return 1;
    }

    libNumerics::matrix<double> intrinsics_source(3, 3), // The intrinsic parameters of the source camera
    intrinsics_destination(3, 3); // The intrinsic parameters of the destination camera
    if (modelToAnalyse == 2) {
        if (utility::loadCalibration(pathInCalib, intrinsics_source, intrinsics_destination))
            std::cout << "Read " << pathInCalib << " calibration file." << std::endl;
        else {
            std::cerr << "Failed reading calib from " << pathInCalib << std::endl;
            return 1;
        }
    }
    // Generic algorithm pointer that will be used to call each algorithm.
    orsa::RansacAlgorithm *algorithm;

    // All metrics are stored across each algorithm across each run.
    int numberAlgorithm = algorithmNames.size();
    std::vector<std::vector<double>> precisionVectVect(numberAlgorithm);
    std::vector<std::vector<double>> recallVectVect(numberAlgorithm);
    std::vector<std::vector<double>> runtimeVectVect(numberAlgorithm);
    std::vector<std::vector<double>> sigmaEstimatedVectVect(numberAlgorithm);

    // MAGSAC specific outputs:
    std::vector<std::vector<int>> possibleInliersVect;
    std::vector<std::vector<double>> weightsVect;
    std::vector<std::vector<double>> errorsVect;
    std::vector<std::vector<double>> errorsAllVect;

    // Generate multiple semi-artificial dataset.
    for (int gen = 0; gen < nGen; gen++) {
        std::vector<Match> matchesToEvaluate = pointsAll[gen];
        std::vector<int> labels = groundTruthLabelsAll[gen];
        std::vector<Match> noisyInlierMatches = inliersAll[gen];

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
                if (verbose) std::cout << "Running a Homography estimation..." << std::endl;
                model = new orsa::HomographyModel(matchesToEvaluate, w1, h1, w2, h2, true);
                break;
            case 1:
                if (verbose) std::cout << "Running a Fundamental matrix estimation..." << std::endl;
                model = new orsa::FundamentalModel(matchesToEvaluate, w1, h1, w2, h2, true);
                break;
            case 2:
                if (verbose) std::cout << "Running a Essential matrix estimation..." << std::endl;
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
//        orsa::Orsa OrsaAlgorithm = orsa::Orsa(model);
//        OrsaAlgorithm.setHyperParameters(maxSigmaOrsa);
        orsa::OrsaFast FastOrsaAlgorithm = orsa::OrsaFast(model);
        FastOrsaAlgorithm.setHyperParameters(maxSigmaOrsa);
//        orsa::StarSac StarSacAlgorithm = orsa::StarSac(model);
//        StarSacAlgorithm.setHyperParameters(cpII, maxSigmaStarSac);
        orsa::Muse museAlgorithm = orsa::Muse(model);
        museAlgorithm.setHyperParameters(cpII);
        orsa::MAGSAC magsacAlgorithm = orsa::MAGSAC(model);
        magsacAlgorithm.setHyperParameters(cpII);
        orsa::MAGSAC magsacplusplusAlgorithm = orsa::MAGSAC(model, orsa::MAGSAC::Version::MAGSAC_PLUS_PLUS);
        magsacplusplusAlgorithm.setHyperParameters(cpII);

        // For each dataset, multiple run are possible:
        for (int run = 0; run < nRun; run++) {
            if (!verbose) {
                std::cout << "\rDataset " << gen + 1 << " out of " << nGen << " - Experiment " << run + 1 << " out of "
                          << nRun << std::flush;
            }

//            // Running all algorithms:
            int algoNumber = 0;
            algorithm = &FastOrsaAlgorithm;
            runExperiment(model, algorithm, labels, noisyInlierMatches, algoNumber++,
                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect, algorithmNames,
                          iterMax, verbose);

//            algorithm = &StarSacAlgorithm;
//            runExperiment(model, algorithm, labels, noisyInlierMatches, algoNumber++,
//                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect, algorithmNames,
//                          iterMax, verbose);

            algorithm = &museAlgorithm;
            runExperiment(model, algorithm, labels, noisyInlierMatches, algoNumber++,
                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect, algorithmNames,
                          iterMax, verbose);

            algorithm = &magsacAlgorithm;
            runExperiment(model, algorithm, magsacAlgorithm, labels, noisyInlierMatches, algoNumber++,
                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect,
                          possibleInliersVect, weightsVect, errorsVect, errorsAllVect,
                          algorithmNames, iterMax, verbose);

            algorithm = &magsacplusplusAlgorithm;
            runExperiment(model, algorithm, magsacplusplusAlgorithm, labels, noisyInlierMatches, algoNumber++,
                          sigmaEstimatedVectVect, runtimeVectVect, precisionVectVect, recallVectVect,
                          possibleInliersVect, weightsVect, errorsVect, errorsAllVect,
                          algorithmNames, iterMax, verbose);

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
        std::cout << "For " << algorithmNames[i] << ": " << std::endl;
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

    saveVectOfVect(pathToOutLabels, groundTruthLabelsAll);
    saveVectOfVect(pathToOutWeights, weightsVect);
    saveVectOfVect(pathToOutPossibleInliers, possibleInliersVect);
    saveVectOfVect(pathToOutErrors, errorsVect);
    saveVectOfVect(pathToOutErrorsAll, errorsAllVect);

    return 0;
}
