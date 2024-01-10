//
// Created by riuclement on 3/17/21.
//

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "libImage/image_io.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/model_estimator.hpp"
#include "libOrsa/homography_model.hpp"
#include "libOrsa/fundamental_model.hpp"
#include "libOrsa/essential_model.hpp"

#include "libOrsa/lrtsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/siftMatch.hpp"
#include "utilities/usac_file_handler.hpp"
#include "utilities/data_handler.hpp"

void runExperiment(const int nRun, const int iterMax, const bool verbose, const orsa::LRTSac &lrtAlgorithm,
                   const int algoNumber,
                   std::vector<std::vector<int>> &TVectVect,
                   std::vector<std::vector<double>> &vpmVectVect,
                   std::vector<double> &runtimeVect) {
    clock_t begin = std::clock();
    for (int i = 0; i < nRun; i++) {
        if (!verbose) {
            std::cout << "Run " << i + 1 << " out of " << nRun << std::flush;
        }
        orsa::RansacAlgorithm::RunResult res;
        lrtAlgorithm.run(res, iterMax, verbose);
        TVectVect[algoNumber].push_back(res.T);
        vpmVectVect[algoNumber].push_back(res.vpm);
    }
    if (!verbose) {
        std::cout << std::endl;
    }
    clock_t end = std::clock();
    runtimeVect.push_back((double) (end - begin) / CLOCKS_PER_SEC);
    std::cout << "Done." << std::endl;
}

int main(int argc, char **argv) {
    int iterMax = 1000; // Maximum number of iterations of the algorithm.
    int nRun = 100; // Number of run over which to average the runtime, vpm and number of iterations.

    double precision = 16; // Default pixel precision of the algorithm.
    double cpIIT = 0.99; // Default confidence wrt type II error (iterations)

    double maxSigmaLRT = precision; // LRT specific precision /// Redundant initialisation for clarity.
    double cpI = 0.99;
    double cpIIB = 0.95;

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    int modelToAnalyse = 0; // 0 for Homography, 1 for Fundamental, 2 for Essential.

    utility::CmdLine cmd;
    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential."));

    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of the algorithm."));
    cmd.add(utility::make_option('n', nRun, "num-run")
                    .doc("Number of run of the algorithm."));

    cmd.add(utility::make_option('p', precision, "precision")
                    .doc("Max precision (in pixels) of registration of LRT."));
    cmd.add(utility::make_option(0, cpIIT, "cpIIT")
                    .doc("Value of confidence proba for LRT"));

    cmd.add(utility::make_option(0, cpI, "cpI")
                    .doc("Confidence proba wrt type I error (LRT)"));
    cmd.add(utility::make_option(0, cpIIB, "cpIIB")
                    .doc("Confidence proba wrt bailout (LRT)"));

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

    int minArgc = 5;
    int minArgcE = 6;

    if (((argc < minArgc) && modelToAnalyse < 2) ||
        (((argc < minArgcE) && modelToAnalyse == 2))) {
        std::cerr << "Usage: " << argv[0] << " [options] imgInA imgInB "
                  << "[calibrationMatrix] "
                  << "allInOutMatches.txt inlierOutMatches.txt "
                  << "[imgOutInliers imgOutOutliers [imgOutMosaic/imgOutEpi "
                  << "[imgOutMosaicA imgOutMosaicA]]\n"
                  << "- imgInA, imgInB: the two input images (JPG or PNG format)\n"
                  << "- calibrationMatrix (only if Essential): one file containing the two calibration matrixes.\n"
                  << "- allInOutMatches.txt: output (input if option -r) text file of "
                     "format \"x1 y1 x2 y2\"\n"
                  << "- inlierOutMatches.txt: output, but only with inliers.\n"
                  << "- imgOutInliers (optional): output image showing inliers\n"
                  << "- imgOutOutliers (optional): output image showing outliers and "
                     "their error\n"
                  << "- imgOutMosaic/imgOutEpi (optional): output mosaic image if homography, output epipolar image if Fundamental or Essential\n"
                  << "- imgOutMosaicA, imgOutMosaicA (optional): registered images if homography\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }

    const char *pathInImage1 = argv[1]; // First image path.
    const char *pathInImage2 = argv[2]; // Second image path.
    const char *pathInCalib; // Essential only: calibration matrix of both images.
    const char *pathInMatches; // Path to all matches as input.
    const char *pathOutMetrics; // Path to only good matches.

    if (modelToAnalyse < 2) {
        pathInMatches = argv[3];
        pathOutMetrics = argv[4];
    }
    if (modelToAnalyse == 2) {
        pathInCalib = argv[3];
        pathInMatches = argv[4];
        pathOutMetrics = argv[5];
    }

    bool forceCompute = cmd.used('f');
    bool verbose = cmd.used('v');
    // Init random seed
    srand(seed);

    std::ifstream test_f(pathOutMetrics);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    maxSigmaLRT = precision; // LRT specific precision

    Image<RGBColor> image1, image2;
    if (!libs::ReadImage(pathInImage1, &image1))
        return 1;
    if (!libs::ReadImage(pathInImage2, &image2))
        return 1;

    Image<unsigned char> image1Gray, image2Gray;
    libs::convertImage(image1, &image1Gray);
    libs::convertImage(image2, &image2Gray);

    // Read correspondence file
    std::vector<Match> vec_matchings;

    if (Match::loadMatch(pathInMatches, vec_matchings))
        std::cout << "Read " << vec_matchings.size() << " matches" << std::endl;
    else {
        std::cerr << "Failed reading matches from " << pathInMatches << std::endl;
        return 1;
    }


    // Remove duplicates (frequent with SIFT)
    utility::rm_duplicates(vec_matchings);

    const int w1 = int(image1Gray.Width()), h1 = int(image1Gray.Height());
    const int w2 = int(image2Gray.Width()), h2 = int(image2Gray.Height());

    // Creation of the generic model and generic algorithm:
    orsa::ModelEstimator *model;

    std::vector<Match> vec_matchingsNormalised;
    libNumerics::matrix<double> intrinsics_source(3, 3), // The intrinsic parameters of the source camera
    intrinsics_destination(3, 3); // The intrinsic parameters of the destination camera

    // Creation of the model:
    switch (modelToAnalyse) {
        case 0:
            std::cout << "Running a Homography estimation..." << std::endl;
            model = new orsa::HomographyModel(vec_matchings, w1, h1, w2, h2, true);
            break;
        case 1:
            std::cout << "Running a Fundamental matrix estimation..." << std::endl;
            model = new orsa::FundamentalModel(vec_matchings, w1, h1, w2, h2, true);
            break;
        case 2:
            std::cout << "Running a Essential matrix estimation..." << std::endl;
            if (utility::loadCalibration(pathInCalib, intrinsics_source, intrinsics_destination))
                std::cout << "Read " << pathInCalib << " calibration" << std::endl;
            else {
                std::cerr << "Failed reading calibration matrixes from " << pathInCalib << std::endl;
                return 1;
            }

            orsa::normalisePoints(vec_matchings, intrinsics_source, intrinsics_destination, vec_matchingsNormalised);

            model = new orsa::EssentialModel(vec_matchings, vec_matchingsNormalised, w1, h1, w2, h2, intrinsics_source,
                                             intrinsics_destination, true);
            break;
        default:
            std::cerr << "Model number not supported. Choose in {0, 1, 2}" << std::endl;
            return 1;
            break;
    }

    if (model->NbData() < model->SizeSample()) {
        std::cerr << "Error: The algorithm needs " << model->SizeSample()
                  << " matches or more to proceed" << std::endl;
        delete model;
        return 1;
    }

    std::vector<const char *> algorithmNames = {
            "LRT with bailout, with update of T, with reduce Sigma",          // 000
            "LRT with bailout, without update of T, with reduce Sigma",       // 010
            "LRT without bailout, with update of T, with reduce Sigma",       // 100
            "LRT without bailout, without update of T, with reduce Sigma",    // 110
            "LRT without bailout, without update of T, without reduce Sigma", // 111
    };

    int numberAlgorithms = algorithmNames.size();

    std::vector<std::vector<int>> TVectVect(numberAlgorithms);
    std::vector<std::vector<double>> vpmVectVect(numberAlgorithms);
    std::vector<double> runtimeVect;

    // Definition of all algorithms:
    orsa::LRTSac lrtAlgorithm = orsa::LRTSac(model);
    std::cout << "... with algorithm LRTSac." << std::endl;
    int indexAlgorithm = 0;
    lrtAlgorithm.setHyperParameters(cpI, cpIIB, cpIIT, maxSigmaLRT, true);    // 000
    runExperiment(nRun, iterMax, verbose, lrtAlgorithm, indexAlgorithm++, TVectVect, vpmVectVect, runtimeVect);
    lrtAlgorithm.setHyperParameters(cpI, cpIIB, 1.0, maxSigmaLRT, true);   // 010
    runExperiment(nRun, iterMax, verbose, lrtAlgorithm, indexAlgorithm++, TVectVect, vpmVectVect, runtimeVect);
    lrtAlgorithm.setHyperParameters(cpI, 1.0, cpIIT, maxSigmaLRT, true);   // 100
    runExperiment(nRun, iterMax, verbose, lrtAlgorithm, indexAlgorithm++, TVectVect, vpmVectVect, runtimeVect);
    lrtAlgorithm.setHyperParameters(cpI, 1.0, 1.0, maxSigmaLRT, true);  // 110
    runExperiment(nRun, iterMax, verbose, lrtAlgorithm, indexAlgorithm++, TVectVect, vpmVectVect, runtimeVect);
    lrtAlgorithm.setHyperParameters(cpI, 1.0, 1.0, maxSigmaLRT, false); // 111
    runExperiment(nRun, iterMax, verbose, lrtAlgorithm, indexAlgorithm++, TVectVect, vpmVectVect, runtimeVect);

    // Analysis of results.
    std::vector<double> TMeanVect;
    utility::meanOfVect(TVectVect, TMeanVect);
    std::vector<double> vpmMeanVect;
    utility::meanOfVect(vpmVectVect, vpmMeanVect);

    std::vector<double> TStdVect;
    utility::standardDeviation(TMeanVect, TVectVect, TStdVect);
    std::vector<double> vpmStdVect;
    utility::standardDeviation(vpmMeanVect, vpmVectVect, vpmStdVect);

    std::vector<double> meanRuntime;
    std::vector<double>::const_iterator itRuntimeVect = runtimeVect.begin();
    for (; itRuntimeVect != runtimeVect.end(); itRuntimeVect++) {
        meanRuntime.push_back(*itRuntimeVect / nRun);
    }

    for (size_t i = 0; i < TMeanVect.size(); i++) {
        std::cout << algorithmNames[i] << "\n";
        std::cout << "Over " << nRun << " iterations:\n"
                  << "\tNumber of iterations: " << TMeanVect[i] << " +/- " << TStdVect[i] << "\n"
                  << "\tVPM: " << vpmMeanVect[i] << " +/- " << vpmStdVect[i] << "\n"
                  << "\tRuntime: " << meanRuntime[i] << "s\n"
                  << std::endl;
    }

    utility::saveExpInfo(pathOutMetrics, seed, nRun, pathInMatches, algorithmNames,
                         TMeanVect, TStdVect, TVectVect,
                         vpmMeanVect, vpmStdVect, vpmVectVect,
                         meanRuntime,
                         true);

    delete model;

    return 0;
}
