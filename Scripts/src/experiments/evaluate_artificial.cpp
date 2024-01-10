//
// Created by riuclement on 4/7/21.
//

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "libImage/image_io.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/libNumerics/homography.h"

#include "libOrsa/homography_model.hpp"
#include "libOrsa/fundamental_model.hpp"
#include "libOrsa/essential_model.hpp"

#include "libOrsa/lrtsac.hpp"
#include "libOrsa/magsac.hpp"
#include "libOrsa/magsac_metrics.hpp"
#include "libOrsa/muse.hpp"
#include "libOrsa/orsa.hpp"
#include "libOrsa/orsa_fast.hpp"
#include "libOrsa/ransac.hpp"
#include "libOrsa/starsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/siftMatch.hpp"
#include "utilities/metrics.hpp"
#include "utilities/warping.hpp"
#include "utilities/homography_graphical_output.hpp"
#include "utilities/fundamental_graphical_output.hpp"
#include "utilities/usac_file_handler.hpp"
#include "utilities/data_handler.hpp"

// TODO

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

int main(int argc, char **argv) {
    float fSiftRatio = 0.6f; // SIFT parameter. Only used if -r is not used.
    utility::Geometry region;
    region.x0 = region.y0 = region.x1 = region.y1 = 0;
    int nGen = 1; // Number of dataset available in the given matches files. Only the first one will be read.

    int iterMax = 1000; // Maximum number of iterations of the algorithm.
    double precision = 3; // Default pixel precision of the algorithm.
    double cpIIT = 0.99; // Default confidence of the algorithm.

    double precisionRANSAC = precision; // RANSAC specific precision. /// Redundant initialisation for clarity.
    int nModelMinRansac = 5;

    double maxSigmaOrsa = precision; // ORSA specific precision /// Redundant initialisation for clarity.

    double maxSigmaLRT = precision; // LRT specific precision /// Redundant initialisation for clarity.
    double cpI = 0.99;
    double cpIIB = 0.95;

    double maxSigmaMagsac = precision; // Magsac specific precision /// Redundant initialisation for clarity.
    size_t partitionNumber = 5; // Magsac partition number.

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    int modelToAnalyse = 0; // 0 for Homography, 1 for Fundamental, 2 for Essential.
    int algoToUse = 0; // 0 for Ransac, 1 for AC-RANSAC, 2 for LRT.

    utility::CmdLine cmd;
    cmd.add(make_option('c', region, "cut")
                    .doc("cut region of imagInA: wxh+x+y = rect [x,x+w] x [y,y+h]"));
    cmd.add(utility::make_option('s', fSiftRatio, "sift")
                    .doc("SIFT distance ratio of descriptors"));

    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential."));
    cmd.add(utility::make_option('u', algoToUse, "used-algo")
                    .doc("Algorithm to use: 0 for Ransac, 1 for AC-RANSAC, 2 for LRT, 3 for Fast-AC-RANSAC, 4 for StarSaC, 5 for MUSE."));
    cmd.add(utility::make_option('g', nGen, "num-gen-exp")
                    .doc("Number of noisy datasets generated for the experiment."));

    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of the algorithm."));

    cmd.add(utility::make_option('p', precision, "precision")
                    .doc("Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of chosen algorithm."));
    cmd.add(utility::make_option(0, cpIIT, "cpIIT")
                    .doc("Confidence against type II error (Ransac/LRT)"));

    cmd.add(utility::make_option(0, nModelMinRansac, "nModelMinRansac")
                    .doc("Min number of model to evaluate before terminating ransac"));

    cmd.add(utility::make_option(0, cpI, "cpI")
                    .doc("Confidence proba wrt type I error (LRT)"));
    cmd.add(utility::make_option(0, cpIIB, "cpIIB")
                    .doc("Confidence proba wrt bailout (LRT)"));

    cmd.add(utility::make_option(0, partitionNumber, "partition-number")
                    .doc("Number of partition of the Sigma consensus algorithm."));

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

    int minArgc = 6;
    int minArgcE = 7;

    if (((argc < minArgc && argc != minArgc + 2 && argc != minArgc + 3 && argc != minArgc + 5) && modelToAnalyse < 2) ||
        (((argc < minArgcE && argc != minArgcE + 2 && argc != minArgcE + 3) && modelToAnalyse == 2))) {
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
    const char *pathInGoodMatches; // Path to all matches, either as input or as output.
    const char *pathInBadMatches; // Path to all matches, either as input or as output.
    const char *pathOutGoodMatches; // Path to only good matches.

    /// Optionnal arguments :
    const char *pathOutInlierMatchesImage; // All models: image with inlier matches.
    const char *pathOutOutlierMatchesImage; // All models: image with outlier matches.
    const char *pathOutEpiImage; // Fundamental and Essential only: image with epipolar lines.
    const char *pathOutMosaicGlobal; // Homography only: image with reconstructed mosaic.
    const char *pathOutMosaicImg1; // Homography only: warped first image.
    const char *pathOutMosaicImg2; // Homography only: warped second image.
    if (modelToAnalyse < 2) {
        pathInGoodMatches = argv[3];
        pathInBadMatches = argv[4];
        pathOutGoodMatches = argv[5];
        if (argc >= minArgc + 2) {
            pathOutInlierMatchesImage = argv[6];
            pathOutOutlierMatchesImage = argv[7];
        }
        if (argc >= minArgc + 3) {
            pathOutMosaicGlobal = argv[8];
            pathOutEpiImage = argv[8];
        }
        if (argc >= minArgc + 5) {
            pathOutMosaicImg1 = argv[9];
            pathOutMosaicImg2 = argv[10];
        }
    }
    if (modelToAnalyse == 2) {
        pathInCalib = argv[3];
        pathInGoodMatches = argv[4];
        pathInBadMatches = argv[5];
        pathOutGoodMatches = argv[6];
        if (argc >= minArgcE + 2) {
            pathOutInlierMatchesImage = argv[7];
            pathOutOutlierMatchesImage = argv[8];
        }
        if (argc >= minArgcE + 3) {
            pathOutEpiImage = argv[9];
        }
    }

    bool forceCompute = cmd.used('f');
    bool verbose = cmd.used('v');
    // Init random seed
    srand(seed);

    std::ifstream test_f(pathOutGoodMatches);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    precisionRANSAC = precision; // RANSAC specific precision.

    maxSigmaOrsa = precision; // ORSA specific precision

    maxSigmaLRT = precision; // LRT specific precision

    maxSigmaMagsac = precision; // Magsac specific precision

    Image<RGBColor> image1, image2;
    if (!libs::ReadImage(pathInImage1, &image1))
        return 1;
    if (!libs::ReadImage(pathInImage2, &image2))
        return 1;

    bool bUseRegion = cmd.used('c');
    if (bUseRegion) { // Sanity check
        utility::Geometry zone;
        zone.x0 = zone.y0 = 0;
        zone.x1 = int(image1.Width());
        zone.y1 = int(image1.Height());
        if (!(region.x0 < region.x1 && region.y0 < region.y1 &&
              zone.inside(region.x0, region.y0) &&
              zone.inside(region.x1 - 1, region.y1 - 1))) {
            std::cout << "Invalid cut region " << region << std::endl;
            return 1;
        }
    }

    Image<unsigned char> image1Gray, image2Gray;
    libs::convertImage(image1, &image1Gray);
    libs::convertImage(image2, &image2Gray);

    // Find matches with SIFT or read correspondence file
    std::vector<std::vector<Match>> pointsAll; // The point correspondences, each is of format "x1 y1 x2 y2"
    std::vector<std::vector<int>> groundTruthLabelsAll; // The ground truth labeling provided in the dataset
    std::vector<std::vector<Match>> inliersAll; // The point correspondences, each is of format "x1 y1 x2 y2"

    std::cout << "\nReading " << pathInGoodMatches << " and " << pathInBadMatches << std::endl;
    if (!ReadPoints(pathInGoodMatches, pathInBadMatches, nGen, pointsAll, groundTruthLabelsAll, inliersAll, true)) {
        std::cerr << "Problem loading points !" << std::endl;
        return 1;
    }
    std::vector<Match> vec_matchings = pointsAll[0];
    std::vector<int> labels = groundTruthLabelsAll[0];
    std::vector<Match> noisyInlierMatches = inliersAll[0];

    if (bUseRegion) {
        Image<RGBColor> image = image1;
        Crop(image, region.x0, region.y0,
             region.x1 - region.x0, region.y1 - region.y0, image1);
        libs::convertImage(image1, &image1Gray);
    }

    // Remove duplicates (frequent with SIFT)
    utility::rm_duplicates(vec_matchings);

    const int w1 = int(image1Gray.Width()), h1 = int(image1Gray.Height());
    const int w2 = int(image2Gray.Width()), h2 = int(image2Gray.Height());

    // Definition of output values:
    double runtime;

    // Creation of the generic model and generic algorithm:
    orsa::ModelEstimator *model;
    orsa::RansacAlgorithm *algorithm;

    std::vector<Match> vec_matchingsNormalised;
    libNumerics::matrix<double> intrinsics_source(3,3), // The intrinsic parameters of the source camera
    intrinsics_destination(3,3); // The intrinsic parameters of the destination camera

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

    // Definition of all algorithms:
    orsa::Ransac ransacAlgorithm = orsa::Ransac(model);
    ransacAlgorithm.setHyperParameters(precisionRANSAC, cpIIT, nModelMinRansac);

    orsa::Orsa orsaAlgorithm = orsa::Orsa(model);
    orsaAlgorithm.setHyperParameters(maxSigmaOrsa);

    orsa::LRTSac lrtAlgorithm = orsa::LRTSac(model);
    lrtAlgorithm.setHyperParameters(cpI, cpIIB, cpIIT, maxSigmaLRT, true);

    orsa::OrsaFast FastOrsaAlgorithm = orsa::OrsaFast(model);
    FastOrsaAlgorithm.setHyperParameters(maxSigmaOrsa);

    orsa::StarSac StarSacAlgorithm = orsa::StarSac(model);
    StarSacAlgorithm.setHyperParameters(cpIIT, maxSigmaOrsa);

    orsa::Muse museAlgorithm = orsa::Muse(model);
    museAlgorithm.setHyperParameters(cpIIT);

    orsa::MAGSAC magsacAlgorithm = orsa::MAGSAC(model);
    magsacAlgorithm.setHyperParameters(cpIIT, maxSigmaMagsac, partitionNumber);

    orsa::MAGSAC magsacplusplusAlgorithm = orsa::MAGSAC(model, orsa::MAGSAC::Version::MAGSAC_PLUS_PLUS);
    magsacplusplusAlgorithm.setHyperParameters(cpIIT, maxSigmaMagsac, partitionNumber);

    switch (algoToUse) {
        case 0:
            std::cout << "... with algorithm Ransac." << std::endl;
            algorithm = &ransacAlgorithm;
            break;
        case 1:
            std::cout << "... with algorithm AC-Ransac." << std::endl;
            algorithm = &orsaAlgorithm;
            break;
        case 2:
            std::cout << "... with algorithm LRTSac." << std::endl;
            algorithm = &lrtAlgorithm;
            break;
        case 3:
            std::cout << "... with algorithm Fast-AC-Ransac." << std::endl;
            algorithm = &FastOrsaAlgorithm;
            break;
        case 4:
            std::cout << "... with algorithm StarSac." << std::endl;
            algorithm = &StarSacAlgorithm;
            break;
        case 5:
            std::cout << "... with algorithm MUSE." << std::endl;
            algorithm = &museAlgorithm;
            break;
        case 6:
            std::cout << "... with algorithm MAGSAC." << std::endl;
            algorithm = &magsacAlgorithm;
            break;
        case 7:
            std::cout << "... with algorithm PlusPlusMAGSAC." << std::endl;
            algorithm = &magsacplusplusAlgorithm;
            break;
        default:
            std::cerr << "Algorithm number not supported. Choose in integer between 0 and 6" << std::endl;
            delete model;
            return 1;
            break;
    }

    // Estimation of the model parameters:
    orsa::RansacAlgorithm::RunResult res;
    bool ok = algorithm->evalModel(res, runtime, iterMax, verbose);

    double error = -1;
    int numTruePositives = 0;
    double computedPrecision = 0;
    double computedRecall = 0;
    int inlierCount = res.vInliers.size();

    // Evaluate the estimated model:
    if (ok) {
        error = model->ErrorStats(res.vInliers, res.model).first;
        numTruePositives = utility::computeTruePositive(res.vInliers, labels);
        computedPrecision = utility::computePrecision(numTruePositives, res.vInliers.size());
        computedRecall = utility::computeRecall(numTruePositives, noisyInlierMatches.size());
    }

    if (ok && verbose) {
        std::cout << "\tModelParam=" << res.model << std::endl;
        std::cout << "\tInlier error: " << error << "\n\tComputed sigma: " << res.sigma
                  << "\n\tNumber of inliers: " << inlierCount
                  << "\n\tNumber of iterations: " << res.T << "\n\tVPM: " << res.vpm
                  << "\n\tRuntime: " << runtime
                  << std::endl;
        std::cout << "\tNumber true positives: " << numTruePositives
                  << "\n\tPrecision: " << computedPrecision
                  << "\n\tRecall: " << computedRecall
                  << std::endl;
    }

    std::vector<Match> good_match;
    std::vector<int>::const_iterator it = res.vInliers.begin();
    for (; it != res.vInliers.end(); it++)
        good_match.push_back(vec_matchings[*it]);

    if (!Match::saveMatch(pathOutGoodMatches, good_match)) {
        std::cerr << "Failed saving matches into " << pathOutGoodMatches << std::endl;
        delete model;
        return 1;
    }
    if (!ok) {
        std::cerr << "Failed to estimate a model" << std::endl;
        delete model;
        return 0;
    }

    // Output images:
    switch (modelToAnalyse) {
        case 0:
            // Sift de-duplicated output display
            if (argc >= minArgc + 2)
                utility::homography_matches_output(image1Gray, image2Gray,
                                                   vec_matchings, res.vInliers,
                                                   ok ? &res.model : 0,
                                                   pathOutInlierMatchesImage, pathOutOutlierMatchesImage);
            // Mosaics
            if (argc >= minArgc + 3) {
                utility::Rect intersection;
                if (utility::IntersectionBox(w1, h1, w2, h2, res.model, intersection) &&
                    intersection.Width() > 0 && intersection.Height() > 0) {
                    const char *fileReg1 = 0, *fileReg2 = 0;
                    if (argc >= minArgc + 5) {
                        fileReg1 = pathOutMosaicImg1;
                        fileReg2 = pathOutMosaicImg2;
                    }
                    utility::homography_registration_output(image1, image2, res.model, intersection,
                                                            pathOutMosaicGlobal, fileReg1, fileReg2, true);
                }
            }
            break;
        case 1:
            if (argc >= minArgc + 2) // Output images
            {
                const char *fileIn = 0, *fileOut = 0, *fileEpi = 0;
                fileIn = pathOutInlierMatchesImage;
                fileOut = pathOutOutlierMatchesImage;
                if (argc >= minArgc + 3)
                    fileEpi = pathOutEpiImage;
                utility::fundamental_graphical_output(image1Gray, image2Gray,
                                                      vec_matchings, res.vInliers, ok ? &res.model : 0,
                                                      fileIn, fileOut, fileEpi);
            }
            break;
        case 2:
            if (argc >= minArgcE + 2) {// Output images
                libNumerics::matrix<double> modelParamF = model->toPixelSpace(res.model);
                const char *fileIn = 0, *fileOut = 0, *fileEpi = 0;
                fileIn = pathOutInlierMatchesImage;
                fileOut = pathOutOutlierMatchesImage;
                if (argc >= minArgcE + 3)
                    fileEpi = pathOutEpiImage;
                utility::fundamental_graphical_output(image1Gray, image2Gray,
                                                      vec_matchings, res.vInliers, &modelParamF,
                                                      fileIn, fileOut, fileEpi);
            }
            break;
        default:
            std::cerr << "Model number not supported. Choose in {0, 1, 2}" << std::endl;
            delete model;
            return 1;
            break;
    }
    delete model;

    return 0;
}
