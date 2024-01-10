//
// Created by clementriu on 2/17/21.
//

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <fstream>

#include "libImage/image_io.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/libNumerics/homography.h"

#include "libOrsa/homography_model.hpp"
#include "libOrsa/fundamental_model.hpp"
#include "libOrsa/essential_model.hpp"

#include "libOrsa/orsa.hpp"

#include "utilities/cmdLine.h"
#include "utilities/Rect.hpp"
#include "utilities/siftMatch.hpp"
#include "utilities/metrics.hpp"
#include "utilities/warping.hpp"
#include "utilities/homography_graphical_output.hpp"
#include "utilities/fundamental_graphical_output.hpp"
#include "utilities/usac_file_handler.hpp"

#include "apply_transform.hpp"
#include "read_write_experiment.hpp"


/// File to select a set of inliers and correct them with a model estimated with AC-RANSAC.

int main(int argc, char **argv) {
    int iterMax = 10000; // Maximum number of iterations of the algorithm.
    double precision = 0; // Default pixel precision of the algorithm.

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    int modelToAnalyse = 0; // 0 for Homography, 1 for Fundamental, 2 for Essential.

    utility::CmdLine cmd;
    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential."));

    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of AC-RANSAC."));

    cmd.add(utility::make_option('p', precision, "precision")
                    .doc("Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of RANSAC, AC-RANSAC and LRT."));

    cmd.add(utility::make_switch('g', "generate-sift")
                    .doc("Generate SIFT matches instead of reading match file."));
    cmd.add(utility::make_switch('f', "force-compute")
                    .doc("Run the computation even if an output file already exists."));
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

    bool generateSift = cmd.used('g');
    bool forceCompute = cmd.used('f');
    bool verbose = cmd.used('v');

    int modifier = 0;
    if (generateSift) {
        modifier = -1;
    }

    int minArgcH = 8 + modifier;
    int minArgcF = 8 + modifier;
    int minArgcE = 9 + modifier;

    if (
            (modelToAnalyse == 0 &&
             !(argc == minArgcH || argc == minArgcH + 1))
            || (modelToAnalyse == 1 &&
                !(argc == minArgcF || argc == minArgcF + 2 || argc == minArgcF + 3))
            || (modelToAnalyse == 2 &&
                !(argc == minArgcE || argc == minArgcE + 2 || argc == minArgcE + 3))
            ) {
        std::cerr << "Usage: " << argv[0] << " [options] [inMatches] imgInA imgInB "
                  << "[calibrationMatrix] "
                  << "OutGoodMatches.txt outInfo.txt "
                  << "[imgWarped]/[imgOutInliers imgOutOutliers [imgOutEpi]] "
                  << "\n"
                  << "- inMatches: file containing the matches to read.\n"
                  << "- imgInA, imgInB: the two input images (JPG or PNG format)\n"
                  << "- calibrationMatrix (only if Essential): one file containing the two calibration matrixes.\n"
                  << "- OutGoodMatches.txt: output good matches text file of "
                     "format \"x1 y1 x2 y2\"\n"
                  << "- outInfo.txt: info about the run.\n"
                  << "- imgOutInliers (optional): output image showing inliers if Fundamental or Essential\n"
                  << "- imgOutOutliers (optional): output image showing outliers and "
                     "their error if Fundamental or Essential\n"
                  << "- imgOutEpi (optional): output epipolar image if Fundamental or Essential\n"
                  << "- imgWarped (optional): registered image if homography\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }

    int readIndex = 1;
    const char *pathInMatches = argv[readIndex++]; // USAC format in file containing matches.
    const char *pathInImage1 = argv[readIndex++ + modifier]; // First image.
    const char *pathInImage2 = argv[readIndex++ + modifier]; // Second image.
    const char *pathOutInlierError = argv[readIndex++ + modifier];
    const char *pathOutOutlierError = argv[readIndex++ + modifier];
    const char *pathInCalib; // Essential only: calibration matrix of both images.
    const char *pathOutGoodMatches; // Output file with good matches coordinate.
    const char *pathOutInfo; // Output file with info about the run.


    /// Optionnal arguments :
    const char *pathOutInlierMatchesImage; // Fundamental and Essential only: image with inlier matches.
    const char *pathOutOutlierMatchesImage; // Fundamental and Essential only: image with outlier matches.
    const char *pathOutEpiImage; // Fundamental and Essential only: image with epipolar lines.
    const char *pathOutWarpedImage; // Homography only: image with wraped image.

    if (modelToAnalyse < 2) {
        pathOutGoodMatches = argv[readIndex++ + modifier];
        pathOutInfo = argv[readIndex++ + modifier];
        if (modelToAnalyse == 0 && argc == minArgcH + 1) {
            pathOutWarpedImage = argv[readIndex++ + modifier];
        }
        if (modelToAnalyse == 1 && argc >= minArgcF + 2) {
            pathOutInlierMatchesImage = argv[readIndex++ + modifier];
            pathOutOutlierMatchesImage = argv[readIndex++ + modifier];
            if (argc == minArgcF + 3) {
                pathOutEpiImage = argv[readIndex++ + modifier];
            }
        }
    }
    if (modelToAnalyse == 2) {
        pathInCalib = argv[readIndex++ + modifier];
        pathOutGoodMatches = argv[readIndex++ + modifier];
        pathOutInfo = argv[readIndex++ + modifier];
        if (argc >= minArgcE + 2) {
            pathOutInlierMatchesImage = argv[readIndex++ + modifier];
            pathOutOutlierMatchesImage = argv[readIndex++ + modifier];
            if (argc == minArgcE + 3) {
                pathOutEpiImage = argv[readIndex++ + modifier];
            }
        }
    }

    // Init random seed
    srand(seed);

    std::ifstream test_f(pathOutGoodMatches);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    Image<RGBColor> image1, image2;
    if (!libs::ReadImage(pathInImage1, &image1))
        return 1;
    if (!libs::ReadImage(pathInImage2, &image2))
        return 1;

    Image<unsigned char> image1Gray, image2Gray;
    libs::convertImage(image1, &image1Gray);
    libs::convertImage(image2, &image2Gray);

    std::vector<Match> vec_matchings;
    if (generateSift) {
        utility::SIFT(image1Gray,image2Gray, vec_matchings, 0.6f, 0);
        // Remove duplicates (frequent with SIFT)
        utility::rm_duplicates(vec_matchings);
    } else {
        // Read correspondence file
        if (Match::loadMatch(pathInMatches, vec_matchings))
            std::cout << "Read " << vec_matchings.size() << " matches" << std::endl;
        else {
            std::cerr << "Failed reading matches from " << pathInMatches << std::endl;
            return 1;
        }
    }

    // Remove duplicates (frequent with SIFT)
    utility::rm_duplicates(vec_matchings);

    const int w1 = int(image1Gray.Width()), h1 = int(image1Gray.Height());
    const int w2 = int(image2Gray.Width()), h2 = int(image2Gray.Height());

    // Estimation of homography with LRT
    bool ok = false;

    orsa::ModelEstimator *model;

    std::vector<Match> vec_matchingsNormalised; // Matches in camera coordinates.
    libNumerics::matrix<double> intrinsics_source(3, 3), // The intrinsic parameters of the source camera
    intrinsics_destination(3, 3); // The intrinsic parameters of the destination camera

    // Creating the model to estimate:
    switch (modelToAnalyse) {
        case 0:
            if (verbose) {
                std::cout << "Running a Homography estimation..." << std::endl;
            }
            model = new orsa::HomographyModel(vec_matchings, w1, h1, w2, h2, true);
            break;
        case 1:
            if (verbose) {
                std::cout << "Running a Fundamental matrix estimation..." << std::endl;
            }
            model = new orsa::FundamentalModel(vec_matchings, w1, h1, w2, h2, true);
            break;
        case 2:
            if (verbose) {
                std::cout << "Running a Essential matrix estimation..." << std::endl;
            }
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

    // We use the AC-RANSAC algorithm to estimate the model.
    orsa::Orsa algo(model);
    algo.setHyperParameters(precision);

    if (verbose) {
        std::cout << "... with algorithm AC-Ransac." << std::endl;
    }

    orsa::RansacAlgorithm::RunResult res;
    double runtime;
    ok = algo.evalModel(res, runtime, iterMax, verbose);

    std::vector<double> errorsInliers, errorsOutliers;
    std::vector<int> indexOutliers;
    for (int i = 0; i < vec_matchings.size(); i++){
        if (std::find(res.vInliers.begin(), res.vInliers.end(), i) == res.vInliers.end()) {
            indexOutliers.push_back(i);
        }
    }
    utility::computeModelError(res.vInliers, model, res.model, errorsInliers);
    utility::computeModelError(indexOutliers, model, res.model, errorsOutliers);
    std::ofstream fin(pathOutInlierError);
    if (fin.is_open()) {
        std::vector<double>::const_iterator VectIt = errorsInliers.begin();
        for (; VectIt != errorsInliers.end(); ++VectIt) {
            fin << *VectIt << " ";
        }
    }
    fin.close();
    std::ofstream fout(pathOutOutlierError);
    if (fout.is_open()) {
        std::vector<double>::const_iterator VectIt = errorsOutliers.begin();
        for (; VectIt != errorsOutliers.end(); ++VectIt) {
            fout << *VectIt << " ";
        }
    }
    fout.close();

    double err = ok ? model->ErrorStats(res.vInliers, res.model).first : 0;

    int inlierCount = res.vInliers.size();

    if (ok && verbose) {
        std::cout << "Model Params=" << res.model << std::endl;
        std::cout << "Inlier error: " << err << "\nComputed Sigma: " << res.sigma
                  << "\nNumber of inliers: " << inlierCount
                  << "\nNumber of iterations: " << res.T << "\nVPM: " << res.vpm << "\nRuntime: " << runtime
                  << std::endl;
    }

    if (!ok && verbose) {
        std::cerr << "Failed to estimate a model" << std::endl;
    }

    std::vector<Match> good_match;
    std::vector<int>::const_iterator it = res.vInliers.begin();
    for (; it != res.vInliers.end(); it++)
        good_match.push_back(vec_matchings[*it]);

    // Output info :
    utility::Rect intersection;
    std::vector<Match> newMatchings;
    libNumerics::matrix<double> modelParamF = model->toPixelSpace(res.model);

    // Preparation of the random generator:
    std::default_random_engine generator;
    generator.seed(seed);
    if (ok) {
        switch (modelToAnalyse) {
            case 0:
                if (utility::IntersectionBox(w1, h1, w2, h2, res.model, intersection) &&
                    intersection.Width() > 0 && intersection.Height() > 0) {
                    Image<RGBColor> warpedImage = apply_homography(image1, image2, res.model, intersection,
                                                                   res.vInliers, vec_matchings, newMatchings, generator,
                                                                   0, 0);

                    if (!Match::saveMatch(pathOutGoodMatches, newMatchings)) {
                        std::cerr << "Failed saving inliers into " << pathOutGoodMatches << std::endl;
                    }

                    saveArtificialGenerationParams(pathOutInfo, seed, precision, iterMax, 0, 0, 0, 0, 0, 0, w1, h1,
                                                   warpedImage.Width(), warpedImage.Height(), res.model, err,
                                                   inlierCount, res.sigma);

                }
                break;
            case 1:
            case 2:
                apply_fundamental(image2, modelParamF, res.vInliers, vec_matchings, newMatchings);

                if (!Match::saveMatch(pathOutGoodMatches, newMatchings)) {
                    std::cerr << "Failed saving inliers into " << pathOutGoodMatches << std::endl;
                }

                saveArtificialGenerationParams(pathOutInfo, seed, precision, iterMax, 0, 0, 0, 0, 0, 0, w1, h1, w2,
                                               h2, res.model, err, inlierCount, res.sigma);

                if (argc > 7) // Output images TODO
                {
                    const char *fileIn = 0, *fileOut = 0, *fileEpi = 0;
                    fileIn = pathOutInlierMatchesImage;
                    fileOut = pathOutOutlierMatchesImage;
                    if (argc != 7)
                        fileEpi = pathOutEpiImage;
//                    utility::fundamental_graphical_output(image1Gray, image2Gray,
//                                                          vec_matchings, res.vInliers, &modelParamF,
//                                                          fileIn, fileOut, fileEpi);
                }
                break;
            default:
                std::cerr << "Model number not supported. Choose in {0, 1, 2}" << std::endl;
                delete model;
                return 1;
                break;
        }
    }
    delete model;

    return 0;
}
