//
// Created by riuclement on 9/9/21.
//

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "libImage/image_io.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/pnp_model.hpp"

#include "libOrsa/orsa.hpp"

#include "utilities/cmdLine.h"
#include "utilities/Rect.hpp"
#include "utilities/metrics.hpp"
#include "utilities/siftMatch.hpp"
#include "utilities/warping.hpp"
#include "utilities/usac_file_handler.hpp"

#include "apply_transform.hpp"
#include "read_write_experiment.hpp"


int main(int argc, char **argv) {
    int iterMax = 10000; // Maximum number of iterations of the algorithm.
    double precision = 0; // Default pixel precision of the algorithm.

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    int modelToAnalyse = 0; // 0 for PNP.

    utility::CmdLine cmd;
    cmd.add(utility::make_option('m', modelToAnalyse, "model-analyse")
                    .doc("Model to analyse: 0 for PNP."));

    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of AC-RANSAC."));

    cmd.add(utility::make_option('p', precision, "precision")
                    .doc("Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of RANSAC, AC-RANSAC and LRT."));

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

    bool forceCompute = cmd.used('f');
    bool verbose = cmd.used('v');

    int minArgc = 8;

    if (modelToAnalyse != 0 || argc != minArgc) { // TODO
        std::cerr << "Usage: " << argv[0] << " [options] inMatches imgIn "
                  << "calibrationMatrix "
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

    const char *pathInMatches = argv[1]; // 2D3D matches file.
    const char *pathInImage = argv[2]; // Image.
    const char *pathInCalib = argv[3]; // Calibration matrix.
    const char *pathOutGoodMatches = argv[4]; // Output file with good matches coordinate.
    const char *pathOutInfo = argv[5]; // Output file with info about the run.
    const char *pathOutInlierError = argv[6];
    const char *pathOutOutlierError = argv[7];

    // Init random seed
    srand(seed);

    std::ifstream test_f(pathOutGoodMatches);
    if (!forceCompute && test_f.good()) {
        std::cout << "Already computed!" << std::endl;
        return 0;
    }

    Image<RGBColor> image;
    if (!libs::ReadImage(pathInImage, &image))
        return 1;

    Image<unsigned char> imageGray;
    libs::convertImage(image, &imageGray);

    std::vector<Match2D3D> vec_matchings;
    // Read correspondence file
    if (Match2D3D::loadMatch2D3D(pathInMatches, vec_matchings))
        std::cout << "Read " << vec_matchings.size() << " matches" << std::endl;
    else {
        std::cerr << "Failed reading matches from " << pathInMatches << std::endl;
        return 1;
    }

    // Remove duplicates (frequent with SIFT)
    utility::rm_duplicates(vec_matchings);

    const int w = int(imageGray.Width()), h = int(imageGray.Height());

    // Estimation of homography with LRT
    bool ok = false;

    orsa::ModelEstimator *model;

    std::vector<Match2D3D> vec_matchingsNormalised; // Matches in camera coordinates.
    libNumerics::matrix<double> intrinsics(3, 3); // The intrinsic parameters of the source camera

    // Creating the model to estimate:
    switch (modelToAnalyse) {
        case 0:
            if (verbose) {
                std::cout << "Running PNP estimation..." << std::endl;
            }
            if (utility::loadCalibration(pathInCalib, intrinsics))
                std::cout << "Read " << pathInCalib << " calibration" << std::endl;
            else {
                std::cerr << "Failed reading calibration matrix from " << pathInCalib << std::endl;
                return 1;
            }

            orsa::normalisePoints(vec_matchings, intrinsics, vec_matchingsNormalised);

            model = new orsa::PnPModel(vec_matchings, vec_matchingsNormalised, intrinsics, w, h);
            break;
        default:
            std::cerr << "Model number not supported. Choose in {0}" << std::endl;
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

    double err = ok ? model->ErrorStats(res.vInliers, res.model).first : 0;

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

    std::vector<Match2D3D> good_match;
    std::vector<int>::const_iterator it = res.vInliers.begin();
    for (; it != res.vInliers.end(); it++)
        good_match.push_back(vec_matchings[*it]);

    // Output info :
    utility::Rect intersection;
    std::vector<Match2D3D> newMatchings;
    libNumerics::matrix<double> modelParamF = model->toPixelSpace(res.model);
    if (ok) {
        switch (modelToAnalyse) {
            case 0:
                applyPnP(image, res.model, res.vInliers, vec_matchingsNormalised, intrinsics, newMatchings);

                if (!Match2D3D::saveMatch2D3D(pathOutGoodMatches, newMatchings)) {
                    std::cerr << "Failed saving inliers into " << pathOutGoodMatches << std::endl;
                }

                saveArtificialGenerationParams(pathOutInfo, seed, precision, iterMax, 0, 0, 0, 0, 0, 0, w, h,
                                               0, 0, res.model, err, inlierCount, res.sigma);
                break;
            default:
                std::cerr << "Model number not supported. Choose in {0}" << std::endl;
                delete model;
                return 1;
                break;
        }
    }
    delete model;

    return 0;
}
