/**
* @file demo.cpp
* @brief Launch 4 Ransac variants and 3 models
* @author Clement Riu, Pascal Monasse
*
* Copyright (c) 2021 Clement Riu
* Copyright (c) 2021 Pascal Monasse
* All rights reserved.
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "libImage/image_io.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/homography_model.hpp"
#include "libOrsa/fundamental_model.hpp"
#include "libOrsa/essential_model.hpp"

#include "libOrsa/ransac.hpp"
#include "libOrsa/orsa.hpp"
#include "libOrsa/orsa_fast.hpp"
#include "libOrsa/lrtsac.hpp"
#include "libOrsa/magsac.hpp"
#include "libOrsa/muse.hpp"
#include "libOrsa/starsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/siftMatch.hpp"
#include "utilities/warping.hpp"
#include "utilities/homography_graphical_output.hpp"
#include "utilities/fundamental_graphical_output.hpp"
#include "utilities/usac_file_handler.hpp"

/// Put full name, since the input name could be shortened.
bool expand_prefix(std::string &name, const std::vector<std::string> &V) {
    std::string out;
    std::vector<std::string>::const_iterator it = V.begin();
    for (; it != V.end(); ++it)
        if (name.size() <= it->size() && name == it->substr(0, name.size())) {
            if (out.empty()) out = *it;
            else {
                std::cerr << "Prefix " << name << " is ambiguous: could be "
                          << out << " or " << *it << std::endl;
                return false;
            }
        }
    if (out.empty()) {
        std::cerr << "Unknown model/algo " << name << ". Choices are: ";
        for (it = V.begin(); it != V.end(); ++it) std::cerr << *it << ' ';
        std::cerr << std::endl;
        return false;
    }
    name = out;
    return true;
}

// Check region is inside image frame
bool check_region(const utility::Geometry &region, const Image<RGBColor> &im) {
    utility::Geometry zone;
    zone.x0 = zone.y0 = 0;
    zone.x1 = int(im.Width());
    zone.y1 = int(im.Height());
    if (!(region.x0 < region.x1 && region.y0 < region.y1 &&
          zone.inside(region.x0, region.y0) &&
          zone.inside(region.x1 - 1, region.y1 - 1)))
        return false;
    return true;
}

// Demo file that can run any algorithm on any model.
int main(int argc, char **argv) {
    float fSiftRatio = 0.6f; // SIFT parameter. Only used if -r is not used.
    utility::Geometry region;
    region.x0 = region.y0 = region.x1 = region.y1 = 0;

    int iterMax = 50000; // Maximum number of iterations of the algorithm
    double precision = 3; // Default pixel precision of the algorithm
    double cpIIT = 0.99; // Default confidence wrt type II error (adjust #iters)

    int nModelMinRansac = 1;

    double cpI = 0.99;
    double cpIIB = 0.95;

    size_t partitionNumber = 5; // Magsac partition number.

    std::string sInK; // Essential only: K matrix file of both images

    bool readMatches = false, verbose = false;
    unsigned int seed = (unsigned) time(0); // Use -t for a reproducible run

    std::string modelClass="Homography"; // Homography/Fundamental/Essential
    std::string algo="Ransac"; // Ransac/AC-Ransac/fast-AC-Ransac/LRT/StarSac/Muse/MAGSAC/PlusPlusMAGSAC

    utility::CmdLine cmd;
    cmd.add(make_option('c', region, "cut")
                    .doc("cut region of imagInA: wxh+x+y = rect [x,x+w] x [y,y+h]"));
    cmd.add(utility::make_option('r', readMatches, "read")
                    .doc("Read file of matches allMatches.txt, do not use SIFT"));
    cmd.add(utility::make_option('s', fSiftRatio, "sift")
                    .doc("SIFT distance ratio of descriptors"));

    cmd.add(utility::make_option('m', modelClass, "model")
                    .doc("Model class: Homography/Fundamental/Essential"));
    cmd.add(utility::make_option('a', algo, "algo")
            .doc("Algorithm to use: Ransac/AC-Ransac/fast-AC-Ransac/LRT/StarSac/Muse/MAGSAC/PlusPlusMAGSAC"));

    cmd.add(utility::make_option('i', iterMax, "iterMax")
                    .doc("Number of iterations of the algorithm."));

    cmd.add(utility::make_option('p', precision, "precision")
                    .doc("Max precision (in pixels) of registration (0=arbitrary)"));
    cmd.add(utility::make_option(0, cpIIT, "cpIIT")
                    .doc("Confidence against type II error (Ransac/LRT)"));

    cmd.add(utility::make_option('n', nModelMinRansac, "nModelMinRansac")
                    .doc("Min number of models before terminating ransac"));

    cmd.add(utility::make_option(0, cpI, "cpI")
                    .doc("Confidence proba wrt type I error (LRT)"));
    cmd.add(utility::make_option(0, cpIIB, "cpIIB")
                    .doc("Confidence proba wrt bailout (LRT)"));

    cmd.add(utility::make_option(0, partitionNumber, "partition-number")
                    .doc("Number of partition of the Sigma consensus algorithm."));

    cmd.add(utility::make_option('k', sInK, "Kfile")
                    .doc("File for calibration matrices (Essential only)"));

    cmd.add(utility::make_option('v', verbose, "verbose")
                    .doc("Print info during the run."));
    cmd.add(utility::make_option('t', seed, "time-seed")
                    .doc("Use value instead of time for random seed."));

    try {
        cmd.process(argc, argv);
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return 1;
    }

    std::string M[] = {"Homography", "Fundamental", "Essential"};
    size_t n = sizeof(M) / sizeof(M[0]);
    if (!expand_prefix(modelClass, std::vector<std::string>(M, M + n)))
        return 1;
    std::string A[] = {"Ransac", "AC-Ransac", "fast-AC-Ransac", "LRT", "StarSac", "Muse", "MAGSAC", "PlusPlusMAGSAC"};
    n = sizeof(A)/sizeof(A[0]);
    if(! expand_prefix(algo, std::vector<std::string>(A,A+n)))
        return 1;
    if (argc!=5 && argc!=7 && argc!=8 && (modelClass!="Homography"||argc!=10)) {
        std::cerr << "Usage: " << argv[0] << " [options] imgInA imgInB "
                  << "allInOutMatches.txt inlierOutMatches.txt [optImgOut]\n"
                  << "- imgInA, imgInB: the two input images (JPG/PNG format)\n"
                  << "- allInOutMatches.txt: output (input if -r) text file of "
                     "format \"x1 y1 x2 y2\"\n"
                  << "- inlierOutMatches.txt: output, but only with inliers.\n"

                  << "\t[optImgOut] (output images): "
                  << "inliers outliers [mosaic/epi [regA regB]]\n"
                  << "- inliers, outliers: lines for inliers, "
                  << "lines for outliers and their error\n"
                  << "- mosaic/epi: mosaic if Homography else epipolar lines\n"
                  << "- regA, regB: registered images (Homography only)\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }
    if (cmd.used('k')) {
        if (modelClass != "Essential") {
            std::cerr << "Option -k only used for Essential" << std::endl;
            return 1;
        }
    } else {
        if (modelClass == "Essential") {
            std::cerr << "Option -k required for Essential" << std::endl;
            return 1;
        }
    }
    bool bUseRegion = cmd.used('c');

    // Required program arguments
    const char *sInImage1 = argv[1]; // First image path
    const char *sInImage2 = argv[2]; // Second image path
    const char *sInOutMatches = argv[3]; // Path to all matches, input or output
    const char *sOutGoodMatches = argv[4]; // Path to only good matches

    // Optional output images:
    const char *sOutImgInliers = 0;  // All models: image with inliers
    const char *sOutImgOutliers = 0; // All models: image with outliers
    const char *sOutImgEpi = 0;      // F/E only: image with epipolar lines
    const char *sOutImgMosaic = 0;   // H only: image with reconstructed mosaic
    const char *sOutImgReg1 = 0;     // H only: warped first image
    const char *sOutImgReg2 = 0;     // H only: warped second image
    if (argc >= 7) {
        sOutImgInliers = argv[5];
        sOutImgOutliers = argv[6];
    }
    if (argc >= 8) sOutImgMosaic = sOutImgEpi = argv[7];
    if (argc >= 10) {
        sOutImgReg1 = argv[8];
        sOutImgReg2 = argv[9];
    }

    // Init random seed
    srand(seed);
    if (verbose && !cmd.used('t'))
        std::cout << "Rerun with \"-t " << seed << "\" to reproduce" << std::endl;

    Image<RGBColor> image1, image2;
    if (!libs::ReadImage(sInImage1, &image1))
        return 1;
    if (!libs::ReadImage(sInImage2, &image2))
        return 1;

    if (bUseRegion && !check_region(region, image1)) { // Sanity check
        std::cout << "Invalid cut region " << region << std::endl;
        return 1;
    }

    Image<unsigned char> image1Gray, image2Gray;
    libs::convertImage(image1, &image1Gray);
    libs::convertImage(image2, &image2Gray);

    // Find matches with SIFT or read correspondence file
    std::vector<Match> vMatch;
    if (readMatches) {
        if (Match::loadMatch(sInOutMatches, vMatch))
            std::cout << "Read " << vMatch.size() << " matches" << std::endl;
        else {
            std::cerr << "Failed reading from " << sInOutMatches << std::endl;
            return 1;
        }
    } else
        SIFT(image1Gray, image2Gray, vMatch, fSiftRatio, bUseRegion ? &region : 0);

    if (bUseRegion) {
        Image<RGBColor> image = image1;
        Crop(image, region.x0, region.y0,
             region.x1 - region.x0, region.y1 - region.y0, image1);
        libs::convertImage(image1, &image1Gray);
    }

    // Remove duplicates (frequent with SIFT)
    utility::rm_duplicates(vMatch);

    // Save match files
    if (!readMatches && !Match::saveMatch(sInOutMatches, vMatch)) {
        std::cerr << "Failed saving matches into " << sInOutMatches << std::endl;
        return 1;
    }

    const int w1 = int(image1Gray.Width()), h1 = int(image1Gray.Height());
    const int w2 = int(image2Gray.Width()), h2 = int(image2Gray.Height());

    // Creation of the generic model and generic algorithm:
    orsa::ModelEstimator *model = 0;
    orsa::RansacAlgorithm *algorithm = 0;

    if (modelClass == "Homography")
        model = new orsa::HomographyModel(vMatch, w1, h1, w2, h2, true);
    else if (modelClass == "Fundamental")
        model = new orsa::FundamentalModel(vMatch, w1, h1, w2, h2, true);
    else if(modelClass == "Essential") {
        libNumerics::matrix<double> K1(3,3), K2(3,3); // Calibration matrices of cameras
        if (! utility::loadCalibration(sInK.c_str(), K1, K2)) {
            std::cerr << "Failed reading K matrices from " << sInK <<std::endl;
            return 1;
        }
        std::cout << "Read calibration file " << sInK << std::endl;
        std::vector<Match> vMatchNormalised;
        orsa::normalisePoints(vMatch, K1, K2, vMatchNormalised);
        model = new orsa::EssentialModel(vMatch, vMatchNormalised,
                                         w1, h1, w2, h2, K1, K2, true);
    }
    std::cout << "Model: " << modelClass;

    if (algo == "Ransac") {
        orsa::Ransac *a = new orsa::Ransac(model);
        a->setHyperParameters(precision, cpIIT, nModelMinRansac);
        algorithm = a;
    } else if (algo == "AC-Ransac") {
        orsa::Orsa *a = new orsa::Orsa(model);
        a->setHyperParameters(precision);
        algorithm = a;
    } else if(algo == "fast-AC-Ransac") {
        orsa::OrsaFast* a = new orsa::OrsaFast(model);
        a->setHyperParameters(precision);
        algorithm = a;
    } else if(algo == "LRT") {
        orsa::LRTSac* a = new orsa::LRTSac(model);
        a->setHyperParameters(cpI,cpIIB,cpIIT,precision,true);
        algorithm = a;
    } else if(algo == "StarSac") {
        orsa::StarSac* a = new orsa::StarSac(model);
        a->setHyperParameters(cpIIT);
        algorithm = a;
    } else if(algo == "Muse") {
        orsa::Muse* a = new orsa::Muse(model);
        a->setHyperParameters(cpIIT);
        algorithm = a;
    } else if(algo == "MAGSAC") {
        orsa::MAGSAC *a = new orsa::MAGSAC(model);
        a->setHyperParameters(cpIIT, precision, partitionNumber);
        algorithm = a;
    } else if(algo == "PlusPlusMAGSAC") {
        orsa::MAGSAC *a = new orsa::MAGSAC(model, orsa::MAGSAC::Version::MAGSAC_PLUS_PLUS);
        a->setHyperParameters(cpIIT, precision, partitionNumber);
        algorithm = a;
    }
    std::cout << ", Algorithm: " << algo << std::endl;

    // Estimation of the model parameters:
    orsa::RansacAlgorithm::RunResult res;
    double runtime;
    bool ok = algorithm->evalModel(res, runtime, iterMax, verbose);

    if (ok) {
        std::cout << "Result=" << res.model << std::endl;
        std::cout << "Inliers: " << res.vInliers.size() << "/" << vMatch.size() << "="
                  << (100 * res.vInliers.size()) / vMatch.size() << "%" << std::endl;
    }
    if (ok && verbose) {
        std::cout << "Sigma: " << res.sigma << "\nIterations: " << res.T
                  << "\nVerif/model: " << res.vpm
                  << "\nRuntime (s): " << runtime
                  << std::endl;
    }

    // Save output file of inliers
    std::vector<Match> inliers;
    std::vector<int>::const_iterator it = res.vInliers.begin();
    for (; it != res.vInliers.end(); it++)
        inliers.push_back(vMatch[*it]);
    if (!Match::saveMatch(sOutGoodMatches, inliers)) {
        std::cerr << "Failed saving matches into " << sOutGoodMatches << std::endl;
        delete model;
        delete algorithm;
        return 1;
    }

    if (!ok)
        std::cerr << "Failed to estimate a model." << std::endl;

    // Optional output images:
    if (modelClass == "Homography") {
        if (sOutImgInliers && sOutImgOutliers)
            utility::homography_matches_output(image1Gray, image2Gray,
                                               vMatch, res.vInliers,
                                               &res.model,
                                               sOutImgInliers, sOutImgOutliers);
        if (sOutImgMosaic) {
            utility::Rect intersect;
            if (utility::IntersectionBox(w1, h1, w2, h2, res.model, intersect) &&
                intersect.Width() > 0 && intersect.Height() > 0) {
                utility::homography_registration_output(image1, image2,
                                                        res.model, intersect,
                                                        sOutImgMosaic,
                                                        sOutImgReg1,
                                                        sOutImgReg2, true);
            }
        }
    }
    if ((modelClass == "Fundamental" || modelClass == "Essential") &&
        sOutImgInliers && sOutImgOutliers) { // Output images
        libNumerics::matrix<double> modelF = model->toPixelSpace(res.model);
        utility::fundamental_graphical_output(image1Gray, image2Gray,
                                              vMatch, res.vInliers, &modelF,
                                              sOutImgInliers, sOutImgOutliers,
                                              sOutImgEpi);
    }

    delete model;
    delete algorithm;
    return 0;
}
