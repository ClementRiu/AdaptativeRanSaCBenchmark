/**
* @file demo_PNP.cpp
* @brief Launch 4 Ransac variants and 1 model
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

#include "libOrsa/pnp_model.hpp"

#include "libOrsa/ransac.hpp"
#include "libOrsa/orsa.hpp"
#include "libOrsa/orsa_fast.hpp"
#include "libOrsa/lrtsac.hpp"
#include "libOrsa/starsac.hpp"
#include "libOrsa/muse.hpp"
#include "libOrsa/magsac.hpp"

#include "utilities/cmdLine.h"
#include "utilities/siftMatch.hpp"
#include "utilities/warping.hpp"
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
//    utility::Geometry region;

    int iterMax = 50000; // Maximum number of iterations of the algorithm
    double precision = 3; // Default pixel precision of the algorithm
    double cpIIT = 0.99; // Default confidence wrt type II error (adjust #iters)

    int nModelMinRansac = 1;

    double cpI = 0.99;
    double cpIIB = 0.95;

    size_t partitionNumber = 5; // Magsac partition number.

    bool verbose = false;
    unsigned int seed = (unsigned) time(0); // Use -t for a reproducible run

    std::string modelClass = "PNP"; // PNP
    std::string algo = "Ransac"; // Ransac/AC-Ransac/fast-AC-Ransac/LRT/StarSac/Muse/MAGSAC/PlusPlusMAGSAC

    utility::CmdLine cmd;

    cmd.add(utility::make_option('m', modelClass, "model")
                    .doc("Model class: PNP"));
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

    std::string M[] = {"PNP"};
    size_t n = sizeof(M) / sizeof(M[0]);
    if (!expand_prefix(modelClass, std::vector<std::string>(M, M + n)))
        return 1;
    std::string A[] = {"Ransac", "AC-Ransac", "fast-AC-Ransac", "LRT", "StarSac", "Muse", "MAGSAC", "PlusPlusMAGSAC"};
    n = sizeof(A) / sizeof(A[0]);
    if (!expand_prefix(algo, std::vector<std::string>(A, A + n)))
        return 1;

    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " [options] images.jpg Calibration.txt "
                  << "allInMatches.txt inlierOutMatches.txt\n"
                  << "- allInMatches.txt: input text file of "
                     "format \"x1 y1 x2 y2 z2\"\n"
                  << "- inlierOutMatches.txt: output, but only with inliers.\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }


    // Required program arguments
    const char *sInImage = argv[1]; // Image path
    const char *kIn = argv[2]; // Calibration matrix path
    const char *sInMatches = argv[3]; // Path to all matches, input
    const char *sOutGoodMatches = argv[4]; // Path to only good matches

    // Init random seed
    srand(seed);
    if (verbose && !cmd.used('t'))
        std::cout << "Rerun with \"-t " << seed << "\" to reproduce" << std::endl;

    Image<RGBColor> image;
    if (!libs::ReadImage(sInImage, &image))
        return 1;

    Image<unsigned char> imageGray;
    libs::convertImage(image, &imageGray);

    // Find matches with SIFT or read correspondence file
    std::vector<Match2D3D> vMatch;
    if (Match2D3D::loadMatch2D3D(sInMatches, vMatch))
        std::cout << "Read " << vMatch.size() << " matches" << std::endl;
    else {
        std::cerr << "Failed reading from " << sInMatches << std::endl;
        return 1;
    }

    // Remove duplicates (frequent with SIFT)
    utility::rm_duplicates(vMatch);

    const int w = int(imageGray.Width()), h = int(imageGray.Height());

    // Creation of the generic model and generic algorithm:
    orsa::ModelEstimator *model = 0;
    orsa::RansacAlgorithm *algorithm = 0;

    if (modelClass == "PNP") {
        libNumerics::matrix<double> K(3, 3); // Calibration matrices of cameras
        if (!utility::loadCalibration(kIn, K)) {
            std::cerr << "Failed reading K matrix from " << kIn << std::endl;
            return 1;
        }
        std::cout << "Read calibration file " << kIn << std::endl;
        std::vector<Match2D3D> vMatchNormalised;
        orsa::normalisePoints(vMatch, K, vMatchNormalised);
        model = new orsa::PnPModel(vMatch, vMatchNormalised, K, w, h);
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
    } else if (algo == "fast-AC-Ransac") {
        orsa::OrsaFast *a = new orsa::OrsaFast(model);
        a->setHyperParameters(precision);
        algorithm = a;
    } else if (algo == "LRT") {
        orsa::LRTSac *a = new orsa::LRTSac(model);
        a->setHyperParameters(cpI, cpIIB, cpIIT, precision, true);
        algorithm = a;
    } else if (algo == "StarSac") {
        orsa::StarSac *a = new orsa::StarSac(model);
        a->setHyperParameters(cpIIT);
        algorithm = a;
    } else if (algo == "Muse") {
        orsa::Muse *a = new orsa::Muse(model);
        a->setHyperParameters(cpIIT);
        algorithm = a;
    } else if (algo == "MAGSAC") {
        orsa::MAGSAC *a = new orsa::MAGSAC(model);
        a->setHyperParameters(cpIIT, precision, partitionNumber);
        algorithm = a;
    } else if (algo == "PlusPlusMAGSAC") {
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
    std::vector<Match2D3D> inliers(res.vInliers.size());
    std::vector<int>::const_iterator it = res.vInliers.begin();
    for (int i = 0; it != res.vInliers.end(); it++, i++) {
        inliers[i] = vMatch[*it];
    }
    if (!Match2D3D::saveMatch2D3D(sOutGoodMatches, inliers)) {
        std::cerr << "Failed saving matches into " << sOutGoodMatches << std::endl;
        delete model;
        delete algorithm;
        return 1;
    }

    if (!ok) {
        std::cerr << "Failed to estimate a model" << std::endl;
        delete model;
        delete algorithm;
        return 0;
    }

    delete model;
    delete algorithm;
    return 0;
}
