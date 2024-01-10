//
// Created by clementriu on 10/26/20.
//

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "libImage/image_io.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/homography_model.hpp"
#include "libOrsa/libNumerics/homography.h"

#include "utilities/cmdLine.h"
#include "utilities/data_handler.hpp"
#include "utilities/siftMatch.hpp"
#include "utilities/warping.hpp"
#include "utilities/homography_graphical_output.hpp"

template<typename T>
void printVector(const std::vector<T> vectToPrint) {
    typename std::vector<T>::const_iterator itVectToPrint = vectToPrint.begin();
    for (; itVectToPrint != vectToPrint.end(); ++itVectToPrint) {
        std::cout << *itVectToPrint << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char **argv) {
    int errorThreshold = 0;

    unsigned int seed = (unsigned) time(0); // Use option -t for a reproducible run

    utility::CmdLine cmd;
    utility::Geometry region;
    region.x0 = region.y0 = region.x1 = region.y1 = 0;
    cmd.add(utility::make_option('c', region, "cut")
                    .doc("cut region of imagInA: wxh+x+y = rect [x,x+w] x [y,y+h]"));

    cmd.add(utility::make_switch('o', "ignoreOutlier")
                    .doc("Does not read outliers if set."));
    cmd.add(utility::make_option('e', errorThreshold, "errorThreshold")
                    .doc("Error threshold for the inliers matches print. If set > 0 other 3 matches will be outputed."));

    cmd.add(utility::make_option('t', seed, "time-seed")
                    .doc("Use value instead of time for random seed (for debug)"));

    try {
        cmd.process(argc, argv);
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return 1;
    }

    const int maxArg = 11;
    if (argc != maxArg) { // TODO
        std::cerr << "Usage: " << argv[0] << " [options] homography.txt imgInA.jpg imgInB.jpg "
                  << "imgMosaic imgMosaicA imgMosaicB\n"
                  << "- imgInA, imgInB: the two input images (JPG or PNG format)\n"
                  << "- allMatches.txt: output (input if option -r) text file of "
                     "format \"x1 y1 x2 y2\"\n"
                  << "- orsaMatches.txt: output, but only with inliers.\n"
                  << "- imgInliers (optional): output image showing inliers\n"
                  << "- imgOutliers (optional): output image showing outliers and "
                     "their error\n"
                  << "- imgMosaic (optional): output mosaic image\n"
                  << "- imgMosaicA,imgMosaicB (optional): registered images\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }
    bool readOutliers = !cmd.used('o');

    // Init random seed
    srand(seed);

    const char *inPathToHomography = argv[1];

    const char *inPathToImg1 = argv[2];
    const char *inPathToImg2 = argv[3];

    const char *inPathToInliers = argv[4];
    const char *inPathToOutliers = argv[5];

    const char *outPathToInliers = argv[6];
    const char *outPathToOutliers = argv[7];

    const char *outPathToMosaic = argv[8];
    const char *outPathToMosaic1 = argv[9];
    const char *outPathToMosaic2 = argv[10];

    Image<RGBColor> image1, image2;
    if (!libs::ReadImage(inPathToImg1, &image1))
        return 1;
    if (!libs::ReadImage(inPathToImg2, &image2))
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

    if (bUseRegion) {
        Image<RGBColor> image = image1;
        Crop(image, region.x0, region.y0,
             region.x1 - region.x0, region.y1 - region.y0, image1);
        libs::convertImage(image1, &image1Gray);
    }

    const int w1 = int(image1Gray.Width()), h1 = int(image1Gray.Height());
    const int w2 = int(image2Gray.Width()), h2 = int(image2Gray.Height());

    // Estimation of homography with ORSA
    libNumerics::matrix<double> H(3, 3);

    utility::readHomography(inPathToHomography, H);
    std::cout << "H=" << H << std::endl;

    H /= H(2, 2);
    std::cout << "H=" << H << std::endl;

    std::vector<Match> ptsInliers;
    std::vector<Match> ptsOutliers;

    if (!Match::loadMatch(inPathToInliers, ptsInliers)) {
        std::cerr << "Problem loading inliers: " << inPathToInliers << std::endl;
        return false;
    }
    if (readOutliers) {
        if (!Match::loadMatch(inPathToOutliers, ptsOutliers)) {
            std::cerr << "Problem loading outliers: " << inPathToOutliers << std::endl;
            return false;
        }
    }
    std::vector<Match> ptsMixed;
    std::vector<int> labels;

    utility::randomDataFusionWithMemory(ptsInliers, ptsOutliers, ptsMixed, 1, 0, labels);

    orsa::ModelEstimator *tempModel = new orsa::HomographyModel(ptsInliers, w1, h1, w2, h2);

    std::vector<double> errors;
    double maxError = 0;
    for (size_t i = 0; i < ptsInliers.size(); ++i) {
        errors.push_back(std::sqrt(tempModel->Error(H, i)));
        if (errors.back() > maxError) {
            maxError = errors.back();
        }
    }
    printVector(errors);
    std::cout << maxError << std::endl;

    std::vector<int> vec_inliers;
    for (size_t i = 0; i < labels.size(); i++) {
        if (labels[i]) {
            vec_inliers.push_back(i);
        }
    }

    // Sift de-duplicated output display
    utility::homography_matches_output(image1Gray, image2Gray,
                                       ptsMixed, vec_inliers,
                                       &H,
                                       outPathToInliers, outPathToOutliers);

    // Mosaics
    utility::Rect intersection;
    if (utility::IntersectionBox(w1, h1, w2, h2, H, intersection) &&
        intersection.Width() > 0 && intersection.Height() > 0) {
        utility::homography_registration_output(image1, image2, H, intersection,
                                                outPathToMosaic, outPathToMosaic1, outPathToMosaic2, true);
    }
    delete tempModel;
    return 0;
}
