//
// Created by clementriu on 8/31/20.
//

#include "read_write_experiment.hpp"

#include "utilities/data_handler.hpp"

#include <fstream>
#include <vector>


//// Save all parameters of an artificial match generation to a file.
//// \param[in] fileName: Emplacement of the file.
//// \param[in] seed: Value of the random seed.
//// \param[in] iterMax: Maximum number of iteration allowed during computation.
//// \param[in] noiseType: When inlier noise is applied : 0 for no noise, 1 after model, -1 before model.
//// \param[in] stdNoise: Std of the inlier noise.
//// \param[in] outlierType: Outlier type : 0 for no outliers, 1 for uniform outliers, 2 for uniform outliers that can't be inliers.
//// \param[in] outlierRatio: Ratio of outliers.
//// \param[in] outlierThreshold: Inlier/Outlier threshold.
//// \param[in] maxOutlier: Maximum number of iterations when adding outliers.
//// \param[in] w1: h1: w2: h2: Dimensions of images 1 and 2.
//// \param[in] modelParams: Used model matrix.
//// \param[in] modelError: Error of the used model.
//// \param[in] inlierCount: Number of inliers of the model.
//// \param[in] computedSigma: Estimated inlier/outlier threshold.
//// \return True if file was successfully opened, False, otherwise.
bool saveArtificialGenerationParams(const char *fileName, const unsigned int seed,
                                    const double precision, const int iterMax,
                                    const int noiseType, const double stdNoise,
                                    const int outlierType, const double outlierRatio, const double outlierThreshold,
                                    const int maxOutlier,
                                    const int w1, const int h1, const int w2, const int h2,
                                    const libNumerics::matrix<double> &modelParams, const double modelError,
                                    const int inlierCount, const double computedSigma) {
    std::ofstream f(fileName);

    if (f.is_open()) {
        f << "Seed " << seed << "\n";
        f << "Max precision " << precision << "\n";
        f << "Max iteration " << iterMax << "\n";
        f << "Noise type " << noiseType << "\n";
        f << "Noise std " << stdNoise << "\n";
        f << "Outlier type " << outlierType << "\n";
        f << "Outlier ratio " << outlierRatio << "\n";
        f << "Outlier threshold " << outlierThreshold << "\n";
        f << "Max outlier rejected " << maxOutlier << "\n";
        f << "w1 " << w1 << "\n";
        f << "h1 " << h1 << "\n";
        f << "w2 " << w2 << "\n";
        f << "h2 " << h2 << "\n";
        f << "Estimated model params " << modelParams << "\n";
        f << "Inlier model error " << modelError << "\n";
        f << "Inlier count " << inlierCount << "\n";
        f << "Computed Sigma " << computedSigma << "\n";
    }
    return f.is_open();
}

//// Save all parameters of an artificial match generation to a file.
//// \param[in] fileName: Emplacement of the file.
//// \param[out] w1: h1: w2: h2: Dimensions of images 1 and 2.
//// \param[out] modelParams: Used model matrix.
//// \param[out] inlierCount: Number of inliers of the model.
//// \return True if file was successfully opened, False, otherwise.
bool readArtificialGenerationParams(const char *fileName,
                                    int &w1, int &h1, int &w2, int &h2,
                                    libNumerics::matrix<double> &modelParams, int &inlierCount) {
    std::fstream file;
    std::string word;
    file.open(fileName);

    std::vector<std::string> readParamsStr;

    while (file >> word) { //take word and print
        if (word == "w1") {
            int value;
            file >> value;
            w1 = value;
        }
        if (word == "h1") {
            int value;
            file >> value;
            h1 = value;
        }
        if (word == "w2") {
            int value;
            file >> value;
            w2 = value;
        }
        if (word == "h2") {
            int value;
            file >> value;
            h2 = value;
        }
        if (word == "[") {
            file >> word;
            while (word != "]") {
                readParamsStr.push_back(word);
                file >> word;
            }
        }
        if (word == "Inlier") {
            file >> word;
            if (word == "count") {
                int value;
                file >> value;
                inlierCount = value;
            }
        }
    }

    std::vector<double> readParamsDouble;
    std::vector<std::string>::const_iterator itReadParamsStr = readParamsStr.begin();
    for (; itReadParamsStr != readParamsStr.end(); itReadParamsStr++) {
        std::string strValue = *itReadParamsStr;
        if (strValue.back() == ';') {
            strValue = strValue.substr(0, strValue.size() - 1);
        }
        readParamsDouble.push_back(atof(strValue.c_str()));
    }
    int nRow, nCol;
    nRow = modelParams.nrow();
    nCol = modelParams.ncol();
    for (int ind = 0; ind < readParamsDouble.size(); ind++) {
        int iRow = ind / nCol;
        int iCol = ind % nCol;
        modelParams(iRow, iCol) = readParamsDouble[ind];
    }
    return true;
}

//// Read a set of points and mixes them.
//// \param[in] fileInliers: Path to the inlier matches.
//// \param[in] fileOutliers: Path to the oulier matches.
//// \param[in] nGen: Number of dataset to read.
//// \param[out] pointsAll: All points read.
//// \param[out] groundTruthLabelsAll: All labels: 0 if from outlier, 1 if from inlier.
//// \param[out] inliersAll: Only inliers.
//// \param[out] outliersAll: Only outliers.
//// \param[in] readOutliers: If False, no outliers will be read.
bool ReadPoints(const char *fileInliers, const char *fileOutliers, const int nGen,
                std::vector<std::vector<Match>> &pointsAll, std::vector<std::vector<int>> &groundTruthLabelsAll,
                std::vector<std::vector<Match>> &inliersAll, std::vector<std::vector<Match>> &outliersAll,
                const bool readOutliers) {
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
    const int numPoints = ptsInliersAll.size();

    assert(numPoints % nGen == 0);
    const int numPointsPerGen = numPoints / nGen;

    for (int gen = 0; gen < nGen; gen++) {
        std::vector<Match> ptsInliers;
        std::vector<Match> ptsOutliers;

        for (int i = gen * numPointsPerGen; i < (gen + 1) * numPointsPerGen; i++) {
            ptsInliers.push_back(ptsInliersAll[i]);
            if (readOutliers) {
                ptsOutliers.push_back(ptsOutliersAll[i]);
            }
        }

        std::vector<Match> ptsMixed;
        std::vector<int> groundTruthLabels;

        utility::randomDataFusionWithMemory(ptsInliers, ptsOutliers, ptsMixed, 1, 0, groundTruthLabels);

        pointsAll.push_back(ptsMixed);
        inliersAll.push_back(ptsInliers);
        outliersAll.push_back(ptsOutliers);
        groundTruthLabelsAll.push_back(groundTruthLabels);
    }
    return true;
}