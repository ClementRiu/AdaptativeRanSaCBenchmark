//
// Created by clementriu on 6/1/20.
//

#include "usac_file_handler.hpp"

#include <cstring>
#include <fstream>
#include <iostream>

namespace utility {

    //// Save experiments results to a file. Saves info, precisions, recalls, runtimes and estimated sigmas.
    //// Can take 4 different algorithms: Ransac, Ransac, AC-RANSAC, LRT.
    //// \param nameFile: Emplacement of the file.
    //// \param seed: The random seed used.
    //// \param nGen: Number of different dataset generated saved.
    //// \param nRun: Number of runs saved.
    //// \param pathToIn: Path to the input data of the experiment.
    //// \param algorithmNames: Name of the different algorithms for printing.
    //// \param <precision/recall/runtime/sigmasEstimated><Mean/STD/Vect>Vect: <Mean/Std/Vector> of values of the <precision/recall/runtime/estimated sigma> of the different algorithms in a vector.
    //// \param fullSave: if fullSave is true, the vectors are saved. Otherwise just general info, mean and std.
    bool saveExpInfo(const char *nameFile, const unsigned int seed, const int nGen, const int nRun,
                     const char *pathToIn, const std::vector<const char *> &algorithmNames,
                     const std::vector<double> &precisionMeanVect, const std::vector<double> &precisionStdVect,
                     const std::vector<std::vector<double>> &precisionVectVect,
                     const std::vector<double> &recallMeanVect, const std::vector<double> &recallStdVect,
                     const std::vector<std::vector<double>> &recallVectVect,
                     const std::vector<double> &runtimeMeanVect, const std::vector<double> &runtimeStdVect,
                     const std::vector<std::vector<double>> &runtimeVectVect,
                     const std::vector<std::vector<double>> &sigmaEstimatedVectVect,
                     const bool fullSave) {

        std::ofstream f(nameFile);

        if (f.is_open()) {
            f << "Input folder read: " << pathToIn << "\n";
            f << "Datasets generated with seed: " << seed << "\n";
            f << "\nResults over " << nGen << " random datasets and " << nRun << " runs.\n";

            for (size_t i = 0; i < algorithmNames.size(); i++) {
                f << algorithmNames[i] << "\n";
                f << " p= " << precisionMeanVect[i] << " | " << precisionStdVect[i]
                  << " r= " << recallMeanVect[i] << " | " << recallStdVect[i]
                  << " t= " << runtimeMeanVect[i] << " | " << runtimeStdVect[i] << "\n";
            }
            for (size_t i = 0; i < algorithmNames.size(); i++) {
                f << algorithmNames[i] << "\n";
                if (fullSave) {
                    f << "\n";
                    f << ":\n\tp:\n";
                    std::vector<double>::const_iterator itPrecision = precisionVectVect[i].begin();
                    for (; itPrecision != precisionVectVect[i].end(); ++itPrecision) {
                        f << *itPrecision << " ";
                    }
                    f << "\n\tr:\n";
                    std::vector<double>::const_iterator itRecall = recallVectVect[i].begin();
                    for (; itRecall != recallVectVect[i].end(); ++itRecall) {
                        f << *itRecall << " ";
                    }
                    f << "\n\tt:\n";
                    std::vector<double>::const_iterator itRuntime = runtimeVectVect[i].begin();
                    for (; itRuntime != runtimeVectVect[i].end(); ++itRuntime) {
                        f << *itRuntime << " ";
                    }
                    f << "\n\ts:\n";
                    std::vector<double>::const_iterator itSigma = sigmaEstimatedVectVect[i].begin();
                    for (; itSigma != sigmaEstimatedVectVect[i].end(); ++itSigma) {
                        f << *itSigma << " ";
                    }
                    f << "\n";
                }
            }
        }
        return f.is_open();
    }

    //// Save experiments results to a file. Saves info, number of iterations, vpm and runtimes.
    //// Takes up to 3 differents runs: LRT with all options, LRT without early bailout, LRT without any options.
    //// \param nameFile: Emplacement of the file.
    //// \param seed: The random seed used.
    //// \param nRun: Number of runs saved.
    //// \param pathToIn: Path to the input data of the experiment.
    //// \param algorithmNames: Name of the different algorithms for printing.
    //// \param <T/vpm><Mean/STD/Vect>Vect: <Mean/Std/Vector> of values of the <number of iterations/vpm> of the different algorithms in a vector.
    //// \param meanRuntime: mean runtime of each algorithm over the nRun runs.
    //// \param fullSave: if fullSave is true, the vectors are saved. Otherwise just general info, mean and std.
    bool saveExpInfo(const char *nameFile, const unsigned int seed, const int nRun,
                     const char *pathToIn, const std::vector<const char *> &algorithmNames,
                     const std::vector<double> &TMeanVect, const std::vector<double> &TStdVect,
                     const std::vector<std::vector<int>> &TVectVect,
                     const std::vector<double> &vpmMeanVect, const std::vector<double> &vpmStdVect,
                     const std::vector<std::vector<double>> &vpmVectVect,
                     const std::vector<double> &meanRuntime,
                     const bool fullSave) {

        std::ofstream f(nameFile);

        if (f.is_open()) {
            f << "Input folder read: " << pathToIn << "\n";
            f << "Datasets generated with seed: " << seed << "\n";
            f << "\nResults over " << nRun << " runs.\n";

            for (size_t i = 0; i < TMeanVect.size(); i++) {
                f << algorithmNames[i] << "\n";
                f << "\tT= " << TMeanVect[i] << " | " << TStdVect[i]
                  << "\tVPM= " << vpmMeanVect[i] << " | " << vpmStdVect[i]
                  << "\tt= " << meanRuntime[i] << "s\n";
            }
            f << "\n";
            for (size_t i = 0; i < TMeanVect.size(); i++) {
                f << algorithmNames[i] << "\n";
                if (fullSave) {
                    f << "\tT:\n";
                    std::vector<int>::const_iterator itT = TVectVect[i].begin();
                    for (; itT != TVectVect[i].end(); ++itT) {
                        f << *itT << " ";
                    }
                    f << "\n\tVPM:\n";
                    std::vector<double>::const_iterator itVpm = vpmVectVect[i].begin();
                    for (; itVpm != vpmVectVect[i].end(); ++itVpm) {
                        f << *itVpm << " ";
                    }
                    f << "\n";
                }
                f << "\n";
            }
        }
        return f.is_open();
    }

//    //// Save matches to a file.
//    //// \param pathToMatches: Emplacement of the file.
//    //// \param matchesVect: Vector of vectors of matches.
//    bool saveMatches(const char *pathToMatches, const std::vector<std::vector<Match>> &matchesVect) {
//        std::ofstream f(pathToMatches);
//        if (f.is_open()) {
//            std::vector<std::vector<Match>>::const_iterator matchesVectIt = matchesVect.begin();
//            for (; matchesVectIt != matchesVect.end(); ++matchesVectIt) {
//                std::vector<Match>::const_iterator matchIt = (*matchesVectIt).begin();
//                for (; matchIt != (*matchesVectIt).end(); ++matchIt) {
//                    f << *matchIt;
//                }
//                f << "\n";
//            }
//        }
//        return f.is_open();
//    }

    //// Read a calibration matrix, of the USAC format.
    //// \param[in] path: path to the calibration matrix.
    //// \param[out] K1, K2: both calibration matrix.
    //// \return True if 10 values where extracted, False otherwise.
    bool loadCalibration(const char *path, libNumerics::matrix<double> &K1, libNumerics::matrix<double> &K2) {
        std::ifstream infile(path);

        if (!infile.is_open())
            return false;

        double element;

        K1(1, 0) = 0;
        K1(2, 0) = 0;
        K1(2, 1) = 0;
        K1(2, 2) = 1;

        K2(1, 0) = 0;
        K2(2, 0) = 0;
        K2(2, 1) = 0;
        K2(2, 2) = 1;

        int count = 0;
        while (infile >> element) {
            if (count < 3) {
                K1(0, count) = element;
            }
            if (count >= 3 && count < 5) {
                K1(1, count - 2) = element;
            }
            if (count >= 5 && count < 8) {
                K2(0, count - 5) = element;
            }
            if (count >= 8) {
                K2(1, count - 7) = element;
            }
            count++;
        }

        infile.close();

        return count == 10;
    }

    //// Read a calibration matrix, of the USAC format.
    //// \param[in] path: path to the calibration matrix.
    //// \param[out] K1: calibration matrix.
    //// \return True if 5 values where extracted, False otherwise.
    bool loadCalibration(const char *path, libNumerics::matrix<double> &K) {
        std::ifstream infile(path);

        if (!infile.is_open())
            return false;

        double element;

        K(1, 0) = 0;
        K(2, 0) = 0;
        K(2, 1) = 0;
        K(2, 2) = 1;

        int count = 0;
        while (infile >> element) {
            if (count < 3) {
                K(0, count) = element;
            }
            if (count >= 3 && count < 5) {
                K(1, count - 2) = element;
            }
            count++;
        }

        infile.close();

        return count == 5;
    }

} // namespace utility
