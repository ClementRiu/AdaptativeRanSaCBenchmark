//
// Created by clementriu on 6/1/20.
//

#ifndef MMM_ORSA_USAC_FILE_HANDLER_HPP
#define MMM_ORSA_USAC_FILE_HANDLER_HPP

#include <cmath>
#include <vector>

#include "libOrsa/match.hpp"


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
    bool saveExpInfo(const char *nameFile, unsigned int seed, int nGen, int nRun, const char *pathToIn,
                     const std::vector<const char *> &algorithmNames,
                     const std::vector<double> &precisionMeanVect, const std::vector<double> &precisionStdVect,
                     const std::vector<std::vector<double>> &precisionVectVect,
                     const std::vector<double> &recallMeanVect, const std::vector<double> &recallStdVect,
                     const std::vector<std::vector<double>> &recallVectVect,
                     const std::vector<double> &runtimeMeanVect, const std::vector<double> &runtimeStdVect,
                     const std::vector<std::vector<double>> &runtimeVectVect,
                     const std::vector<std::vector<double>> &sigmaEstimatedVectVect,
                     bool fullSave = false);

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
    bool saveExpInfo(const char *nameFile, unsigned int seed, int nRun, const char *pathToIn,
                     const std::vector<const char *> &algorithmNames,
                     const std::vector<double> &TMeanVect, const std::vector<double> &TStdVect,
                     const std::vector<std::vector<int>> &TVectVect,
                     const std::vector<double> &vpmMeanVect, const std::vector<double> &vpmStdVect,
                     const std::vector<std::vector<double>> &vpmVectVect,
                     const std::vector<double> &meanRuntime,
                     bool fullSave);

    //// Save matches to a file.
    //// \param pathToMatches: Emplacement of the file.
    //// \param matchesVect: Vector of vectors of matches.
    template <typename M>
    bool saveMatches(const char *pathToMatches, const std::vector<std::vector<M>> &matchesVect) {
        std::ofstream f(pathToMatches);
        if (f.is_open()) {
            typename std::vector<std::vector<M>>::const_iterator matchesVectIt = matchesVect.begin();
            for (; matchesVectIt != matchesVect.end(); ++matchesVectIt) {
                typename std::vector<M>::const_iterator matchIt = (*matchesVectIt).begin();
                for (; matchIt != (*matchesVectIt).end(); ++matchIt) {
                    f << *matchIt;
                }
                f << "\n";
            }
        }
        return f.is_open();
    }

    //// Read a calibration matrix, of the USAC format.
    //// \param[in] path: path to the calibration matrix.
    //// \param[out] K1, K2: both calibration matrix.
    //// \return True if 10 values where extracted, False otherwise.
    bool loadCalibration(const char *path, libNumerics::matrix<double> &K1, libNumerics::matrix<double> &K2);

    //// Read a calibration matrix, of the USAC format.
    //// \param[in] path: path to the calibration matrix.
    //// \param[out] K1: calibration matrix.
    //// \return True if 5 values where extracted, False otherwise.
    bool loadCalibration(const char *path, libNumerics::matrix<double> &K1);
} // namespace utility

#endif //MMM_ORSA_USAC_FILE_HANDLER_HPP
