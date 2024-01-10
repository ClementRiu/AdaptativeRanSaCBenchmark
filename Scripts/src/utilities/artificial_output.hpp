//
// Created by clementriu on 5/18/20.
//

#ifndef MMM_ORSA_ARTIFICIAL_OUTPUT_HPP
#define MMM_ORSA_ARTIFICIAL_OUTPUT_HPP

#include <vector>

#include "libOrsa/libNumerics/matrix.h"


namespace utility {

    //// Save data to a text file for visualisation.
    //// \param fileName Emplacement of the file.
    //// \param points The dataset.
    //// \param vecInliers Index of the inliers.
    //// \param model Parameters of the estimated model.
    //// \param width, height Range of the dataset.
    void
    saveArtificialData(const char *fileName, const libNumerics::matrix<double> &points,
                       const std::vector<int> &vecInliers,
                       const libNumerics::matrix<double> &model, int width, int height, bool verbose = false);


} // namespace utility
#endif //MMM_ORSA_ARTIFICIAL_OUTPUT_HPP
