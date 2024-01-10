//
// Created by clementriu on 5/18/20.
//

#include "artificial_output.hpp"

#include <fstream>
#include <iostream>

namespace utility {

    //// Save data to a text file for visualisation.
    //// \param fileName Emplacement of the file.
    //// \param points The dataset.
    //// \param vecInliers Index of the inliers.
    //// \param model Parameters of the estimated model.
    //// \param width, height Range of the dataset.
    void saveArtificialData(const char *fileName,
                            const libNumerics::matrix<double> &points,
                            const std::vector<int> &vecInliers,
                            const libNumerics::matrix<double> &model,
                            const int width, const int height,
                            bool verbose) {
        if (verbose) {
            std::cout << "Saving data to " << fileName << std::endl;
        }

        std::ofstream dataSave;
        dataSave.open(fileName);

        dataSave << width << " " << height << "\n";

        dataSave << model << "\n";

        std::vector<int>::const_iterator it = vecInliers.begin();
        for (; it != vecInliers.end(); ++it) {
            dataSave << *it << " ";
        }
        dataSave << "\n";

        dataSave << points << "\n";

        dataSave.close();

        if (verbose) {
            std::cout << "Data saved." << std::endl;
        }
    }

} // namespace utility
