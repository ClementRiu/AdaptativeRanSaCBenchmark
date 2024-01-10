//
// Created by clementriu on 8/28/20.
//

#ifndef MMM_ORSA_DATA_HANDLER_HPP
#define MMM_ORSA_DATA_HANDLER_HPP

#include <vector>

#include "libOrsa/match.hpp"


namespace utility {

    /// Randomly merges two vectors and save the origin of each entry.
    /// \param[in] vect0: First vector to merge, associated with label0.
    /// \param[in] vect1: Second vector to merge, associated with label1.
    /// \param[out] mixedVector: Resulting vector.
    /// \param[in] label0: Label to give elements from the first vector.
    /// \param[in] label1: Label to give elements from the second vector.
    /// \param[out] labels: Resulting label vector.
    /// If value is label0, the corresponding element is from vect0.
    /// If value is label1, the corresponding element is from vect1.
    template<typename T, typename U>
    void
    randomDataFusionWithMemory(const std::vector<T> &vect0, const std::vector<T> &vect1, std::vector<T> &mixedVector,
                               const U label0, const U label1, std::vector<U> &labels) {
        typename std::vector<T>::const_iterator it0 = vect0.begin();
        typename std::vector<T>::const_iterator it1 = vect1.begin();
        size_t size0 = vect0.size();
        size_t size1 = vect1.size();
        while (mixedVector.size() < size0 + size1) {
            if (it0 != vect0.end() && it1 != vect1.end()) {
                int chosenDataSet = std::rand() % 2;
                if (chosenDataSet == 0) {
                    mixedVector.push_back(*it0);
                    labels.push_back(label0);
                    it0++;
                } else {
                    mixedVector.push_back(*it1);
                    labels.push_back(label1);
                    it1++;
                }
            } else {
                if (it0 == vect0.end() && it1 != vect1.end()) {
                    mixedVector.push_back(*it1);
                    labels.push_back(label1);
                    it1++;
                } else {
                    if (it1 == vect1.end() && it0 != vect0.end()) {
                        mixedVector.push_back(*it0);
                        labels.push_back(label0);
                        it0++;
                    } else {
                        throw std::runtime_error("Problem during fusion of dataset.");
                    }
                }
            }
        }

    }


    //// Compute mean of a vect.
    //// \param vect [in] vector on which to compute the mean.
    template<typename T>
    double meanOfVect(const std::vector<T> &vect) {
        double mean = 0;
        typename std::vector<T>::const_iterator it = vect.begin();
        for (; it != vect.end(); it++) {
            mean += *it;
        }
        return mean / vect.size();
    }

    //// Compute linewise mean of a vect of vect.
    //// \param vect [in] vector of vectors on which to compute the mean.
    //// \param meanVect [out] vector of the computed means.
    template<typename T>
    void meanOfVect(const std::vector<std::vector<T>> &vect, std::vector<double> &meanVect) {
        typename std::vector<std::vector<T>>::const_iterator itVect = vect.begin();
        for (; itVect != vect.end(); itVect++) {
            double mean = 0;
            typename std::vector<T>::const_iterator it = (*itVect).begin();
            for (; it != (*itVect).end(); it++) {
                mean += *it;
            }
            meanVect.push_back(mean / (*itVect).size());
        }
    }

    //// Compute STD of a vect if mean is given.
    //// \param mean [in] mean of the vect.
    //// \param vect [in] vector on which to compute the std.
    template<typename T>
    double standardDeviation(const double mean, const std::vector<T> &vect) {
        double var = 0;
        typename std::vector<T>::const_iterator it = vect.begin();
        for (; it != vect.end(); it++) {
            var += (*it - mean) * (*it - mean);
        }
        return std::sqrt(var / vect.size());
    }

    //// Compute linewise STD of a vect of vect.
    //// \param meanVect [in] vector of the vectors' means.
    //// \param vect [in] vector of vectors on which to compute the STD.
    //// \param meanVect [out] vector of the computed STDs.
    template<typename T>
    void standardDeviation(const std::vector<double> &meanVect, const std::vector<std::vector<T>> &vect,
                           std::vector<double> &stdVect) {
        for (size_t i = 0; i < vect.size(); i++) {
            std::vector<T> locVect = vect[i];
            double locMean = meanVect[i];
            double var = 0;
            typename std::vector<T>::const_iterator it = locVect.begin();
            for (; it != locVect.end(); it++) {
                var += (*it - locMean) * (*it - locMean);
            }
            stdVect.push_back(std::sqrt(var / locVect.size()));
        }
    }

} // namespace utility


#endif //MMM_ORSA_DATA_HANDLER_HPP
