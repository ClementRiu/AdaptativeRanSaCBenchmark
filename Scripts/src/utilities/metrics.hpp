//
// Created by clementriu on 8/31/20.
//

#ifndef MMM_ORSA_METRICS_HPP
#define MMM_ORSA_METRICS_HPP

#include <vector>

#include "libOrsa/model_estimator.hpp"
#include "libOrsa/libNumerics/matrix.h"


namespace utility {

    /// Compute number of true positives from a list of estimated positive indexes and a label vector of 1 for relevant elements and 0 for irrelevant elements.
    /// \param[in] positiveIdx: Vector of indexes of the elements estimated positives.
    /// \param[in] labels: Vector of labels of all data.
    /// Must be 1 for a relevant element and 0 for an irrelevant one.
    /// \return The number of true positives.
    int computeTruePositive(const std::vector<int> &positiveIdx, const std::vector<int> &labels);

    /// Compute the precision: tp / (tp + fp)
    /// \param[in] numTruePositive: Number of true positives.
    /// \param[in] numEstimatedPositive: Number of estimated positives.
    /// \return The precision.
    double computePrecision(int numTruePositive, int numEstimatedPositive);

    /// Compute the recall. tp / (tp + fn)
    /// \param[in] numTruePositive: Number of true positives.
    /// \param[in] numEstimatedPositive: Number of true inliers.
    /// \return The recall.
    double computeRecall(int numTruePositive, int numInlier);

    // TODO
    void computeModelError(const std::vector<int> &indexesToCompute,
                           const orsa::ModelEstimator *model,
                           const orsa::ModelEstimator::Model &modelParams,
                           std::vector<double> &errors);

   // TODO
    void computeModelError(const orsa::ModelEstimator *model,
                           const orsa::ModelEstimator::Model &modelParams,
                           std::vector<double> &errors);

} // namespace utility

#endif //MMM_ORSA_METRICS_HPP
