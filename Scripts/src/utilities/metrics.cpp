//
// Created by clementriu on 8/31/20.
//

#include "metrics.hpp"


namespace utility {

    /// Compute number of true positives from a list of estimated positive indexes and a label vector of 1 for relevant elements and 0 for irrelevant elements.
    /// \param[in] positiveIdx: Vector of indexes of the elements estimated positives.
    /// \param[in] labels: Vector of labels of all data.
    /// Must be 1 for a relevant element and 0 for an irrelevant one.
    /// \return The number of true positives.
    int computeTruePositive(const std::vector<int> &estimatedPositiveIdx, const std::vector<int> &labels) {
        int numTruePositives = 0;
        std::vector<int>::const_iterator itEstimatedPositiveIdx = estimatedPositiveIdx.begin();
        for (; itEstimatedPositiveIdx != estimatedPositiveIdx.end(); itEstimatedPositiveIdx++) {
            numTruePositives += labels[*itEstimatedPositiveIdx];
        }
        return numTruePositives;
    }

    /// Compute the precision: tp / (tp + fp)
    /// \param[in] numTruePositive: Number of true positives.
    /// \param[in] numEstimatedPositive: Number of estimated positives.
    /// \return The precision.
    double computePrecision(const int numTruePositive, const int numEstimatedPositive) {
        return numTruePositive / (double) numEstimatedPositive;
    }

    /// Compute the recall. tp / (tp + fn)
    /// \param[in] numTruePositive: Number of true positives.
    /// \param[in] numEstimatedPositive: Number of true inliers.
    /// \return The recall.
    double computeRecall(const int numTruePositive, const int numInlier) {
        return numTruePositive / (double) numInlier;
    }

    void computeModelError(const std::vector<int> &indexesToCompute,
                           const orsa::ModelEstimator *model,
                           const orsa::ModelEstimator::Model &modelParams,
                           std::vector<double> &errors) {
        std::vector<int>::const_iterator itIndexesToCompute = indexesToCompute.begin();
        for (; itIndexesToCompute != indexesToCompute.end(); itIndexesToCompute++){
            errors.push_back(std::sqrt(model->Error(modelParams, *itIndexesToCompute)));
        }
    }

    void computeModelError(const orsa::ModelEstimator *model,
                           const orsa::ModelEstimator::Model &modelParams,
                           std::vector<double> &errors) {
        for (size_t i = 0; i < model->NbData(); i++){
            errors.push_back(std::sqrt(model->Error(modelParams, i)));
        }
    }

} // namespace utility