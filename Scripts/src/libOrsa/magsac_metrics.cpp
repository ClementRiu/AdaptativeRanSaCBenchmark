//
// Created by riuclement on 12/8/21.
//

#include "magsac_metrics.hpp"

#include <algorithm>

namespace orsa {

void createProposedInliers(int numPoints,
                           const std::vector<int> &estimatedInliers,
                           std::vector<bool> &proposedInliers) {
    proposedInliers.clear();
    for (int i = 0; i < numPoints; i++) {
        proposedInliers.push_back(false);
    }

    std::vector<int>::const_iterator itEstimatedInliers = estimatedInliers.begin();
    for (; itEstimatedInliers != estimatedInliers.end(); itEstimatedInliers++) {
        proposedInliers[*itEstimatedInliers] = true;
    }
}

bool isEstimatedInliers(const size_t index, const double threshold,
                        const std::vector<double> &errors,
                        const std::vector<bool> &proposedInliers) {
    return errors[index] < threshold && proposedInliers[index];
}

bool isTruePositive(size_t index, double threshold,
                    const std::vector<double> &errors,
                    const std::vector<bool> &groundTruth,
                    const std::vector<bool> &possibleInliers) {
    return isEstimatedInliers(index, threshold, errors, possibleInliers) && groundTruth[index];
}

void computeAllMetrics(const int numPoints, const int numGTInliers,
                       const std::vector<bool> &groundTruth,
                       const std::vector<double> &errors,
                       const std::vector<bool> &possibleInliers,
                       std::vector<double> &thresholds,
                       std::vector<double> &precisions,
                       std::vector<double> &recalls) {
    int numTruePositives = 0;
    int numEstimatedPositives = 0;

    std::vector<size_t> sortedIndexes = argsort(errors);

    thresholds.clear();
    precisions.clear();
    recalls.clear();
    for (int i = 0; i < numPoints; i++) {
        thresholds.push_back(0);
        precisions.push_back(0);
        recalls.push_back(0);
    }

    std::vector<size_t>::const_iterator itSortedIndex = sortedIndexes.begin();
    for (; itSortedIndex != sortedIndexes.end(); itSortedIndex++) {
        double threshold = errors[*itSortedIndex];
        bool estimatedInlier = isEstimatedInliers(*itSortedIndex, threshold, errors, possibleInliers);
        numTruePositives += (int) (estimatedInlier && groundTruth[*itSortedIndex]);
        numEstimatedPositives += (int) estimatedInlier;

        thresholds[*itSortedIndex] = threshold;
        if (numEstimatedPositives > 0) {
            precisions[*itSortedIndex] = numTruePositives / numEstimatedPositives;
        }
        if (numGTInliers > 0) {
            recalls[*itSortedIndex] = numTruePositives / numGTInliers;
        }
    }
}

void getMetricAFromMetricB(const int numPoints,
                           const double recallToMatch,
                           const std::vector<double> &thresholds,
                           const std::vector<double> &metricsA,
                           const std::vector<double> &metricsB,
                           int &bestThreshold, int &bestMetricA, int &bestMetricB) {
    int bestIndex = -1;
    bestMetricA = -1;
    for (int i = 0; i < numPoints; i++) {
        if (metricsB[i] >= recallToMatch && metricsA[i] > bestMetricA) {
            bestIndex = i;
            bestMetricA = metricsA[bestIndex];
        }
    }
    bestMetricB = metricsB[bestIndex];
    bestThreshold = thresholds[bestIndex];
}

void getWeightedMetric(const int numGTInliers,
                       const std::vector<int> &estimatedInliers,
                       const std::vector<int> &groundTruthLabels,
                       const std::vector<double> &weights,
                       double &weightedPrecision, double &weightedRecall) {
    double numerator = 0;
    double denominatorPrecision = 0;
    std::vector<int>::const_iterator itEstimatedInliers = estimatedInliers.begin();
    for (; itEstimatedInliers != estimatedInliers.end(); itEstimatedInliers++) {
        double weight = weights[*itEstimatedInliers];
        int label = groundTruthLabels[*itEstimatedInliers];
        numerator += label * weight;
        denominatorPrecision += weight;
    }
    double maxWeight = -1;
    std::vector<double>::const_iterator itWeights = weights.begin();
    for (; itWeights != weights.end(); itWeights++) {
        if (*itWeights > maxWeight) maxWeight = *itWeights;
    }
    double denominatorRecall = numGTInliers * maxWeight;

    if (denominatorPrecision > 0) {
        weightedPrecision = numerator / denominatorPrecision;
    }
    if (denominatorRecall > 0) {
        weightedRecall = numerator / denominatorRecall;
    }
}

} // Namespace orsa