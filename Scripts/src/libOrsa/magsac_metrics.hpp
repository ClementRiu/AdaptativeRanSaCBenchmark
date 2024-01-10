//
// Created by riuclement on 12/8/21.
//

#ifndef MMM_ORSA_MAGSAC_METRICS_H
#define MMM_ORSA_MAGSAC_METRICS_H

#include <vector>
#include <algorithm>
#include <numeric>

namespace orsa {

template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] < array[right];
              });

    return indices;
}

void createProposedInliers(int numPoints,
                           const std::vector<int> &estimatedInliers,
                           std::vector<bool> &proposedInliers);

bool isEstimatedInliers(size_t index, double threshold,
                        const std::vector<double> &errors,
                        const std::vector<bool> &possibleInliers);

bool isTruePositive(size_t index, double threshold,
                    const std::vector<double> &errors,
                    const std::vector<bool> &groundTruth,
                    const std::vector<bool> &possibleInliers);

void computeAllMetrics(int numPoints, int numGTInliers,
                       const std::vector<bool> &groundTruth,
                       const std::vector<double> &errors,
                       const std::vector<bool> &possibleInliers,
                       std::vector<double> &thresholds,
                       std::vector<double> &precisions,
                       std::vector<double> &recalls);

void getMetricAFromMetricB(int numPoints,
                           double recallToMatch,
                           const std::vector<double> &thresholds,
                           const std::vector<double> &metricsA,
                           const std::vector<double> &metricsB,
                           int &bestThreshold, int &bestMetricA, int &bestMetricB);

void getWeightedMetric(const int numGTInliers,
                       const std::vector<int> &estimatedInliers,
                       const std::vector<int> &groundTruthLabels,
                       const std::vector<double> &weights,
                       double &weightedPrecision, double &weightedRecall);

} // Namespace orsa

#endif //MMM_ORSA_MAGSAC_METRICS_H
