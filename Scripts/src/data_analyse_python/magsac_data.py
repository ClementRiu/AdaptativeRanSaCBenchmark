import matplotlib.pyplot as plt
import numpy as np


class Magsac_Data:
    colors = ['b', 'r', 'k']
    def __init__(self, labels, errors, prop_inliers, weights):
        self._labels = labels
        self._errors = errors
        self._prop_inliers = prop_inliers
        self._weights = weights

        self._bad_data = False
        if not isinstance(labels, np.ndarray):
            if labels == None:
                self._bad_data = True

        self._all_computed = False

    def get_GT_inliers(self, threshold = 0):
        return np.where(self._labels)[0]

    def get_estimated_inliers(self, threshold):
        return np.intersect1d(np.where(self._errors <= threshold)[0], self._prop_inliers)

    def get_true_positive(self, threshold):
        return np.intersect1d(self.get_GT_inliers(), self.get_estimated_inliers(threshold))

    def is_GT_inliers(self, index, threshold = 0):
        return bool(self._labels[index])

    def is_estimated_inliers(self, index, threshold):
        return self._errors[index] <= threshold and index in self._prop_inliers

    def is_true_positive(self, index, threshold):
        return self.is_estimated_inliers(index, threshold) and self.is_GT_inliers(index)

    def get_threshold_from_metric(self, metric_value, metric_type):
        print("DEPRECATED")
        computed_metric = metric_value
        old_metric = metric_value

        threshold = 0
        saved_threshold = []

        num_true_positives = 0
        denominator = 0

        errors_index_sorted = np.argsort(self._errors)

        for index in errors_index_sorted:
            if ((old_metric - metric_value) * (computed_metric - metric_value) < 0):
                saved_threshold.append(threshold)

            threshold = self._errors[index]
            num_true_positives += self.is_true_positive(index, threshold)
            if metric_type == 'precision':
                denominator += self.is_estimated_inliers(index, threshold)
            if metric_type == 'recall':
                number_GT_inliers = len(self.get_GT_inliers())
                if number_GT_inliers > 0:
                    denominator = number_GT_inliers

            if denominator > 0:
                old_metric = computed_metric
                computed_metric = num_true_positives / denominator
        return saved_threshold

    def _compute_all_threshold_metric(self):
        if self._all_computed:
            return None

        self._all_computed = True

        compute_recall = False
        number_GT_inliers = len(self.get_GT_inliers())
        if number_GT_inliers > 0:
            compute_recall = True

        num_true_positives = 0
        num_estimated_positives = 0

        errors_index_sorted = np.argsort(self._errors)

        num_pts = len(self._labels)
        self._thresholds = np.zeros(num_pts)
        self._precisions = np.zeros(num_pts)
        self._recalls = np.zeros(num_pts)

        for index in errors_index_sorted:
            threshold = self._errors[index]
            num_true_positives += self.is_true_positive(index, threshold)
            num_estimated_positives += self.is_estimated_inliers(index, threshold)

            self._thresholds[index] = threshold
            if num_estimated_positives > 0:
                self._precisions[index] = num_true_positives / num_estimated_positives
            if compute_recall:
                self._recalls[index] = num_true_positives / number_GT_inliers

        return None

    def _get_metrics_from_metric(self, value_required, values, values_other):
        possible_indexes = np.where(values >= value_required)[0]
        if len(possible_indexes) == 0:
            return -1, 0, 0
        possible_other = values_other[possible_indexes]
        best_index = possible_indexes[np.argmax(possible_other)]

        best_threshold = self._thresholds[best_index]
        best_precision = self._precisions[best_index]
        best_recall = self._recalls[best_index]

        return best_threshold, best_precision, best_recall

    def get_metrics_from_name(self, value_requiered, metric_name):
        if self._bad_data:
            return -1, -1, -1
        if metric_name == 'p':
            return self.get_metrics_from_precision(value_requiered)
        if metric_name == 'r':
            return self.get_metrics_from_recall(value_requiered)
        return -1, -1, -1

    def get_metrics_from_precision(self, precision_requiered):
        if self._bad_data:
            return -1, -1, -1
        self._compute_all_threshold_metric()
        return self._get_metrics_from_metric(precision_requiered, self._precisions, self._recalls)

    def get_metrics_from_recall(self, recall_requiered):
        if self._bad_data:
            return -1, -1, -1
        self._compute_all_threshold_metric()
        return self._get_metrics_from_metric(recall_requiered, self._recalls, self._precisions)


    def get_precision_from_threshold(self, threshold):
        if self._bad_data:
            return -1
        number_estimated_inliers = len(self.get_estimated_inliers(threshold))
        if number_estimated_inliers > 0:
            return len(self.get_true_positive(threshold)) / number_estimated_inliers
        return 0

    def get_recall_from_threshold(self, threshold):
        if self._bad_data:
            return -1
        number_GT_inliers = len(self.get_GT_inliers())
        if number_GT_inliers > 0:
            return len(self.get_true_positive(threshold)) / number_GT_inliers
        return 0

    def get_weighted_metrics(self):
        if self._bad_data:
            return -1, -1
        numerator = 0
        denominator_precision = 0
        for index_inliers, index_global in enumerate(self._prop_inliers):
            weight = self._weights[index_inliers]
            label = self._labels[index_global]
            numerator += weight * label
            denominator_precision += weight
        denominator_recall = np.sum(self._labels) * np.ceil(self._weights.max())
        if denominator_precision > 0:
            precision = numerator / denominator_precision
        else:
            precision = 0
        if denominator_recall > 0:
            recall = numerator / denominator_recall
        else:
            recall = 0
        return precision, recall

    def show_weights_errors(self):
        points_col = []
        for index in self._prop_inliers:
            points_col.append(self.colors[self._labels[index]])

        plt.figure(figsize=(20, 20))
        plt.scatter(self._weights, self._errors[self._prop_inliers], c=np.array(points_col), alpha=0.2)
        plt.xlabel("Weights")
        plt.ylabel("Errors")
        plt.show()
