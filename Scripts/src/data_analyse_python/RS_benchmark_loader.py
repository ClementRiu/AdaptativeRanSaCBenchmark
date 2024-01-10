import glob
import matplotlib.pyplot as plt
import numpy as np

from magsac_data import Magsac_Data


class RS_Data_Loader:
    def __init__ (self,
                  default_file_paths_,
                  magsac_file_paths_,
                  magsac_file_names_,
                  algorithm_names_,
                  algorithm_names_and_number_,
                  metric_names_,
                  dataset_size_,
                  sigma_values_,
                  ratio_values_,
                  num_runs_,
                  dataset_name_):
        self._default_file_paths = default_file_paths_
        self._magsac_file_paths = magsac_file_paths_
        self._magsac_file_names = magsac_file_names_

        self._algorithm_names = algorithm_names_
        self._algorithm_names_and_number = algorithm_names_and_number_
        self._metric_names = metric_names_
        self._dataset_size = dataset_size_

        self._num_metrics = len(self._metric_names)

        self._sigma_values = sigma_values_
        self._ratio_values = ratio_values_
        self._num_runs = num_runs_

        self._all_values = np.empty((len(self._algorithm_names),
                                     self._num_metrics,
                                     self._dataset_size,
                                     len(self._sigma_values),
                                     len(self._ratio_values),
                                     self._num_runs,
                                     ))
        self._all_values[:] = np.NaN

        if len(self._magsac_file_paths) > 0:
            self._magsac_data_raw = {dataset_index: {sigma_value: {ratio_value: []
                                              for ratio_value in self._ratio_values}
                                  for sigma_value in self._sigma_values}
                    for dataset_index in range(self._dataset_size)}

        self._dataset_name = dataset_name_


    def _format_value_for_filename(self, numerical_value_):
        if numerical_value_ == 0:
            return "0.0"
        else:
            formated_value = "%.1f" % numerical_value_
            if formated_value.startswith("0."):
                return formated_value[1:]
            return formated_value

    def _parse_values(self, values_):
        parsed_values = np.empty((len(self._algorithm_names), len(values_), self._num_runs))
        parsed_values[:] = np.NaN
        for index_metric, value_metric in enumerate(values_):
            for index_algorithm, value_algorithm in enumerate(value_metric):
                for index_run, value_run in enumerate(value_algorithm.rstrip().split()):
                    parsed_values[index_algorithm, index_metric, index_run] = value_run
        return parsed_values


    def _open_and_parse(self, path_):
        metrics = [[] for i in range(self._num_metrics)]
        with open(path_, 'r') as file:
            line = file.readline()
            while line:
                if 'p:' in line:
                    line = file.readline()
                    metrics[self._metric_names['p']].append(line)
                elif 'r:' in line:
                    line = file.readline()
                    metrics[self._metric_names['r']].append(line)
                elif 's:' in line:
                    line = file.readline()
                    metrics[self._metric_names['s']].append(line)
                elif 't:' in line:
                    line = file.readline()
                    metrics[self._metric_names['t']].append(line)
                else:
                    line = file.readline()
        return self._parse_values(metrics)

    def _load_normal_data(self):
        for dataset_index in range(self._dataset_size):
            for sigma_index, sigma_value in enumerate(self._sigma_values):
                for ratio_index, ratio_value in enumerate(self._ratio_values):
                    path = glob.glob(self._default_file_paths.format(dataset_index + 1,
                                                         self._format_value_for_filename(sigma_value),
                                                         self._format_value_for_filename(ratio_value))
                                     )
                    if len(path) == 1:
                        parsed_values = self._open_and_parse(path[0])
                        self._all_values[:, :, dataset_index, sigma_index, ratio_index, :] = parsed_values

        return None

    def _open_magsac(self, path_, type_used_=float):
        data = []
        with open(path_, "r") as file:
            for line in file:
                element = np.array([value for value in line.split(", ")[:-1]]).astype(type_used_)
                data.append(element)
        return data

    def _load_magsac_data(self):
        for dataset_index in range(self._dataset_size):
            for sigma_value in self._sigma_values:
                for ratio_value in self._ratio_values:
                    magsac_data = [None] * len(self._magsac_file_names)
                    for file_index, (file_name, file_type) in enumerate(self._magsac_file_names):
                        path = glob.glob(self._magsac_file_paths.format(dataset_index + 1,
                                                              self._format_value_for_filename(sigma_value),
                                                              self._format_value_for_filename(ratio_value),
                                                              file_name)
                                        )
                        if len(path) == 1:
                            magsac_data[file_index] = self._open_magsac(path[0], file_type)
                    if magsac_data[0] != None:
                        data_to_analyse = [Magsac_Data(
                            magsac_data[1][i // 5], ## Labels
                            magsac_data[3][i],      ## Full errors
                            magsac_data[4][i],      ## Inliers index
                            magsac_data[0][i])      ## Weights
                                for i in range(len(magsac_data[0]))]
                    else:
                        data_to_analyse = [Magsac_Data(None, None, None, None)]
                    self._magsac_data_raw[dataset_index][sigma_value][ratio_value] = data_to_analyse

        return None

    def _threshold_values_and_mean(self, threshold_ = 0):
        self._all_values_thresholded = self._all_values.copy()

        bad_data = self._all_values_thresholded <= threshold_
        self._all_values_thresholded[bad_data] = np.NaN
        self._all_values_thresholded[:, self._metric_names['t'], :, :, :, :] = self._all_values[:, self._metric_names['t'], :, :, :, :]

        self._all_values_mean = np.nanmean(self._all_values_thresholded, axis=5)
        self._all_values_mean = np.nan_to_num(self._all_values_mean, nan=-1)
        self._other_failed = np.sum(bad_data, axis=5)

        return None

    def _true_threshold(self, dataset_index_, sigma_value_, ratio_value_, sigma_index_, ratio_index_):
        return sigma_value_

    def _orsa_threshold(self, dataset_index_, sigma_value_, ratio_value_, sigma_index_, ratio_index_):
        return self._all_values_mean[self._algorithm_names_and_number["AC-RanSaC"], self._metric_names['s'], dataset_index, sigma_index, ratio_index]

    def _get_metrics_from_fixed_threshold(self, threshold_function_, algorithm_name_, threshold_):
        for dataset_index in range(self._dataset_size):
            for sigma_index, sigma_value in enumerate(self._sigma_values):
                for ratio_index, ratio_value in enumerate(self._ratio_values):
                    cutoff_threshold = threshold_function_(dataset_index, sigma_value, ratio_value, sigma_index, ratio_index)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['s'], dataset_index, sigma_index, ratio_index] = cutoff_threshold

                    values_to_transform = self._magsac_data_raw[dataset_index][sigma_value][ratio_value]

                    magsac_computed_p = np.array([value_to_transform.get_precision_from_threshold(cutoff_threshold)
                                         for value_to_transform in values_to_transform], dtype=np.float)
                    magsac_computed_p[magsac_computed_p <= threshold_] = np.NaN
                    magsac_computed_p = np.nan_to_num(np.nanmean(magsac_computed_p), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['p'], dataset_index, sigma_index, ratio_index] = magsac_computed_p

                    magsac_computed_r = np.array([value_to_transform.get_recall_from_threshold(cutoff_threshold)
                                         for value_to_transform in values_to_transform], dtype=np.float)
                    magsac_computed_r[magsac_computed_r <= threshold_] = np.NaN
                    magsac_computed_r = np.nan_to_num(np.nanmean(magsac_computed_r), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['r'], dataset_index, sigma_index, ratio_index] = magsac_computed_r

        return None

    def _get_metrics_from_other_algo(self, other_algo_name_, other_var_name_, algorithm_name_, threshold_):
        for dataset_index in range(self._dataset_size):
            for sigma_index, sigma_value in enumerate(self._sigma_values):
                for ratio_index, ratio_value in enumerate(self._ratio_values):
                    orsa_var = self._all_values_mean[self._algorithm_names_and_number[other_algo_name_], self._metric_names[other_var_name_], dataset_index, sigma_index, ratio_index]

                    values_to_transform = self._magsac_data_raw[dataset_index][sigma_value][ratio_value]

                    magsac_computed_p = []
                    magsac_computed_r = []
                    magsac_computed_s = []

                    for value_to_transform in values_to_transform:
                        try:
                            cutoff_threshold, precision, recall = value_to_transform.get_metrics_from_name(orsa_var, other_var_name_)
                        except IndexError:
                            cutoff_threshold, precision, recall = np.NaN, np.NaN, np.NaN

                        magsac_computed_s.append(cutoff_threshold)
                        magsac_computed_p.append(precision)
                        magsac_computed_r.append(recall)

                    magsac_computed_p = np.array(magsac_computed_p, dtype=np.float)
                    magsac_computed_r = np.array(magsac_computed_r, dtype=np.float)
                    magsac_computed_s = np.array(magsac_computed_s, dtype=np.float)

                    magsac_computed_s[magsac_computed_s <= threshold_] = np.NaN
                    magsac_computed_s = np.nan_to_num(np.nanmean(magsac_computed_s), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['s'], dataset_index, sigma_index, ratio_index] = magsac_computed_s

                    magsac_computed_p[magsac_computed_p <= threshold_] = np.NaN
                    magsac_computed_p = np.nan_to_num(np.nanmean(magsac_computed_p), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['p'], dataset_index, sigma_index, ratio_index] = magsac_computed_p

                    magsac_computed_r[magsac_computed_r <= threshold_] = np.NaN
                    magsac_computed_r = np.nan_to_num(np.nanmean(magsac_computed_r), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['r'], dataset_index, sigma_index, ratio_index] = magsac_computed_r

        return None

    def _get_weighted_metrics(self, algorithm_name_, threshold_):
        for dataset_index in range(self._dataset_size):
            for sigma_index, sigma_value in enumerate(self._sigma_values):
                for ratio_index, ratio_value in enumerate(self._ratio_values):
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['s'], dataset_index, sigma_index, ratio_index] = -1

                    values_to_transform = self._magsac_data_raw[dataset_index][sigma_value][ratio_value]

                    magsac_computed_prs = [value_to_transform.get_weighted_metrics()
                                         for value_to_transform in values_to_transform]


                    magsac_computed_p = np.array([magsac_computed_pr[0] for magsac_computed_pr in magsac_computed_prs], dtype=np.float)
                    magsac_computed_p[magsac_computed_p <= threshold_] = np.NaN
                    magsac_computed_p = np.nan_to_num(np.nanmean(magsac_computed_p), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['p'], dataset_index, sigma_index, ratio_index] = magsac_computed_p

                    magsac_computed_r = np.array([magsac_computed_pr[1] for magsac_computed_pr in magsac_computed_prs], dtype=np.float)
                    magsac_computed_r[magsac_computed_r <= threshold_] = np.NaN
                    magsac_computed_r = np.nan_to_num(np.nanmean(magsac_computed_r), nan=-1)
                    self._all_values_mean[self._algorithm_names_and_number[algorithm_name_], self._metric_names['r'], dataset_index, sigma_index, ratio_index] = magsac_computed_r

        return None

    def _get_magsac_metric(self, threshold_ = 0):
        if "MAGSAC-True-Threshold" in self._algorithm_names_and_number.keys():
            self._get_metrics_from_fixed_threshold(self._true_threshold, "MAGSAC-True-Threshold", threshold_)
        if "MAGSAC-AC-RanSaC-Threshold" in self._algorithm_names_and_number.keys() and "AC-RanSaC" in self._algorithm_names_and_number.keys():
            self._get_metrics_from_fixed_threshold(self._orsa_threshold, "MAGSAC-AC-RanSaC-Threshold", threshold_)
        if "MAGSAC-P" in self._algorithm_names_and_number.keys() and "AC-RanSaC" in self._algorithm_names_and_number.keys():
            self._get_metrics_from_other_algo("AC-RanSaC", 'p', "MAGSAC-P", threshold_)
        if "MAGSAC-R" in self._algorithm_names_and_number.keys() and "AC-RanSaC" in self._algorithm_names_and_number.keys():
            self._get_metrics_from_other_algo("AC-RanSaC", 'r', "MAGSAC-R", threshold_)
        if "MAGSAC-W" in self._algorithm_names_and_number.keys():
            self._get_weighted_metrics("MAGSAC-W", threshold_)

        return None

    def _compute_f1score(self):
        precision = self._all_values_mean[:, self._metric_names['p'], :, : , :]
        recall = self._all_values_mean[:, self._metric_names['r'], :, : , :]
        self._f1_score = (2 * precision * recall) / (precision + recall)

        return None

    def load_data(self):
        self._load_normal_data()
        print("Normal data loaded.")

        if len(self._magsac_file_paths) > 0:
            self._load_magsac_data()
            print("MAGSAC data loaded.")

        return None

    def compute_metrics(self, threshold_ = 0.0):
        self._threshold_values_and_mean(threshold_)
        print("Normal metrics computed.")

        self._get_magsac_metric(threshold_)
        print("MAGSAC metrics computed.")

        self._all_values_mean[self._all_values_mean < 0] = np.NaN

        self._compute_f1score()
        print("F1-Score computed.")

        return None

    def generate_single_image(self,
                              algo_to_remove_,
                              metric_to_print_,
                              dataset_index_,
                              ratio_index_,
                              value_printed_,
                              display_style_,
                              save_=False,
                              threshold_=None):

        legend = "{}."
        xlabel_n = "Noise std value (in pixel)"
        xlabel_o = "Outlier ratio value"
        generic_title = "{} for outlier ratio of {} on dataset {}."
        generic_suptitle = "{} " + "for outlier ratio of {} on artificial {} - dataset {}.".format(self._ratio_values[ratio_index_],
                    self._dataset_name,
                    dataset_index_)

        fig, axs = plt.subplots(1, 1, figsize=(8, 4), sharey=True)
        fig.suptitle(generic_suptitle.format(value_printed_))
        axs.set_xlabel(xlabel_n)
        axs.set_ylabel(value_printed_)
        if threshold_ and threshold_ >= 1:
            axs.set_ylim(-1.1, threshold_)
        if threshold_ and threshold_ < 1:
            axs.set_ylim(threshold_, 1.05)
        lines = []
        legends = []

        for algorithm_index, algorithm_name in enumerate(self._algorithm_names):
            if algorithm_name not in algo_to_remove_:
                if metric_to_print_ != "f1":
                    line, = axs.plot(
                        self._sigma_values,
                        self._all_values_mean[algorithm_index,
                                              self._metric_names[metric_to_print_],
                                              dataset_index_,
                                              :,
                                              ratio_index_],
                        color=display_style_["colors"][algorithm_name],
                        linestyle=display_style_["linestyles"][algorithm_name],
                        alpha=display_style_["alpha"])
                    lines.append(line)
                else:
                    line, = axs.plot(
                        self._sigma_values,
                        self._f1_score[algorithm_index,
                                       dataset_index_,
                                       :,
                                       ratio_index_],
                        color=display_style_["colors"][algorithm_name],
                        linestyle=display_style_["linestyles"][algorithm_name],
                        alpha=display_style_["alpha"])
                    lines.append(line)
                legends.append(legend.format(algorithm_name))

        axs.legend(lines, legends)
        if save_:
            plt.savefig(save_, bbox_inches='tight')
        else:
            plt.show()

    def generate_full_image(self,
                              algo_to_remove_,
                              metric_to_print_,
                              value_printed_,
                              display_style_,
                              save_=False,
                              threshold_=None):

        legend = "{}."
        xlabel_n = "Noise std value (in pixel)"
        xlabel_o = "Outlier ratio value"
        generic_title = "{} for outlier ratio of {} on dataset {}."
        generic_suptitle = "{} of " + "{} for only inlier, with noise, artificial {}.".format(" ".join(self._algorithm_names), self._dataset_name)

        fig, axs = plt.subplots(self._dataset_size, len(self._ratio_values), figsize=(10 * len(self._ratio_values), 10 * self._dataset_size), sharey=True)
        fig.suptitle(generic_suptitle.format(value_printed_))

        for dataset_index in range(self._dataset_size):
            for ratio_index, ratio_value in enumerate(self._ratio_values):
                axs[dataset_index, ratio_index].set_title(generic_title.format(value_printed_, ratio_value, dataset_index))
                axs[dataset_index, ratio_index].set_xlabel(xlabel_n)
                axs[dataset_index, ratio_index].set_ylabel(value_printed_)
                if threshold_ and threshold_ >= 1:
                    axs[dataset_index, ratio_index].set_ylim(-1.1, threshold_)
                if threshold_ and threshold_ < 1:
                        axs[dataset_index, ratio_index].set_ylim(threshold_, 1.05)
                lines = []
                legends = []

                for algorithm_index, algorithm_name in enumerate(self._algorithm_names):
                    if algorithm_name not in algo_to_remove_:
                        if metric_to_print_ != "f1":
                            line, = axs[dataset_index, ratio_index].plot(
                                self._sigma_values,
                                self._all_values_mean[algorithm_index,
                                                      self._metric_names[metric_to_print_],
                                                      dataset_index,
                                                      :,
                                                      ratio_index],
                                color=display_style_["colors"][algorithm_name],
                                linestyle=display_style_["linestyles"][algorithm_name],
                                alpha=display_style_["alpha"])
                            lines.append(line)
                        else:
                            line, = axs[dataset_index, ratio_index].plot(
                                self._sigma_values,
                                self._f1_score[algorithm_index,
                                               dataset_index,
                                               :,
                                               ratio_index],
                                color=display_style_["colors"][algorithm_name],
                                linestyle=display_style_["linestyles"][algorithm_name],
                                alpha=display_style_["alpha"])
                            lines.append(line)
                        legends.append(legend.format(algorithm_name))

                axs[dataset_index, ratio_index].legend(lines, legends)
        if save_:
            plt.savefig(save_, bbox_inches='tight')
        else:
            plt.show()
