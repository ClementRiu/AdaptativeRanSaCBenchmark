import glob
import matplotlib.pyplot as plt
import numpy as np

class Metric_Data_Loader:
    def __init__(self,
                 file_glob_format_,
                 algorithms_,
                 num_metrics_,
                 datasets_,
                 sigma_values_,
                 ratio_values_,
                 num_runs_):
        self._file_glob_format = file_glob_format_

        self._algorithms = algorithms_
        self._num_metrics = num_metrics_
        self._datasets = datasets_
        self._sigma_values = sigma_values_
        self._ratio_values = ratio_values_
        self._num_runs = num_runs_

        self._all_values = {
            dataset_name: np.empty((len(self._algorithms),
                                    self._num_metrics,
                                    dataset_size,
                                    len(self._sigma_values),
                                    len(self._ratio_values),
                                    self._num_runs))
            for dataset_name, dataset_size in self._datasets.items()
        }
        for element in self._all_values.values():
            element[:] = np.NaN

        self._metric_names = {
            "p": 0,
            "r": 1,
            "s": 2,
            "t": 3,
        }


    def _format_value_for_filename(self, numerical_value_):
        if numerical_value_ == 0:
            return "0.0"
        else:
            formated_value = "%.1f" % numerical_value_
            if formated_value.startswith("0."):
                return formated_value[1:]
            return formated_value

    def _parse_values(self, values_):
        parsed_values = np.empty((len(self._algorithms), len(values_), self._num_runs))
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
                    metrics[0].append(line)
                elif 'r:' in line:
                    line = file.readline()
                    metrics[1].append(line)
                elif 's:' in line:
                    line = file.readline()
                    metrics[2].append(line)
                elif 't:' in line:
                    line = file.readline()
                    metrics[3].append(line)
                else:
                    line = file.readline()
        return self._parse_values(metrics)

    def load_data(self):
        wrong_dataset = 0
        all_dataset = 0
        for dataset_name, dataset_size in self._datasets.items():
            for dataset_index in range(dataset_size):
                for sigma_index, sigma_value in enumerate(self._sigma_values):
                    for ratio_index, ratio_value in enumerate(self._ratio_values):
                        path = glob.glob(self._file_glob_format.format(dataset_name,
                                                                       dataset_index + 1,
                                                                       self._format_value_for_filename(sigma_value),
                                                                       self._format_value_for_filename(ratio_value)))
                        all_dataset += 1
                        if len(path) != 1:
                            wrong_dataset += 1
                        else:
                            parsed_values = self._open_and_parse(path[0])
                            self._all_values[dataset_name][:, :, dataset_index, sigma_index, ratio_index, :] = parsed_values

        file_glob = self._file_glob_format.format("*", "*", "*", "*")
        all_possible = len(glob.glob(file_glob))
        all_theorical = len(self._sigma_values) * len(self._ratio_values) * np.sum([dataset_size for dataset_size in self._datasets.values()])
        assert(all_possible == all_dataset - wrong_dataset)
        assert(all_theorical == all_dataset)

    def remove_0_and_mean(self):
        self._all_values_without_0 = {
            dataset_name: self._all_values[dataset_name].copy()
            for dataset_name in self._datasets.keys()
        }
        self._mean_values = {
            dataset_name: np.NaN
            for dataset_name in self._datasets.keys()
        }

        for dataset_name in self._datasets.keys():
            self._all_values_without_0[dataset_name][self._all_values_without_0[dataset_name] == 0] = np.NaN
            self._mean_values[dataset_name] = np.nanmean(self._all_values_without_0[dataset_name], axis=5)

    def compute_sigma_ratios(self):
        self._sigma_ratios = {
            dataset_name: self._mean_values[dataset_name][0, self._metric_names["s"], :, :, :].copy()
            for dataset_name in self._datasets.keys()
        }

        for dataset_name in self._datasets.keys():
            for sigma_index, sigma_value in enumerate(self._sigma_values):
                if sigma_value > 0:
                    self._sigma_ratios[dataset_name][:, sigma_index, :] /= sigma_value
        return

    def load_and_clean_data(self):
        self.load_data()
        self.remove_0_and_mean()
        self.compute_sigma_ratios()
        return

class Time_Data_Loader:
    def __init__(self,
                 file_glob_format_,
                 datasize,
                 datanames,
                 labels,
                 numberAlgo):
        self._file_glob_format = file_glob_format_

        self._datasize = datasize
        self._datanames = datanames
        self._labels = labels

        self._numberAlgo = numberAlgo

    def _parse_TimeMetric_line(self, line_):
        line = line_.strip().split()
        line_values = [
            float(line[1]),
            float(line[5]),
            float(line[9][:-1]),
        ]
        line_stds = [
            float(line[3]),
            float(line[7]),
        ]
        return line_values, line_stds

    def _read_TimeMetric_file(self, path_):
        values = []
        with open(path_, 'r') as file:
            line = file.readline()
            while line:
                if "T=" in line:
                    values.append(self._parse_TimeMetric_line(line))
                line=file.readline()
        return values

    def read_values(self):
        self._files_path = []
        for dataname in self._datanames:
            for dataindex in range(1, self._datasize[dataname] + 1):
                self._files_path.append(glob.glob(self._file_glob_format.format(dataname, dataindex))[0])

        self._all_values = {path: self._read_TimeMetric_file(path) for path in self._files_path}
        return

    def parse_TimeMetric_values(self):
        mean_values = []
        std_values = []
        for path in self._files_path:
            value = self._all_values[path]
            for element in value:
                mean_values.append(element[0])
                std_values.append(element[1])
        mean_values = np.array(mean_values)
        std_values = np.array(std_values)

        self._mean = mean_values.reshape(-1, self._numberAlgo, 3).T
        self._std = std_values.reshape(-1, self._numberAlgo, 2).T
        return

    def _create_ratios(self, values_):
        ratios = []
        for indexDim in range(values_.shape[1]):
            ratios.append(values_[:, indexDim, :] / values_[:, -1, :])
        ratios = np.stack(ratios, axis=1)
        return ratios

    def compute_ratios(self):
        self._ratio = self._create_ratios(self._mean)
        return

    def load_data(self):
        self.read_values()
        self.parse_TimeMetric_values()
        self.compute_ratios()
        return
