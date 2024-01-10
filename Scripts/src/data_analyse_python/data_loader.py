import glob
import matplotlib.pyplot as plt
import numpy as np

# from analyse_utils import parse_path
# from magsac_utils import read_magsac, show_dimension, load_magsac_data
from magsac_data import Magsac_Data

# from analyse_utils import load_paths, create_data_container, load_data, summarize_data, generate_image, create_value_container


class All_Data_Loader:
    def __init__(self,
                 editglob_other_paths_,
                 editglob_magsac_path_,
                 folder_names_,
                 data_numbers_,
                 noise_stds_,
                 outlier_ratios_,
                 file_names_=["Weights", "Labels", "Errors", "ErrorsAll", "Inliers", "PosInl"],
                 file_types_=["float", "int", "float", "float", "int", "int"],
                 other_algorithms_=["Ransac", "AC-Ransac", "LRT"],
                 time_algorithms_=["Ransac_other", "Magsac"],
                 magsac_algorithms_=["Magsac-True-Threshold", "Magsac-ORSA-Threshold", "Magsac-ORSA-Precision", "Magsac-ORSA-Recall", "Magsac-Weighted"],
                 model_type_="homographies"):
        # Paths that need formating of :
        # - folder name and dataset number,
        self._editglob_other_paths = editglob_other_paths_
        # - dataset number and file name.
        self._editglob_magsac_path = editglob_magsac_path_
        # Then, the glob.glob method will find all files and make them correspond to the noise_stds_ and outlier_ratios_.

        self._number_of_datasets = len(self._editglob_other_paths)

        # List of elements to iterate other :
        self._folder_names = folder_names_  # Folder names where data are stored.
        self._data_numbers = data_numbers_  # All possible dataset numbers.
        self._noise_stds = noise_stds_  # All possible values of noise std.
        self._outlier_ratios = outlier_ratios_  # All possible values of outlier ratio.

        self._file_names = file_names_  # All Magsac file names.
        self._file_types = file_types_  # All Magsac file type for reading.

        self._file_index = {file_name: index for index, file_name in enumerate(self._file_names)}  # Correspondance dict between file names and index.

        self._other_algorithms = other_algorithms_  # Name of all algorithms but Magsac.
        self._time_algorithms = time_algorithms_  # Name of additional algorithms.
        self._magsac_algorithms = magsac_algorithms_  # Name of Magsac algorithm.

        self._model_type = model_type_

    def _magsac_container(self):
        '''
        Create a data container for magsac paths, data, ...
        Datas are stored in a list in a dict indexed by
        dataset number, noise std and outlier ratio.
        '''
        return {data_number: {noise_std: {outlier_ratio: []
                                          for outlier_ratio in self._outlier_ratios}
                              for noise_std in self._noise_stds}
                for data_number in self._data_numbers}

    def _data_container(self):
        '''
        Create a data container for other data, ...
        Datas are stored in a list in a dict indexed by
        folder name, dataset number, noise std and outlier ratio.
        '''
        return {folder_name: {data_number: {noise_std: {outlier_ratio: []
                                                        for outlier_ratio in self._outlier_ratios}
                                            for noise_std in self._noise_stds}
                              for data_number in self._data_numbers}
                for folder_name in self._folder_names}

    def _value_container(self, algorithms_):
        '''
        Create a data container for mean and std values.
        Datas are stored in a matrix of size (#noise stds x #outlier ratios) in a dict indexed by
        folder name, dataset number, algorithm.
        '''
        num_noise_stds = len(self._noise_stds)
        num_outlier_ratios = len(self._outlier_ratios)
        return {folder_name: {data_number: {algorithm: np.empty([num_noise_stds, num_outlier_ratios])
                                            for algorithm in algorithms_}
                              for data_number in self._data_numbers}
                for folder_name in self._folder_names}

    def _parse_path(self, path_):
        '''
        Parse a path to extract the value of the noise std and the value of the outlier ratio.
        '''
        parsed_path_std = path_.split("_std")[-1]
        parsed_path_std = parsed_path_std.split("_")[0]

        parsed_path_ratio = path_.split("_ratio")[-1]
        parsed_path_ratio = parsed_path_ratio.split("_")[0]

        return float(parsed_path_std), float(parsed_path_ratio)

    def _create_other_path(self):
        '''
        Load all possible path corresponding to _editglob_other_path.
        '''
        self._other_paths = [{folder_name: {data_number: {noise_std: []
                                                          for noise_std in self._noise_stds}
                                            for data_number in self._data_numbers}
                              for folder_name in self._folder_names}
                             for _ in range(self._number_of_datasets)]

        for index_dataset, editglob_other_path in enumerate(self._editglob_other_paths):
            for folder_name in self._folder_names:
                for data_number in self._data_numbers:
                    other_path_globededited = sorted(glob.glob(editglob_other_path.format(folder_name, data_number)))
                    for path in other_path_globededited:
                        noise_std, _ = self._parse_path(path)
                        self._other_paths[index_dataset][folder_name][data_number][noise_std].append(path)

        return None

    def _create_magsac_path(self):
        '''
        Load all possible path corresponding to _editglob_magsac_path.
        '''
        self._magsac_paths = self._magsac_container()

        for data_number in self._data_numbers:
            for file_name in self._file_names:
                magsac_path_globededited = sorted(glob.glob(self._editglob_magsac_path.format(data_number, file_name)))
                for path in magsac_path_globededited:
                    true_noise_std, true_outlier_ratio = self._parse_path(path)
                    self._magsac_paths[data_number][true_noise_std][true_outlier_ratio].append(path)
        return None

    def create_paths(self):
        '''
        Load both other and Magsac paths.
        '''
        self._create_other_path()
        self._create_magsac_path()
        return None

    def _parse_values(self, values_):
        '''
        Parse a list of string representing a list of vector.
        '''
        return [[float(value) for value in string.rstrip().split() if value != "-nan" and value !="nan"] for string in values_]

    def _load_other_data(self, index_):
        '''
        Load data of other algorithms.
        '''
        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for noise_std in self._noise_stds:
                    for path in self._other_paths[index_][folder_name][data_number][noise_std]:
                        p_values = []
                        r_values = []
                        s_values = ['-1 ']
                        t_values = []
                        _, outlier_ratio = self._parse_path(path)
                        with open(path, 'r') as file:
                            line = file.readline()
                            while line:
                                if 'p:' in line:
                                    line = file.readline()
                                    p_values.append(line)
                                elif 'r:' in line:
                                    line = file.readline()
                                    r_values.append(line)
                                elif 'sigma:' in line:
                                    line = file.readline()
                                    s_values.append(line)
                                elif 't:' in line:
                                    line = file.readline()
                                    t_values.append(line)
                                else:
                                    line = file.readline()
                        self._other_p_datas[index_][folder_name][data_number][noise_std][outlier_ratio] = self._parse_values(p_values)
                        self._other_r_datas[index_][folder_name][data_number][noise_std][outlier_ratio] = self._parse_values(r_values)
                        self._other_s_datas[index_][folder_name][data_number][noise_std][outlier_ratio] = self._parse_values(s_values)
                        self._other_t_datas[index_][folder_name][data_number][noise_std][outlier_ratio] = self._parse_values(t_values)
        return None

    def _read_magsac(self, path_, type_used_=float):
        data = []
        with open(path_, "r") as file:
            for line in file:
                element = np.array([value for value in line.split(", ")[:-1]]).astype(type_used_)
                data.append(element)
        return data

    def _load_magsac_data(self):
        magsac_data = self._magsac_container()

        for data_number in self._data_numbers:
            for noise_std in self._noise_stds:
                for outlier_ratio in self._outlier_ratios:
                    # Add error catching for inexistant file.
                    all_path_names = self._magsac_paths[data_number][noise_std][outlier_ratio]
                    if len(all_path_names) == 6:
                        magsac_data[data_number][noise_std][outlier_ratio] = \
                        [self._read_magsac(all_path_names[i], file_type) for i, file_type in enumerate(self._file_types)]
                    else:
                        magsac_data[data_number][noise_std][outlier_ratio] = None

        for data_number in self._data_numbers:
            for noise_std in self._noise_stds:
                for outlier_ratio in self._outlier_ratios:
                    data_for_process = magsac_data[data_number][noise_std][outlier_ratio]
                    if data_for_process != None:
                        data_to_analyse = [Magsac_Data(
                            data_for_process[self._file_index["Labels"]][i // 5],
                            data_for_process[self._file_index["ErrorsAll"]][i],
                            # data_for_process[self._file_index["FullErrors"]][i],
                            data_for_process[self._file_index["PosInl"]][i],
                            # data_for_process[self._file_index["Inliers"]][i],
                            data_for_process[self._file_index["Weights"]][i])
                                for i in range(len(data_for_process[0]))]
                    else:
                        data_to_analyse = [Magsac_Data(None, None, None, None)]
                    self._magsac_data_for_process[data_number][noise_std][outlier_ratio] = data_to_analyse

        return None

    def load_data(self):
        self._other_p_datas = [self._data_container() for _ in range(self._number_of_datasets)]
        self._other_r_datas = [self._data_container() for _ in range(self._number_of_datasets)]
        self._other_s_datas = [self._data_container() for _ in range(self._number_of_datasets)]
        self._other_t_datas = [self._data_container() for _ in range(self._number_of_datasets)]

        for index in range(self._number_of_datasets):
            self._load_other_data(index)

        self._magsac_data_for_process = self._magsac_container()

        self._load_magsac_data()

        return None

    def _summarize_data(self, algorithms_, data_, threshold_=0):
        values_mean = self._value_container(algorithms_)
        values_std = self._value_container(algorithms_)
        values_kept = self._value_container(algorithms_)

        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for algo_index, algorithm in enumerate(algorithms_):
                    for noise_std_index, noise_std in enumerate(self._noise_stds):
                        for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                            try:
                                values = np.array(data_[folder_name][data_number][noise_std][outlier_ratio][algo_index])
                                values_to_keep = values > threshold_
                                values_mean[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = np.nan_to_num(np.mean(values[values_to_keep]), nan=-1)
                                values_std[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = np.nan_to_num(np.std(values[values_to_keep]), nan=0)
                                values_kept[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = (values_to_keep).sum()
                            except IndexError:
                                kept = 0
                                values_mean[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = -1
                                values_std[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = 0
                                values_kept[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = 0
        return values_mean, values_std, values_kept

    def _true_threshold(self, folder_name_, data_number_, noise_std_, outlier_ratio_, noise_std_index_, outlier_ratio_index_):
        return noise_std_

    def _orsa_threshold(self, folder_name_, data_number_, noise_std_, outlier_ratio_, noise_std_index_, outlier_ratio_index_):
        return self._other_sigma_mean[folder_name_][data_number_]["AC-Ransac"][noise_std_index_, outlier_ratio_index_]

    def _get_metrics_from_fixed_threshold(self, threshold_function_):
        magsac_p_datas = self._data_container()
        magsac_r_datas = self._data_container()
        magsac_s_datas = self._data_container()

        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for noise_std_index, noise_std in enumerate(self._noise_stds):
                    for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                        cutoff_threshold = threshold_function_(folder_name, data_number, noise_std, outlier_ratio, noise_std_index, outlier_ratio_index)
                        magsac_s_datas[folder_name][data_number][noise_std][outlier_ratio] = [[cutoff_threshold]]

                        values_to_transform = self._magsac_data_for_process[data_number][noise_std][outlier_ratio]

                        magsac_computed_p = [value_to_transform.get_precision_from_threshold(cutoff_threshold)
                                             for value_to_transform in values_to_transform]
                        magsac_p_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_p]

                        magsac_computed_r = [value_to_transform.get_recall_from_threshold(cutoff_threshold)
                                             for value_to_transform in values_to_transform]
                        magsac_r_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_r]
        return magsac_p_datas, magsac_r_datas, magsac_s_datas


    def _get_metrics_from_orsa_var(self, orsa_var_mean_, var_name_):
        magsac_p_datas = self._data_container()
        magsac_r_datas = self._data_container()
        magsac_s_datas = self._data_container()

        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for noise_std_index, noise_std in enumerate(self._noise_stds):
                    for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                        orsa_var = orsa_var_mean_[folder_name][data_number]["AC-Ransac"][noise_std_index, outlier_ratio_index]
                        values_to_transform = self._magsac_data_for_process[data_number][noise_std][outlier_ratio]

                        magsac_computed_p = []
                        magsac_computed_r = []
                        magsac_computed_s = []

                        for value_to_transform in values_to_transform:
                            try:
                                cutoff_threshold, precision, recall = value_to_transform.get_metrics_from_name(orsa_var, var_name_)
                            except IndexError:
                                cutoff_threshold, precision, recall = -1, -1, -1

                            magsac_computed_s.append(cutoff_threshold)
                            magsac_computed_p.append(precision)
                            magsac_computed_r.append(recall)

                        magsac_p_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_p]
                        magsac_r_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_r]
                        magsac_s_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_s]

        return magsac_p_datas, magsac_r_datas, magsac_s_datas

    def _get_weighted_metrics(self):
        magsac_p_datas = self._data_container()
        magsac_r_datas = self._data_container()
        magsac_s_datas = self._data_container()

        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for noise_std_index, noise_std in enumerate(self._noise_stds):
                    for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                        magsac_s_datas[folder_name][data_number][noise_std][outlier_ratio] = [[-1]]

                        values_to_transform = self._magsac_data_for_process[data_number][noise_std][outlier_ratio]
                        magsac_computed_prs = [value_to_transform.get_weighted_metrics()
                                             for value_to_transform in values_to_transform]

                        magsac_computed_p = [magsac_computed_pr[0] for magsac_computed_pr in magsac_computed_prs]
                        magsac_computed_r = [magsac_computed_pr[1] for magsac_computed_pr in magsac_computed_prs]

                        magsac_p_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_p]
                        magsac_r_datas[folder_name][data_number][noise_std][outlier_ratio] = [magsac_computed_r]
        return magsac_p_datas, magsac_r_datas, magsac_s_datas

    def _get_summary_other(self, threshold_=0):

        self._other_precision_mean, self._other_precision_std, self._other_precision_kept = self._summarize_data(self._other_algorithms, self._other_p_datas[0], threshold_)
        self._other_recall_mean, self._other_recall_std, self._other_recall_kept = self._summarize_data(self._other_algorithms, self._other_r_datas[0], threshold_)
        self._other_sigma_mean, self._other_sigma_std, _ = self._summarize_data(self._other_algorithms, self._other_s_datas[0])
        self._other_time_mean, self._other_time_std, _ = self._summarize_data(self._other_algorithms, self._other_t_datas[0])

        # self._ransac_precision_mean, self._ransac_precision_std, self._ransac_precision_kept = self._summarize_data([self._time_algorithms[0]], self._other_p_datas[1], threshold_)
        # self._ransac_recall_mean, self._ransac_recall_std, self._ransac_recall_kept = self._summarize_data([self._time_algorithms[0]], self._other_r_datas[1], threshold_)
        # self._ransac_sigma_mean, self._ransac_sigma_std, _ = self._summarize_data([self._time_algorithms[0]], self._other_s_datas[1])
        # self._ransac_time_mean, self._ransac_time_std, _ = self._summarize_data([self._time_algorithms[0]], self._other_t_datas[1])

        self._magsac_time_mean, self._magsac_time_std, _ = self._summarize_data([self._time_algorithms[0]], self._other_t_datas[1])

        return None

    def get_metrics(self, threshold_=0):
        self._get_summary_other(threshold_)

        self._magsac_p_datas_tt, self._magsac_r_datas_tt, self._magsac_s_datas_tt = self._get_metrics_from_fixed_threshold(self._true_threshold)
        self._magsac_p_datas_ot, self._magsac_r_datas_ot, self._magsac_s_datas_ot = self._get_metrics_from_fixed_threshold(self._orsa_threshold)
        self._magsac_p_datas_op, self._magsac_r_datas_op, self._magsac_s_datas_op = self._get_metrics_from_orsa_var(self._other_precision_mean, 'precision')
        self._magsac_p_datas_or, self._magsac_r_datas_or, self._magsac_s_datas_or = self._get_metrics_from_orsa_var(self._other_recall_mean, 'recall')
        self._magsac_p_datas_w, self._magsac_r_datas_w, self._magsac_s_datas_w = self._get_weighted_metrics()

        return None


    def create_summary(self, threshold_=0):
        self._magsac_precision_mean_tt, self._magsac_precision_std_tt, self._magsac_precision_kept_tt = self._summarize_data([self._magsac_algorithms[0]], self._magsac_p_datas_tt, threshold_)
        self._magsac_recall_mean_tt, self._magsac_recall_std_tt, self._magsac_recall_kept_tt = self._summarize_data([self._magsac_algorithms[0]], self._magsac_r_datas_tt, threshold_)
        self._magsac_sigma_mean_tt, self._magsac_sigma_std_tt, _ = self._summarize_data([self._magsac_algorithms[0]], self._magsac_s_datas_tt)

        self._magsac_precision_mean_ot, self._magsac_precision_std_ot, self._magsac_precision_kept_ot = self._summarize_data([self._magsac_algorithms[1]], self._magsac_p_datas_ot, threshold_)
        self._magsac_recall_mean_ot, self._magsac_recall_std_ot, self._magsac_recall_kept_ot = self._summarize_data([self._magsac_algorithms[1]], self._magsac_r_datas_ot, threshold_)
        self._magsac_sigma_mean_ot, self._magsac_sigma_std_ot, _ = self._summarize_data([self._magsac_algorithms[1]], self._magsac_s_datas_ot)

        self._magsac_precision_mean_op, self._magsac_precision_std_op, self._magsac_precision_kept_op = self._summarize_data([self._magsac_algorithms[2]], self._magsac_p_datas_op, threshold_)
        self._magsac_recall_mean_op, self._magsac_recall_std_op, self._magsac_recall_kept_op = self._summarize_data([self._magsac_algorithms[2]], self._magsac_r_datas_op, threshold_)
        self._magsac_sigma_mean_op, self._magsac_sigma_std_op, _ = self._summarize_data([self._magsac_algorithms[2]], self._magsac_s_datas_op)

        self._magsac_precision_mean_or, self._magsac_precision_std_or, self._magsac_precision_kept_or = self._summarize_data([self._magsac_algorithms[3]], self._magsac_p_datas_or, threshold_)
        self._magsac_recall_mean_or, self._magsac_recall_std_or, self._magsac_recall_kept_or = self._summarize_data([self._magsac_algorithms[3]], self._magsac_r_datas_or, threshold_)
        self._magsac_sigma_mean_or, self._magsac_sigma_std_or, _ = self._summarize_data([self._magsac_algorithms[3]], self._magsac_s_datas_or)

        self._magsac_precision_mean_w, self._magsac_precision_std_w, self._magsac_precision_kept_w = self._summarize_data([self._magsac_algorithms[4]], self._magsac_p_datas_w, threshold_)
        self._magsac_recall_mean_w, self._magsac_recall_std_w, self._magsac_recall_kept_w = self._summarize_data([self._magsac_algorithms[4]], self._magsac_r_datas_w, threshold_)
        self._magsac_sigma_mean_w, self._magsac_sigma_std_w, _ = self._summarize_data([self._magsac_algorithms[4]], self._magsac_s_datas_w)

        return None

    def _concatenate_values(self, algorithms_, value1_, value2_, separation=1):

        values_concat = self._value_container(algorithms_)
        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for algo_index, algorithm in enumerate(algorithms_):
                    # for noise_std_index, noise_std in enumerate(self._noise_stds):
                    #     for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                    # print(algo_index, algorithm, separation)
                    if algo_index < separation:
                        exctracted_data = value1_[folder_name][data_number][algorithm]
                    else:
                        exctracted_data = value2_[folder_name][data_number][algorithm]
                    values_concat[folder_name][data_number][algorithm] = exctracted_data
        return values_concat

    def concatenate_all_data(self):
        self._all_algorithms = self._magsac_algorithms[0:2]
        self._precision_mean = self._concatenate_values(self._all_algorithms, self._magsac_precision_mean_tt, self._magsac_precision_mean_ot)
        self._precision_std = self._concatenate_values(self._all_algorithms, self._magsac_precision_std_tt, self._magsac_precision_std_ot)
        self._recall_mean = self._concatenate_values(self._all_algorithms, self._magsac_recall_mean_tt, self._magsac_recall_mean_ot)
        self._recall_std = self._concatenate_values(self._all_algorithms, self._magsac_recall_std_tt, self._magsac_recall_std_ot)
        self._sigma_mean = self._concatenate_values(self._all_algorithms, self._magsac_sigma_mean_tt, self._magsac_sigma_mean_ot)
        self._sigma_std = self._concatenate_values(self._all_algorithms, self._magsac_sigma_std_tt, self._magsac_sigma_std_ot)

        self._precision_kept = self._concatenate_values(self._all_algorithms, self._magsac_precision_kept_tt, self._magsac_precision_kept_ot)
        self._recall_kept = self._concatenate_values(self._all_algorithms, self._magsac_recall_kept_tt, self._magsac_recall_kept_ot)

        self._all_algorithms = self._magsac_algorithms[0:3]
        self._precision_mean = self._concatenate_values(self._all_algorithms, self._precision_mean, self._magsac_precision_mean_op, 2)
        self._precision_std = self._concatenate_values(self._all_algorithms, self._precision_std, self._magsac_precision_std_op, 2)
        self._recall_mean = self._concatenate_values(self._all_algorithms, self._recall_mean, self._magsac_recall_mean_op, 2)
        self._recall_std = self._concatenate_values(self._all_algorithms, self._recall_std, self._magsac_recall_std_op, 2)
        self._sigma_mean = self._concatenate_values(self._all_algorithms, self._sigma_mean, self._magsac_sigma_mean_op, 2)
        self._sigma_std = self._concatenate_values(self._all_algorithms, self._sigma_std, self._magsac_sigma_std_op, 2)

        self._precision_kept = self._concatenate_values(self._all_algorithms, self._precision_kept, self._magsac_precision_kept_op, 2)
        self._recall_kept = self._concatenate_values(self._all_algorithms, self._recall_kept, self._magsac_recall_kept_op, 2)

        self._all_algorithms = self._magsac_algorithms[0:4]
        self._precision_mean = self._concatenate_values(self._all_algorithms, self._precision_mean, self._magsac_precision_mean_or, 3)
        self._precision_std = self._concatenate_values(self._all_algorithms, self._precision_std, self._magsac_precision_std_or, 3)
        self._recall_mean = self._concatenate_values(self._all_algorithms, self._recall_mean, self._magsac_recall_mean_or, 3)
        self._recall_std = self._concatenate_values(self._all_algorithms, self._recall_std, self._magsac_recall_std_or, 3)
        self._sigma_mean = self._concatenate_values(self._all_algorithms, self._sigma_mean, self._magsac_sigma_mean_or, 3)
        self._sigma_std = self._concatenate_values(self._all_algorithms, self._sigma_std, self._magsac_sigma_std_or, 3)

        self._precision_kept = self._concatenate_values(self._all_algorithms, self._precision_kept, self._magsac_precision_kept_or, 3)
        self._recall_kept = self._concatenate_values(self._all_algorithms, self._recall_kept, self._magsac_recall_kept_or, 3)

        self._all_algorithms = self._magsac_algorithms
        self._precision_mean = self._concatenate_values(self._all_algorithms, self._precision_mean, self._magsac_precision_mean_w, 4)
        self._precision_std = self._concatenate_values(self._all_algorithms, self._precision_std, self._magsac_precision_std_w, 4)
        self._recall_mean = self._concatenate_values(self._all_algorithms, self._recall_mean, self._magsac_recall_mean_w, 4)
        self._recall_std = self._concatenate_values(self._all_algorithms, self._recall_std, self._magsac_recall_std_w, 4)
        self._sigma_mean = self._concatenate_values(self._all_algorithms, self._sigma_mean, self._magsac_sigma_mean_w, 4)
        self._sigma_std = self._concatenate_values(self._all_algorithms, self._sigma_std, self._magsac_sigma_std_w, 4)

        self._precision_kept = self._concatenate_values(self._all_algorithms, self._precision_kept, self._magsac_precision_kept_w, 4)
        self._recall_kept = self._concatenate_values(self._all_algorithms, self._recall_kept, self._magsac_recall_kept_w, 4)

        self._all_algorithms = self._magsac_algorithms + self._other_algorithms
        self._precision_mean = self._concatenate_values(self._all_algorithms, self._precision_mean, self._other_precision_mean, 5)
        self._precision_std = self._concatenate_values(self._all_algorithms, self._precision_std, self._other_precision_std, 5)
        self._recall_mean = self._concatenate_values(self._all_algorithms, self._recall_mean, self._other_recall_mean, 5)
        self._recall_std = self._concatenate_values(self._all_algorithms, self._recall_std, self._other_recall_std, 5)
        self._sigma_mean = self._concatenate_values(self._all_algorithms, self._sigma_mean, self._other_sigma_mean, 5)
        self._sigma_std = self._concatenate_values(self._all_algorithms, self._sigma_std, self._other_sigma_std, 5)

        self._precision_kept = self._concatenate_values(self._all_algorithms, self._precision_kept, self._other_precision_kept, 5)
        self._recall_kept = self._concatenate_values(self._all_algorithms, self._recall_kept, self._other_recall_kept, 5)

        # self._all_algorithms = self._all_algorithms + [self._time_algorithms[0]]
        # self._precision_mean = self._concatenate_values(self._all_algorithms, self._precision_mean, self._ransac_precision_mean, 8)
        # self._precision_std = self._concatenate_values(self._all_algorithms, self._precision_std, self._ransac_precision_std, 8)
        # self._recall_mean = self._concatenate_values(self._all_algorithms, self._recall_mean, self._ransac_recall_mean, 8)
        # self._recall_std = self._concatenate_values(self._all_algorithms, self._recall_std, self._ransac_recall_std, 8)
        # self._sigma_mean = self._concatenate_values(self._all_algorithms, self._sigma_mean, self._ransac_sigma_mean, 8)
        # self._sigma_std = self._concatenate_values(self._all_algorithms, self._sigma_std, self._ransac_sigma_std, 8)
        #
        # self._precision_kept = self._concatenate_values(self._all_algorithms, self._precision_kept, self._ransac_precision_kept, 8)
        # self._recall_kept = self._concatenate_values(self._all_algorithms, self._recall_kept, self._ransac_recall_kept, 8)

        # self._all_time_algorithms = self._time_algorithms
        # self._time_mean = self._concatenate_values(self._all_time_algorithms, self._magsac_time_mean)
        # self._time_std = self._concatenate_values(self._all_time_algorithms, self._magsac_time_std)
        self._all_time_algorithms = self._time_algorithms + self._other_algorithms
        self._time_mean = self._concatenate_values(self._all_time_algorithms, self._magsac_time_mean, self._other_time_mean)
        self._time_std = self._concatenate_values(self._all_time_algorithms, self._magsac_time_std, self._other_time_std)

        return None

    def _harmonic_mean(self, value1_, value2_):
        if value1_ > 0 and value2_ > 0:
            return (2 * value1_ * value2_) / (value1_ + value2_)
        else:
            return 0

    def _compute_f1score(self, precision_, recall_):
        f1score = self._data_container()
        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for noise_std_index, noise_std in enumerate(self._noise_stds):
                    for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                        precision = precision_[folder_name][data_number][noise_std][outlier_ratio]
                        recall = recall_[folder_name][data_number][noise_std][outlier_ratio]
                        f1score[folder_name][data_number][noise_std][outlier_ratio] = [[self._harmonic_mean(precision[i][j], recall[i][j])
                                                                                        for j in range(len(precision[i]))]
                                                                                       for i in range(len(precision))]
        return f1score

    def generate_f1score(self, threshold_=0):
        self._other_f1score = [self._compute_f1score(self._other_p_datas[i], self._other_r_datas[i]) for i in range(2)]
        self._magsac_f1score_tt = self._compute_f1score(self._magsac_p_datas_tt, self._magsac_r_datas_tt)
        self._magsac_f1score_ot = self._compute_f1score(self._magsac_p_datas_ot, self._magsac_r_datas_ot)
        self._magsac_f1score_op = self._compute_f1score(self._magsac_p_datas_op, self._magsac_r_datas_op)
        self._magsac_f1score_or = self._compute_f1score(self._magsac_p_datas_or, self._magsac_r_datas_or)
        self._magsac_f1score_w = self._compute_f1score(self._magsac_p_datas_w, self._magsac_r_datas_w)

        self._other_f1score_mean, self._other_f1score_std, self._other_f1score_kept = self._summarize_data(self._other_algorithms, self._other_f1score[0], threshold_)
        self._ransac_f1score_mean, self._ransac_f1score_std, self._ransac_f1score_kept = self._summarize_data([self._time_algorithms[0]], self._other_f1score[1], threshold_)
        self._magsac_f1score_mean_tt, self._magsac_f1score_std_tt, self._magsac_f1score_kept_tt = self._summarize_data([self._magsac_algorithms[0]], self._magsac_f1score_tt, threshold_)
        self._magsac_f1score_mean_ot, self._magsac_f1score_std_ot, self._magsac_f1score_kept_ot = self._summarize_data([self._magsac_algorithms[1]], self._magsac_f1score_ot, threshold_)
        self._magsac_f1score_mean_op, self._magsac_f1score_std_op, self._magsac_f1score_kept_op = self._summarize_data([self._magsac_algorithms[2]], self._magsac_f1score_op, threshold_)
        self._magsac_f1score_mean_or, self._magsac_f1score_std_or, self._magsac_f1score_kept_or = self._summarize_data([self._magsac_algorithms[3]], self._magsac_f1score_or, threshold_)
        self._magsac_f1score_mean_w, self._magsac_f1score_std_w, self._magsac_f1score_kept_w = self._summarize_data([self._magsac_algorithms[4]], self._magsac_f1score_w, threshold_)

        self._all_algorithms_f1 = self._magsac_algorithms[0:2]
        self._f1score_mean = self._concatenate_values(self._all_algorithms_f1, self._magsac_f1score_mean_tt, self._magsac_f1score_mean_ot)
        self._f1score_std = self._concatenate_values(self._all_algorithms_f1, self._magsac_f1score_std_tt, self._magsac_f1score_std_ot)
        self._f1score_kept = self._concatenate_values(self._all_algorithms_f1, self._magsac_f1score_kept_tt, self._magsac_f1score_kept_ot)

        self._all_algorithms_f1 = self._magsac_algorithms[0:3]
        self._f1score_mean = self._concatenate_values(self._all_algorithms_f1, self._f1score_mean, self._magsac_f1score_mean_op, 2)
        self._f1score_std = self._concatenate_values(self._all_algorithms_f1, self._f1score_std, self._magsac_f1score_std_op, 2)
        self._f1score_kept = self._concatenate_values(self._all_algorithms_f1, self._f1score_kept, self._magsac_f1score_kept_op, 2)

        self._all_algorithms_f1 = self._magsac_algorithms[0:4]
        self._f1score_mean = self._concatenate_values(self._all_algorithms_f1, self._f1score_mean, self._magsac_f1score_mean_or, 3)
        self._f1score_std = self._concatenate_values(self._all_algorithms_f1, self._f1score_std, self._magsac_f1score_std_or, 3)
        self._f1score_kept = self._concatenate_values(self._all_algorithms_f1, self._f1score_kept, self._magsac_f1score_kept_or, 3)

        self._all_algorithms_f1 = self._magsac_algorithms
        self._f1score_mean = self._concatenate_values(self._all_algorithms_f1, self._f1score_mean, self._magsac_f1score_mean_w, 4)
        self._f1score_std = self._concatenate_values(self._all_algorithms_f1, self._f1score_std, self._magsac_f1score_std_w, 4)
        self._f1score_kept = self._concatenate_values(self._all_algorithms_f1, self._f1score_kept, self._magsac_f1score_kept_w, 4)

        self._all_algorithms_f1 = self._magsac_algorithms + self._other_algorithms
        self._f1score_mean = self._concatenate_values(self._all_algorithms_f1, self._f1score_mean, self._other_f1score_mean, 5)
        self._f1score_std = self._concatenate_values(self._all_algorithms_f1, self._f1score_std, self._other_f1score_std, 5)
        self._f1score_kept = self._concatenate_values(self._all_algorithms_f1, self._f1score_kept, self._other_f1score_kept, 5)

        # self._all_algorithms_f1 = self._all_algorithms_f1 + [self._time_algorithms[0]]
        # self._f1score_mean = self._concatenate_values(self._all_algorithms_f1, self._f1score_mean, self._ransac_f1score_mean, 9)
        # self._f1score_std = self._concatenate_values(self._all_algorithms_f1, self._f1score_std, self._ransac_f1score_std, 9)
        # self._f1score_kept = self._concatenate_values(self._all_algorithms_f1, self._f1score_kept, self._ransac_f1score_kept, 9)

    def _generate_image(self, data_number_, algorithms_, value_name_, values_mean_, values_std_, std_=True, save_=False, threshold_=None):
        '''
        Just for 1 dataset.
        '''
        alpha = 0.3
        colors = {self._other_algorithms[0]: 'b',
                  self._other_algorithms[1]: 'b',
                  self._other_algorithms[2]: 'r',
                  self._other_algorithms[3]: 'g',
                  self._time_algorithms[0]: 'c',
                  'Magsac-True-Threshold': 'c',
                  'Magsac-AC-Ransac-Threshold': 'c',
                  'Magsac-P': 'm',
                  'Magsac-R': 'm',
                  'Magsac-W': 'c',
                  'MUSE': 'y',
                  'Fast-AC-Ransac': 'r',
                  }
        linestyles = {self._other_algorithms[0]: '-',
                      self._other_algorithms[1]: '--',
                      self._other_algorithms[2]: '-',
                      self._other_algorithms[3]: '-',
                      self._time_algorithms[0]: '-',
                      'Magsac-True-Threshold': '-',
                      'Magsac-AC-Ransac-Threshold': '--',
                      'Magsac-P': ':',
                      'Magsac-R': '-.',
                      'Magsac-W': '-.',
                      'MUSE': '-',
                      'Fast-AC-Ransac': '--',
                      }
        legend = "{} algorithm for {} noise."
        xlabel_n = "Noise std value (in pixel)"
        xlabel_o = "Outlier ratio value"
        generic_title = "{} for outlier ratio of {} on dataset {}."
        generic_suptitle = "{} of " + "{} for only inlier, with noise, artificial {}.".format(" ".join(algorithms_), self._model_type)

        fig, axs = plt.subplots(1, len(self._outlier_ratios), figsize=(10 * len(self._outlier_ratios), 10 * 1), sharey=True)
        fig.suptitle(generic_suptitle.format(value_name_))
        for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
            axs[outlier_ratio_index].set_title(generic_title.format(value_name_, outlier_ratio, data_number_))
            axs[outlier_ratio_index].set_xlabel(xlabel_n)
            axs[outlier_ratio_index].set_ylabel(value_name_)
            if threshold_ and threshold_ >= 1:
                axs[outlier_ratio_index].set_ylim(-1.1, threshold_)
            if threshold_ and threshold_ < 1:
                    axs[outlier_ratio_index].set_ylim(threshold_, 1.05)
            lines = []
            legends = []
            for algo_index, algorithm in enumerate(algorithms_):
                for folder_index, folder_name in enumerate(self._folder_names):
                    if 'Magsac' in algorithm:
                        # print(values_mean_[folder_name][data_number_][algorithm][:, outlier_ratio_index])
                        if std_:
                            line = axs[outlier_ratio_index].errorbar(
                                    self._noise_stds,
                                    values_mean_[folder_name][data_number_][algorithm][:, outlier_ratio_index],
                                    yerr=values_std_[folder_name][data_number_][algorithm][:, outlier_ratio_index],
                                    color=colors[algorithm],
                                    linestyle=linestyles[algorithm],
                                    alpha=alpha
                                    )
                        else:
                            line,  = axs[outlier_ratio_index].plot(
                                    self._noise_stds,
                                    values_mean_[folder_name][data_number_][algorithm][:, outlier_ratio_index],
                                    color=colors[algorithm],
                                    linestyle=linestyles[algorithm],
                                    alpha=alpha
                                    )
                        lines.append(line)
                        legends.append(legend.format(algorithm, folder_name))
                    else:
                        if std_:
                            line  = axs[outlier_ratio_index].errorbar(
                                    self._noise_stds,
                                    values_mean_[folder_name][data_number_][algorithm][:, outlier_ratio_index],
                                    yerr=values_std_[folder_name][data_number_][algorithm][:, outlier_ratio_index],
                                    color=colors[algorithm],
                                    linestyle=linestyles[algorithm],
                                    alpha=alpha
                                    )
                        else:
                            line, =  axs[outlier_ratio_index].plot(
                                    self._noise_stds,
                                    values_mean_[folder_name][data_number_][algorithm][:, outlier_ratio_index],
                                    color=colors[algorithm],
                                    linestyle=linestyles[algorithm],
                                    alpha=alpha
                                    )
                        lines.append(line)
                        legends.append(legend.format(algorithm, folder_name))
            axs[outlier_ratio_index].legend(lines, legends)
        if save_:
            plt.savefig(save_)
        else:
            plt.show()

        return None

    def generate_single_image(self, folder_name_, data_number_, algorithms_, remove_, outlier_ratio_, outlier_ratio_index_, value_name_, values_mean_, values_std_, std_=True, save_=False, threshold_=None):
        '''
        Just for 1 dataset.
        '''
        local_algorithms = [algo for algo in algorithms_ if algo not in remove_]

        alpha = 1.0
        colors = {self._other_algorithms[0]: 'b',
                  self._other_algorithms[1]: 'b',
                  self._other_algorithms[2]: 'r',
                  self._other_algorithms[3]: 'g',
                  self._time_algorithms[0]: 'c',
                  'Magsac-True-Threshold': 'c',
                  'Magsac-AC-Ransac-Threshold': 'c',
                  'Magsac-P': 'm',
                  'Magsac-R': 'm',
                  'Magsac-W': 'c',
                  'MUSE': 'y',
                  'Fast-AC-Ransac': 'r',
                  }
        linestyles = {self._other_algorithms[0]: '-',
                      self._other_algorithms[1]: '--',
                      self._other_algorithms[2]: '-',
                      self._other_algorithms[3]: '-',
                      self._time_algorithms[0]: '-',
                      'Magsac-True-Threshold': '-',
                      'Magsac-AC-Ransac-Threshold': '--',
                      'Magsac-P': ':',
                      'Magsac-R': '-.',
                      'Magsac-W': '-.',
                      'MUSE': '-',
                      'Fast-AC-Ransac': '--',
                      }
        legend = "{}."
        xlabel_n = "Noise std value (in pixel)"
        xlabel_o = "Outlier ratio value"
        generic_title = "{} for outlier ratio of {} on dataset {}."
        generic_suptitle = "{} " + "for outlier ratio of {} on artificial {} - dataset {}.".format(outlier_ratio_, self._model_type, data_number_)

        fig, axs = plt.subplots(1, 1, figsize=(8, 4), sharey=True)
        fig.suptitle(generic_suptitle.format(value_name_))
        # axs.set_title(generic_title.format(value_name_, outlier_ratio_, data_number_))
        axs.set_xlabel(xlabel_n)
        axs.set_ylabel(value_name_)
        if threshold_ and threshold_ >= 1:
            axs.set_ylim(-1.1, threshold_)
        if threshold_ and threshold_ < 1:
                axs.set_ylim(threshold_, 1.05)
        lines = []
        legends = []
        for algo_index, algorithm in enumerate(local_algorithms):
            if 'Magsac' in algorithm:
                # print(values_mean_[folder_name][data_number_][algorithm][:, outlier_ratio_index])
                if std_:
                    line = axs.errorbar(
                            self._noise_stds,
                            values_mean_[folder_name_][data_number_][algorithm][:, outlier_ratio_index_],
                            yerr=values_std_[folder_name][data_number_][algorithm][:, outlier_ratio_index_],
                            color=colors[algorithm],
                            linestyle=linestyles[algorithm],
                            alpha=alpha
                            )
                else:
                    line,  = axs.plot(
                            self._noise_stds,
                            values_mean_[folder_name_][data_number_][algorithm][:, outlier_ratio_index_],
                            color=colors[algorithm],
                            linestyle=linestyles[algorithm],
                            alpha=alpha
                            )
                lines.append(line)
                legends.append(legend.format(algorithm))
            else:
                if std_:
                    line  = axs.errorbar(
                            self._noise_stds,
                            values_mean_[folder_name_][data_number_][algorithm][:, outlier_ratio_index_],
                            yerr=values_std_[folder_name][data_number_][algorithm][:, outlier_ratio_index_],
                            color=colors[algorithm],
                            linestyle=linestyles[algorithm],
                            alpha=alpha
                            )
                else:
                    line, =  axs.plot(
                            self._noise_stds,
                            values_mean_[folder_name_][data_number_][algorithm][:, outlier_ratio_index_],
                            color=colors[algorithm],
                            linestyle=linestyles[algorithm],
                            alpha=alpha
                            )
                lines.append(line)
                legends.append(legend.format(algorithm))
        axs.legend(lines, legends)
        if save_:
            plt.savefig(save_, bbox_inches='tight')
        else:
            plt.show()

        return None

    def generate_precision_image(self, data_number_, std_=True, save_=False, threshold_=None, remove_=[]):
        local_algorithms = [algo for algo in self._all_algorithms if algo not in remove_]
        self._generate_image(data_number_, local_algorithms, 'Precision', self._precision_mean, self._precision_std, std_, save_, threshold_)
        return None

    def generate_recall_image(self, data_number_, std_=True, save_=False, threshold_=None, remove_=[]):
        local_algorithms = [algo for algo in self._all_algorithms if algo not in remove_]
        self._generate_image(data_number_, local_algorithms, 'Recall', self._recall_mean, self._recall_std, std_, save_, threshold_)
        return None

    def generate_f1score_image(self, data_number_, std_=True, save_=False, threshold_=None, remove_=[]):
        local_algorithms = [algo for algo in self._all_algorithms_f1 if algo not in remove_]
        self._generate_image(data_number_, local_algorithms, 'F1-score', self._f1score_mean, self._f1score_std, std_, save_, threshold_)
        return None

    def generate_sigma_image(self, data_number_, std_=True, save_=False, threshold_=None, normalise_=True):
        if normalise_:
            local_algorithms = self._magsac_algorithms[:-1] + self._other_algorithms[1:]
            self._normalised_sigma_mean = self._value_container(local_algorithms)
            self._normalised_sigma_std = self._value_container(local_algorithms)
            for folder_name in self._folder_names:
                for data_number in self._data_numbers:
                    for algorithm in local_algorithms:
                        for noise_std_index, noise_std in enumerate(self._noise_stds):
                            for outlier_ratio_index, outlier_ratio in enumerate(self._outlier_ratios):
                                normalisation = noise_std
                                if normalisation == 0:
                                    normalisation = 1
                                self._normalised_sigma_mean[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = self._sigma_mean[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] / normalisation
                                self._normalised_sigma_std[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = self._sigma_std[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] / normalisation
            self._generate_image(data_number_, local_algorithms, 'Sigma', self._normalised_sigma_mean, self._normalised_sigma_std, std_, save_, threshold_)
        else:
            self._generate_image(data_number_, local_algorithms, 'Sigma', self._sigma_mean, self._sigma_std, save_, threshold_)
        return None

    def generate_time_image(self, data_number_, std_=True, save_=False, threshold_=None):
        self._generate_image(data_number_, self._all_time_algorithms, 'Time (s)', self._time_mean, self._time_std, std_, save_, threshold_)
        return None

    def generate_kept_image(self, data_number_, save_=False, remove_=[]):
        local_algorithms = [algo for algo in self._all_algorithms if algo not in remove_]
        self._kept_data = self._value_container(local_algorithms)
        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for algorithm in local_algorithms:
                    self._kept_data[folder_name][data_number][algorithm] = np.maximum(self._precision_kept[folder_name][data_number][algorithm], self._recall_kept[folder_name][data_number][algorithm])
        self._generate_image(data_number_, local_algorithms, 'Kept', self._kept_data, [], std_=False, save_=save_)
        return None

    def generate_kept_f1_image(self, data_number_, save_=False, remove_=[]):
        local_algorithms = [algo for algo in self._all_algorithms if algo not in remove_]
        self._kept_data = self._value_container(local_algorithms)
        for folder_name in self._folder_names:
            for data_number in self._data_numbers:
                for algorithm in local_algorithms:
                    self._kept_data[folder_name][data_number][algorithm] = np.maximum(self._precision_kept[folder_name][data_number][algorithm], self._recall_kept[folder_name][data_number][algorithm])
        self._generate_image(data_number_, local_algorithms, 'Kept', self._kept_data, [], std_=False, save_=save_)
        return None
