import glob
import numpy as np
from analyse_utils import parse_path
from magsac_data import Magsac_Data

def read_magsac(path, type_used=float):
    data = []
    with open(path, "r") as file:
        for line in file:
            element = np.array([value for value in line.split(", ")[:-1]]).astype(type_used)
            data.append(element)
    return data

def show_dimension(values):
    sizes = []
    for element in values:
        sizes.append(len(element))
    print(len(values), "x ", sizes)

def magsac_container(data_numbers, noise_stds, outlier_ratios):
    return {data_number: {noise_std: {outlier_ratio: []
                                      for outlier_ratio in outlier_ratios}
                          for noise_std in noise_stds}
            for data_number in data_numbers}

def load_magsac_data(data_numbers, noise_stds, outlier_ratios, file_names, file_types, file_index, default_path):
    paths = magsac_container(data_numbers, noise_stds, outlier_ratios)

    for data_number in data_numbers:
        for file_name in file_names:
            default_path_glob = default_path.format(data_number, file_name)
            for path in sorted(glob.glob(default_path_glob)):
                true_noise_std, true_outlier_ratio_ = parse_path(path)
                paths[data_number][true_noise_std][true_outlier_ratio_].append(path)

    magsac_data = magsac_container(data_numbers, noise_stds, outlier_ratios)

    for data_number in data_numbers:
        for noise_std in noise_stds:
            for outlier_ratio in outlier_ratios:
                # Add error catching for inexistant file.
                all_path_names = paths[data_number][noise_std][outlier_ratio]
                magsac_data[data_number][noise_std][outlier_ratio] = \
                [read_magsac(all_path_names[i], file_type) for i, file_type in enumerate(file_types)]

    magsac_data_for_process = magsac_container(data_numbers, noise_stds, outlier_ratios)

    for data_number in data_numbers:
        for noise_std in noise_stds:
            for outlier_ratio in outlier_ratios:
                data_for_process = magsac_data[data_number][noise_std][outlier_ratio]
                data_to_analyse = [Magsac_Data(data_for_process[file_index["Labels"]][i // 5],
                                               data_for_process[file_index["ErrorsAll"]][i],
                                               data_for_process[file_index["PosInl"]][i],
                                               data_for_process[file_index["Weights"]][i])
                                   for i in range(len(data_for_process[0]))]
                magsac_data_for_process[data_number][noise_std][outlier_ratio] = data_to_analyse

    return magsac_data_for_process

def compute_precision(true_positives, estimated_inliers):
    return len(true_positives) / len(estimated_inliers)

def compute_recall(true_positives, GT_inliers):
    return len(true_positives) / len(GT_inliers)

def load_magsac_paths(folder_names, data_numbers, noise_stds, outlier_ratios, file_names, glob_path):
    paths = {folder_name:
                {data_number:
                    {noise_std:
                        {outlier_ratio:
                            {file_name: ''
                            for file_name in file_names}
                        for outlier_ratio in outlier_ratios}
                    for noise_std in noise_stds}
                for data_number in data_numbers}
            for folder_name in folder_names}

    for folder_name in folder_names:
        for data_number in data_numbers:
            for file_name in file_names:
                for path in sorted(glob.glob(glob_path.format(folder_name, data_number, file_name))):
                    noise_std, outlier_ratio = parse_path(path)
                    paths[folder_name][data_number][noise_std][outlier_ratio][file_name] = path
    return paths
