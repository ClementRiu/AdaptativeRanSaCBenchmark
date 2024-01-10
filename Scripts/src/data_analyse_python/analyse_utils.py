import glob
import matplotlib.pyplot as plt
import numpy as np

def parse_path(path):
    parsed_path_std = path.split("_std")[-1]
    parsed_path_std = parsed_path_std.split("_")[0]

    parsed_path_ratio = path.split("_ratio")[-1]
    parsed_path_ratio = parsed_path_ratio.split("_")[0]

    return float(parsed_path_std), float(parsed_path_ratio)

def parse_values(values):
    return [[float(value) for value in string.rstrip().split()] for string in values]


def load_paths(folder_names, data_numbers, noise_stds, glob_path):
    paths = {folder_name: {data_number: {noise_std: []
                                         for noise_std in noise_stds}
                           for data_number in data_numbers}
             for folder_name in folder_names}
    for folder_name in folder_names:
        for data_number in data_numbers:
            for path in sorted(glob.glob(glob_path.format(folder_name, data_number))):
                noise_std, _ = parse_path(path)
                paths[folder_name][data_number][noise_std].append(path)

    return paths

def create_data_container(outlier_ratios, noise_stds, data_numbers, folder_names):
    return {folder_name: {data_number: {noise_std: {outlier_ratio: []
                                                    for outlier_ratio in outlier_ratios}
                                        for noise_std in noise_stds}
                          for data_number in data_numbers}
            for folder_name in folder_names}

def create_value_container(num_noise_stds, num_outlier_ratios, algorithms, data_numbers, folder_names):
    return {folder_name: {data_number: {algorithm: np.empty([num_noise_stds, num_outlier_ratios])
                                        for algorithm in algorithms}
                          for data_number in data_numbers}
            for folder_name in folder_names}

def summarize_data(folder_names, data_numbers, algorithms, noise_stds, outlier_ratios, data):
    values_mean = create_value_container(len(noise_stds), len(outlier_ratios), algorithms, data_numbers, folder_names)
    values_std = create_value_container(len(noise_stds), len(outlier_ratios), algorithms, data_numbers, folder_names)

    for folder_name in folder_names:
        for data_number in data_numbers:
            for algo_index, algorithm in enumerate(algorithms):
                for noise_std_index, noise_std in enumerate(noise_stds):
                    for outlier_ratio_index, outlier_ratio in enumerate(outlier_ratios):
                        try:
                            values = data[folder_name][data_number][noise_std][outlier_ratio][algo_index]
                            values_mean[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = np.mean(values)
                            values_std[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = np.std(values)
                        except IndexError:
                            values_mean[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = -1
                            values_std[folder_name][data_number][algorithm][noise_std_index, outlier_ratio_index] = 0
    return values_mean, values_std


def generate_image(data_numbers, outlier_ratios, algorithms, folder_names, noise_stds, value_name, values_mean, values_std, save=False):

    colors = ['b', 'r', 'g', 'c']
    if len(algorithms) == 2:
        colors = colors[1:3]
    if len(algorithms) == 3:
        colors = colors[0:3]
    linestyles = ['-', '--']
    legend = "{} algorithm for {} noise."
    xlabel_n = "Noise std value (in pixel)"
    xlabel_o = "Outlier ratio value"
    generic_title = "{} for outlier ratio of {} on dataset {}."
    generic_suptitle = "{} of " + "{} for only inlier, with noise, artificial homographies.".format(" ".join(algorithms))

    fig, axs = plt.subplots(len(data_numbers), len(outlier_ratios), figsize=(10 * len(outlier_ratios), 10 * len(data_numbers)), sharey=True)
    fig.suptitle(generic_suptitle.format(value_name))
    if (len(data_numbers) > 1):
        for data_index, data_number in enumerate(data_numbers):
            for outlier_ratio_index, outlier_ratio in enumerate(outlier_ratios):
                axs[data_index][outlier_ratio_index].set_title(generic_title.format(value_name, outlier_ratio, data_number))
                axs[data_index][outlier_ratio_index].set_xlabel(xlabel_n)
                axs[data_index][outlier_ratio_index].set_ylabel(value_name)
                lines = []
                legends = []
                for algo_index, algorithm in enumerate(algorithms):
                    for folder_index, folder_name in enumerate(folder_names):
                        line = axs[data_index][outlier_ratio_index].errorbar(noise_stds,
                                                                             values_mean[folder_name][data_number][algorithm][:, outlier_ratio_index],
                                                                             yerr=values_std[folder_name][data_number][algorithm][:, outlier_ratio_index],
                                                                             color=colors[algo_index],
                                                                             linestyle=linestyles[folder_index],
                                                                             alpha=0.5
                                                                            )
                        lines.append(line)
                        legends.append(legend.format(algorithm, folder_name))
                axs[data_index][outlier_ratio_index].legend(lines, legends)
    else:
        for data_number in data_numbers:
            for outlier_ratio_index, outlier_ratio in enumerate(outlier_ratios):
                axs[outlier_ratio_index].set_title(generic_title.format(value_name, outlier_ratio, data_number))
                axs[outlier_ratio_index].set_xlabel(xlabel_n)
                axs[outlier_ratio_index].set_ylabel(value_name)
                lines = []
                legends = []
                for algo_index, algorithm in enumerate(algorithms):
                    for folder_index, folder_name in enumerate(folder_names):
                        line = axs[outlier_ratio_index].errorbar(noise_stds,
                                                                 values_mean[folder_name][data_number][algorithm][:, outlier_ratio_index],
                                                                 yerr=values_std[folder_name][data_number][algorithm][:, outlier_ratio_index],
                                                                 color=colors[algo_index],
                                                                 linestyle=linestyles[folder_index]
                                                                )
                        lines.append(line)
                        legends.append(legend.format(algorithm, folder_name))
                axs[outlier_ratio_index].legend(lines, legends)
    if save:
        plt.savefig(save)
    else:
        plt.show()

def load_data(p_datas, r_datas, s_datas, t_datas, folder_names, data_numbers, noise_stds, paths):
    for folder_name in folder_names:
        for data_number in data_numbers:
            for noise_std in noise_stds:
                for path in paths[folder_name][data_number][noise_std]:
                    p_values = []
                    r_values = []
                    s_values = []
                    t_values = []
                    _, outlier_ratio = parse_path(path)
                    with open(path, 'r') as file:
                        counter = 0
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
                    p_datas[folder_name][data_number][noise_std][outlier_ratio] = parse_values(p_values)
                    r_datas[folder_name][data_number][noise_std][outlier_ratio] = parse_values(r_values)
                    s_datas[folder_name][data_number][noise_std][outlier_ratio] = parse_values(s_values)
                    t_datas[folder_name][data_number][noise_std][outlier_ratio] = parse_values(t_values)
    return None
