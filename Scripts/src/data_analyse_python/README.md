# Data analyse

To use the python files in this folder, install the required packages with `pip install -r requirements.txt`.

This folder contains files related to the analyse of the outputs of the different experiment files.
Only IPOL related files are ready for use.

* `IPOL_LRT_likelihood.ipynb`: is an analysis of the behavior of the likelihood of LRTSAC to better understand its results;
* `IPOL_LRT_precision_recall.ipynb`: is an analysis of the results of `run_lrt_experiment.sh` file.
* `IPOL_LRT_time_options.ipynb`: is an analysis of the results of `run_time_experiment.sh` file.

The `single_data_loader.py` is used by those file to read the data.
