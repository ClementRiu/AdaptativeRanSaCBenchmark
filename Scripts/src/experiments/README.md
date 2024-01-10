# Experiments

There are multiple available executable files and associated bash scripts.

## The general experiment

There are two executable files and one bash script to handle a full experiment which will launch all available algorithms on all data available for a given model type with different noise level values and outlier ratios values:
* the `./generate_artificial_dataset` executable generates a set of semi-artificial inlier matches;
* the `./general_experiment` executable launches an experiment on a given set of semi-artificial inlier matches by adding noise, outliers and running all available algorithms;
* the `run_general_experiment.sh` bash script that runs a full experiment for a given model on USAC data (need to be downloaded separately).

In order to run a full experiment the `run_general_experiment.sh` bash script must be run for each model type.
**Before running this file `line 57` needs to be set to the actual path to USAC data and `line 142` and `line 233` need to be changed to the actual paths to the executable files.**


To use the bash script:
```
bash run_general_experiment.sh [options]
    Options:
--outputFolder [Requiered] The folder where to write output data.
--inputTimeFolder The folder where to write experiment data in outputFolder, if none given, it will be the timestamp of the experiment.
--genModel Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
--genPrecision Precision of AC-RANSAC for the experiment. Best left at 0. (0)
--genIterMax Maximum number of iterations allowed to AC-Ransac. (100000)
--genForce Run the generation even if a file with the same name already exists if set to 1. (0)
--genVerbose Print info during the generation if set to 1. (0)
--genSeed Use value instead of time for random seed (for debug). Set below 0 to use time. (-1)
--expIterMax Number of iterations of the algorithms. (10000)
--expPrecision Precision (in pixels) of registration of RANSAC. (3)
--expPrecisionRANSACbis Precision (in pixels) of registration of second RANSAC. (9)
--expMaxPrecision Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of AC-RANSAC and LRT. (16)
--expCpIITLRT Value of confidence with respect to early stopping for Ransac algorithm and for LRT. (0.99)
--expNModelMinRansac Min number of model to evaluate before terminating ransac (5)
--expCpILRT Confidence proba wrt type I error (LRT) (0.99)
--expCpIIBLRT Confidence proba wrt type I error (LRT) (0.95)
--expNGen Number of noisy datasets generated for the experiment. (5)
--expNRun Number of run of the algorithm. (5)
--expMaxMatches Maximum number of matches in the dataset. (4000)
--expUniform Use uniform noise - 0 or gaussian noise - 1. (0)
--expOutlierType Add outliers: 0 for no, 1 for uniform outliers, 2 for uniform outliers that can't be inliers. (2)
--expForce Run the experiment even if a file with the same name already exists. If left at -1, the same value as --genForce will be used. (-1)
--expVerbose Print info during the run. If left at -1 the same value as --genVerbose will be used. (-1)
```

Otherwise, the two executable files can be run individually:

To generate a dataset:
```
./generate_artificial_dataset [options] inMatches imgInA imgInB [calibrationMatrix] OutGoodMatches.txt outInfo.txt [imgWarped]/[imgOutInliers imgOutOutliers [imgOutEpi]]
- inMatches: file containing the matches to read.
- imgInA, imgInB: the two input images (JPG or PNG format)
- calibrationMatrix (only if Essential): one file containing the two calibration matrixes.
- OutGoodMatches.txt: output good matches text file of format "x1 y1 x2 y2"
- outInfo.txt: info about the run.
- imgOutInliers (optional): output image showing inliers if Fundamental or Essential
- imgOutOutliers (optional): output image showing outliers and their error if Fundamental or Essential
- imgOutEpi (optional): output epipolar image if Fundamental or Essential
- imgWarped (optional): registered image if homography
	Options:
-m, --model-analyse=ARG Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
-i, --iterMax=ARG Number of iterations of AC-RANSAC. (10000)
-p, --precision=ARG Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of RANSAC, AC-RANSAC and LRT. (0)
-f, --force-compute Run the computation even if an output file already exists.
-v, --verbose Print info during the run.
-t, --time-seed=ARG Use value instead of time for random seed (for debug). (time)
```

To run an experiment:
```
./general_experiment [options] good_matches.txt info.txt calib.txt outputInfo.txt [noisyInliers.txt [artificialOutliers.txt]]
- good_matches.txt: path to the good matches to read.
- info.txt: path to the run info (like dimensions of images).
- calib.txt: path to the calib matrix of both images.
- outputInfo.txt: path to the output file with all the informations on the experiment.
- noisyInliers.txt: path to the output file with the inliers after noise addition.
- artificialOutliers.txt: path to the output file with the artificial outliers.
	Options:
-i, --iterMax=ARG Number of iterations of algorithms. (10000)
-m, --model-analyse=ARG Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
--ransac-precision=ARG Precision (in pixels) of registration of RANSAC. (3)
--ransac-precision-bis=ARG Precision (in pixels) of registration of second RANSAC. (9)
--maximum-precision=ARG Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of AC-RANSAC and LRT. (16)
--cpIIT=ARG Value of confidence (Ransac/LRT) (0.99)
--nModelMinRansac=ARG Min number of model to evaluate before terminating ransac (5)
--cpI=ARG Confidence proba wrt type I error (LRT) (0.99)
--cpIIB=ARG Confidence proba wrt bailout (LRT) (0.95)
-g, --num-gen-exp=ARG Number of noisy datasets generated for the experiment. (1)
-n, --num-run=ARG Number of run of the algorithm. (1)
--max-match=ARG Maximum number of matches in the dataset. (0)
-s, --noise-level=ARG Value of the noise std. If > 0 noise will be added. (1)
-u, --uniform-noise Use uniform noise or gaussian noise.
-o, --outlier-type=ARG Add outliers: 0 for no, 1 for uniform outliers, 2 for uniform outliers that can't be inliers. (0)
-r, --outlier-ratio=ARG Ratio of outlier or number: if in [0, 1[ then it is the ratio, if integer equal or greater than 1 it is the number of outlier. (0)
-f, --force-compute Run the computation even if a file already exists.
-v, --verbose Print info during the run.
-t, --time-seed=ARG Use value instead of time for random seed (for debug). (time)
```

## The evaluate artificial experiment

This experiment computes a given algorithm on semi-artificial data.

* the `./evaluate_artificial` executable reads an inlier and an outlier file, mixes them and run the demo similarly to the `demo` file with addition of precision and recall computing;
* the `run_evaluate_artificial.sh` bash script runs `./evaluate_artificial` on a given dataset with user given inliers and outliers files and USAC data (need to be downloaded separately) with all available algorithms.

To run this experiment select outputs of the bailout quality experiment and run the experiment.  **Before running `run_evaluate_artificial.sh` `line 41` needs to be set to the actual path to USAC data and `line 157`, `line 159`, `line 161` and `line 163` need to be changed to the actual path to the executable file.**

To use the bash script:
```
bash run_evaluate_artificial.sh [options]
    Options:
--inputFileInlier [Requiered] Path to the input file containing the inliers.
--inputFileOutlier [Requiered] Path to the input file containing the outliers.
--inputFileNumber [Requiered] Dataset number the input files where generated from in USAC.
--outputFolder [Requiered] The folder where to write output data.
--expModel Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
--expIterMax Number of iterations of the algorithms. (10000)
--expPrecision Precision (in pixels) of registration of RANSAC. (3)
--expPrecisionRANSACbis Precision (in pixels) of registration of second RANSAC. (9)
--expMaxPrecision Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of AC-RANSAC and LRT. (16)
--expCpIITLRT Value of confidence with respect to early stopping for Ransac algorithm and for LRT. (0.99)
--expNModelMinRansac Min number of model to evaluate before terminating ransac (5)
--expCpILRT Confidence proba wrt type I error (LRT) (0.99)
--expCpIIBLRT Confidence proba wrt type I error (LRT) (0.95)
--expForce Run the experiment even if a file with the same name already exists. If left at -1, the same value as --genForce will be used. (-1)
--expVerbose Print info during the run. If left at -1 the same value as --genVerbose will be used. (-1)
--expSeed Use value instead of time for random seed (for debug). Set below 0 to use time. (-1)
```

Otherwise, the executable file can be run individually:
```
./evaluate_artificial [options] imgInA imgInB [calibrationMatrix] allInOutMatches.txt inlierOutMatches.txt [imgOutInliers imgOutOutliers [imgOutMosaic/imgOutEpi [imgOutMosaicA imgOutMosaicA]]
- imgInA, imgInB: the two input images (JPG or PNG format)
- calibrationMatrix (only if Essential): one file containing the two calibration matrixes.
- allInOutMatches.txt: output (input if option -r) text file of format "x1 y1 x2 y2"
- inlierOutMatches.txt: output, but only with inliers.
- imgOutInliers (optional): output image showing inliers
- imgOutOutliers (optional): output image showing outliers and their error
- imgOutMosaic/imgOutEpi (optional): output mosaic image if homography, output epipolar image if Fundamental or Essential
- imgOutMosaicA, imgOutMosaicA (optional): registered images if homography
	Options:
-c, --cut=ARG cut region of imagInA: wxh+x+y = rect [x,x+w] x [y,y+h] (0x0+0+0)
-s, --sift=ARG SIFT distance ratio of descriptors (0.6)
-m, --model-analyse=ARG Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
-u, --used-algo=ARG Algorithm to use: 0 for Ransac, 1 for AC-RANSAC, 2 for LRT. (0)
-g, --num-gen-exp=ARG Number of noisy datasets generated for the experiment. (1)
-i, --iterMax=ARG Number of iterations of the algorithm. (1000)
-p, --precision=ARG Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of RANSAC, AC-RANSAC and LRT. (3)
--cpIIT=ARG Confidence against type II error (Ransac/LRT) (0.99)
--nModelMinRansac=ARG Min number of model to evaluate before terminating ransac (5)
--cpI=ARG Confidence proba wrt type I error (LRT) (0.99)
--cpIIB=ARG Confidence proba wrt bailout (LRT) (0.95)
-f, --force-compute Run the computation even if a file already exists.
-v, --verbose Print info during the run.
-t, --time-seed=ARG Use value instead of time for random seed (for debug). (time)
```

## The Bailout quality experiment

The bailout quality experiment computes the performance of LRTSAC with and without bailout.

* the `./lrt_experiment` executable runs the LRTSac algorithm with and without early bailout on semi-artificial data generated from inliers of `./generate_artificial_dataset`;
* the `run_lrt_experiment.sh` bash script runs `./generate_artificial_dataset` and `./lrt_experiment` on its output.

In order to run a full experiment the `run_lrt_experiment.sh` bash script must be run for each model type.
**Before running `run_lrt_experiment.sh` `line 53` needs to be set to the actual path to USAC data and `line 138` and `line 229` need to be changed to the actual path to the executable file.**

To use the bash script:
```
bash run_lrt_experiment.sh [options]
    Options:
--outputFolder [Requiered] The folder where to write output data.
--inputTimeFolder The folder where to write experiment data in outputFolder, if none given, it will be the timestamp of the experiment.
--genModel Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
--genPrecision Precision of AC-RANSAC for the experiment. Best left at 0. (0)
--genIterMax Maximum number of iterations allowed to algorithms. (100000)
--genForce Run the generation even if a file with the same name already exists if set to 1. (0)
--genVerbose Print info during the generation if set to 1. (0)
--genSeed Use value instead of time for random seed (for debug). Set below 0 to use time. (-1)
--expIterMax Number of iterations of the algorithms. (10000)
--expMaxPrecision Max precision (in pixels) of registration of LRT. (16)
--expCpIITLRT Value of confidence with respect to early stopping for Ransac algorithm and for LRT. (0.99)
--expCpILRT Confidence proba wrt type I error (LRT) (0.99)
--expCpIIBLRT Confidence proba wrt type I error (LRT) (0.95)
--expNGen Number of noisy datasets generated for the experiment. (5)
--expNRun Number of run of the algorithm. (5)
--expMaxMatches Maximum number of matches in the dataset. (4000)
--expUniform Use uniform noise - 0 or gaussian noise - 1. (0)
--expOutlierType Add outliers: 0 for no, 1 for uniform outliers, 2 for uniform outliers that can't be inliers. (2)
--expForce Run the experiment even if a file with the same name already exists. If left at -1, the same value as --genForce will be used. (-1)
--expVerbose Print info during the run. If left at -1 the same value as --genVerbose will be used. (-1)
```

Otherwise, the executable file can be run individually:
```
./lrt_experiment [options] good_matches.txt info.txt calib.txt outputInfo.txt [noisyInliers.txt [artificialOutliers.txt]]
- good_matches.txt: path to the good matches to read.
- info.txt: path to the run info (like dimensions of images).
- calib.txt: path to the calib matrix of both images.
- outputInfo.txt: path to the output file with all the informations on the experiment.
- noisyInliers.txt: path to the output file with the inliers after noise addition.
- artificialOutliers.txt: path to the output file with the artificial outliers.
	Options:
-i, --iterMax=ARG Number of iterations of algorithms. (10000)
-m, --model-analyse=ARG Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
-p, --maximum-precision=ARG Max precision (in pixels) of registration (0=arbitrary for AC-RANSAC) of AC-RANSAC and LRT. (16)
--cpIIT=ARG Value of confidence (Ransac/LRT) (0.99)
--cpI=ARG Confidence proba wrt type I error (LRT) (0.99)
--cpIIB=ARG Confidence proba wrt to bailout (LRT) (0.95)
-g, --num-gen-exp=ARG Number of noisy datasets generated for the ipol experiment. (1)
-n, --num-run=ARG Number of run of the algorithm. (1)
--max-match=ARG Maximum number of matches in the dataset. (0)
-s, --noise-level=ARG Value of the noise std. If > 0 noise will be added. (1)
-u, --uniform-noise Use uniform noise or gaussian noise.
-o, --outlier-type=ARG Add outliers: 0 for no, 1 for uniform outliers, 2 for uniform outliers that can't be inliers. (0)
-r, --outlier-ratio=ARG Ratio of outlier or number: if in [0, 1[ then it is the ratio, if integer equal or greater than 1 it is the number of outlier. (0)
-f, --force-compute Run the computation even if a file already exists.
-v, --verbose Print info during the run.
-t, --time-seed=ARG Use value instead of time for random seed (for debug). (time)
```

## The time experiment

The time experiment evaluates the runtime of LRTSAC with different options enabled.

* the `./run_time_experiment` executable runs the LRTSac algorithm with different options enabled;
* the `time_experiment.sh` bash script runs `./run_time_experiment` on all USAC data for the given model type.

In order to run a full experiment the `time_experiment.sh` bash script must be run for each model type.
**Before running `time_experiment.sh` `line 35` needs to be set to the actual path to USAC data and `line 100` needs to be changed to the actual path to the executable file.**

To use the bash script:
```
bash run_lrt_experiment.sh [options]
    Options:
--outputFolder [Requiered] The folder where to write output data.
--expModel Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
--expIterMax Maximum number of iterations allowed to algorithms. (100000)
--expNRun Number of run of the algorithm. (100)
--expMaxPrecision Max precision (in pixels) of registration of LRT. (16)
--expCpIITLRT Value of confidence with respect to early stopping for Ransac algorithm and for LRT. (0.99)
--expCpILRT Confidence proba wrt type I error (LRT) (0.99)
--expCpIIBLRT Confidence proba wrt type I error (LRT) (0.95)
--expForce Run the experiment even if a file with the same name already exists. If left at -1, the same value as --genForce will be used. (-1)
--expVerbose Print info during the run. If left at -1 the same value as --genVerbose will be used. (-1)
--expSeed Use value instead of time for random seed (for debug). Set below 0 to use time. (-1)
```

Otherwise, the executable file can be run individually:
```
./time_experiment [options] imgInA imgInB [calibrationMatrix] allInOutMatches.txt inlierOutMatches.txt [imgOutInliers imgOutOutliers [imgOutMosaic/imgOutEpi [imgOutMosaicA imgOutMosaicA]]
- imgInA, imgInB: the two input images (JPG or PNG format)
- calibrationMatrix (only if Essential): one file containing the two calibration matrixes.
- allInOutMatches.txt: output (input if option -r) text file of format "x1 y1 x2 y2"
- inlierOutMatches.txt: output, but only with inliers.
- imgOutInliers (optional): output image showing inliers
- imgOutOutliers (optional): output image showing outliers and their error
- imgOutMosaic/imgOutEpi (optional): output mosaic image if homography, output epipolar image if Fundamental or Essential
- imgOutMosaicA, imgOutMosaicA (optional): registered images if homography
	Options:
-m, --model-analyse=ARG Model to analyse: 0 for Homography, 1 for Fundamental, 2 for Essential. (0)
-i, --iterMax=ARG Number of iterations of the algorithm. (1000)
-n, --num-run=ARG Number of run of the algorithm. (100)
-p, --precision=ARG Max precision (in pixels) of registration of LRT. (16)
--cpIIT=ARG Value of confidence proba for LRT (0.99)
--cpI=ARG Confidence proba wrt type I error (LRT) (0.99)
--cpIIB=ARG Confidence proba wrt bailout (LRT) (0.95)
-f, --force-compute Run the computation even if a file already exists.
-v, --verbose Print info during the run.
-t, --time-seed=ARG Use value instead of time for random seed (for debug). (time)
```
