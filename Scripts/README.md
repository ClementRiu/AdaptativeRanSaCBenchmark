# Improving 3D Reconstruction using Adaptative RanSaC and Benchmarking with Semi-Artificial Data Generation

Clement Riu & Pascal Monasse & Vincent Nozick - LIGM, Université Gustave

# Reviewed files in IPOL:
- src/demo/demo.cpp
- src/libOrsa/lrtsac.{hpp,cpp}

# Description:
Implementation of RANSAC algorithms:
- RANSAC (the classic one)
- AC-RANSAC (a.k.a. ORSA) minimizing Number of False Alarms and Fast-AC-RANSAC
- LRTSAC maximazing likelihood
- MAGSAC and MAGSAC++
- STARSAC
- MUSE

They can solve the following computer vision estimation problems:
- Homography transform
- Fundamental matrix
- Essential matrix
- Perspective from n Points

# Licensing:
See LICENSE.txt file

# Build
This project build relies on CMake (https://cmake.org/).

## LINUX/MAC:
```
  $ cd .../AutoRansac
  $ mkdir Build
  $ cd Build
  $ cmake -D CMAKE_BUILD_TYPE:string=Release ../src
  $ make
```
If you target to use an IDE to compile the code:
```
  $ cmake -G "CodeBlocks - Unix Makefiles" ../src
```

## WINDOWS:
  Launch cmake-gui.exe

  Fill the blank path.
  "Where is the source code :" (where the general CMakeLists.txt is).
    => Scripts/src
  "Where to build the binaries:" (where build project and object files will be).
    => Scripts/build

  Press Configure. (Select your IDE. ie, Visual Studio 10 Win64)
  Press Generate.
  Go to the merged_script/build path.
  Launch the Visual Studio solution and compile in release mode.

# USAGE
Go to the build folder. The executable files are in the `demo`, and `experiments` subfolders.

## Demo

There is a single demo executable file for all models and all available algorithms.
The `-m` or `--model` parameter controls the model and the `-a` or `-algo` parameter controls the algorithm used. Available models are: `Homography/Fundamental/Essential`, available algorithms are: `Ransac/AC-Ransac/fast-AC-Ransac/LRT`.

Usage:
```
./demo/demo [options] imgInA imgInB allInOutMatches.txt inlierOutMatches.txt [optImgOut]
- imgInA, imgInB: the two input images (JPG/PNG format)
- allInOutMatches.txt: output (input if -r) text file of format "x1 y1 x2 y2"
- inlierOutMatches.txt: output, but only with inliers.
        [optImgOut] (output images): inliers outliers [mosaic/epi [regA regB]]
- inliers, outliers: lines for inliers, lines for outliers and their error
- mosaic/epi: mosaic if Homography else epipolar lines
- regA, regB: registered images (Homography only)
        Options:
-c, --cut=ARG cut region of imagInA: wxh+x+y = rect [x,x+w] x [y,y+h] (0x0+0+0)
-r, --read Read file of matches allMatches.txt, do not use SIFT
-s, --sift=ARG SIFT distance ratio of descriptors (0.6)
-m, --model=ARG Model class: Homography/Fundamental/Essential (Homography)
-a, --algo=ARG Algorithm to use: Ransac/AC-Ransac/LRT (Ransac)
-i, --iterMax=ARG Number of iterations of the algorithm. (50000)
-p, --precision=ARG Max precision (in pixels) of registration (0=arbitrary) (3)
--cpIIT=ARG Confidence against type II error (Ransac/LRT) (0.99)
-n, --nModelMinRansac=ARG Min number of models before terminating ransac (1)
--cpI=ARG Confidence proba wrt type I error (LRT) (0.99)
--cpIIB=ARG Confidence proba wrt bailout (LRT) (0.95)
-k, --Kfile=ARG File for calibration matrices (Essential only)
-v, --verbose Print info during the run.
-t, --time-seed=ARG Use value instead of time for random seed. (1620308988)
```

Example of run: registration by homography of two images with LRTSAC, maximum allowable threshold 16 pixels. The algorithm puts a threshold at 2 pixels and finds 97% inliers.
```
./demo/demo ../../data/ramparts[12].jpg all.txt in.txt in.jpg out.jpg mosaic.jpg reg1.jpg reg2.jpg -a LRT -m Homography -v -p 16
Rerun with "-t 1620309519" to reproduce
sift:: 1st image: 629 keypoints
sift:: 2nd image: 695 keypoints
sift:: matches: 177
Remove 20/177 duplicate matches, keeping 157
Model: Homography, Algorithm: LRT
(init)  L=0.073914 inliers=0 precision=16 iter=-1/50000 Sigma={0.25...16}
  L=8.29127 inliers=154 precision=4 iter=0/27 Sigma={0.25...4}
  L=8.97558 inliers=154 precision=2.82843 iter=6/19 Sigma={0.25...2.82843}
  L=9.31486 inliers=150 precision=2 iter=7/16 Sigma={0.25...2}
  L=9.48396 inliers=152 precision=2 iter=10/15 Sigma={0.25...2}
Before refinement: RMSE/max error: 0.808628/1.56892
After  refinement: RMSE/max error: 0.63752/1.55203
Result=[ 1.09241 -0.109282 -162.474; 0.166758 1.04799 47.7185; 0.000210258 -3.15157e-05 1 ]
Inliers: 153/157=97%
Sigma: 2
Iterations: 15
Verif/model: 128.5
Runtime (s): 0.002387
-- Render Mosaic --
-- Render Mosaic - Image 1 --
-- Render Mosaic - Image 2 --
```

## Experiments

The `experiments` subfolder contain the files used to generate the data for ipol paper Automatic RANSAC by Likelihood Maximization https://www.ipol.im/pub/art/2022/357/, VISAPP paper https://www.scitepress.org/Link.aspx?doi=10.5220/0010873300003124 and extension https://link.springer.com/chapter/10.1007/978-3-031-45725-8_1.
There is also a subfolder `data_analyse_python` containing the python files used to analyse the data and create the figures present in the papers.

See specific REAMDE.md files of each subfolder for how to reproduce the experiments.


# Sources:
```
James V. Miller and Charles V. Stewart, “MUSE: Robust Surface Fitting using Unbiased Scale Estimates”, in: Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR) (Nov. 1996), pp. 300–306, doi: 10.1109/ CVPR.1996.517089.
```
```
Jongmoo Choi and Gerard Medioni, “StaRSaC: Stable random sample consensus for parameter estimation”, in: 2009 IEEE Conference on Computer Vision and Pattern Recognition (CVPR), June 2009, pp. 675–682, doi: 10.1109/CVPR.2009.5206678.
```
```
Lionel Moisan, Pierre Moulon, and Pascal Monasse, “Fundamental Matrix of a Stereo Pair, with A Contrario Elimination of Outliers”, in: Image Processing On Line (IPOL) 6 (2016), pp. 89–113.
&
Lionel Moisan and Bérenger Stival, “A probabilistic criterion to detect rigid point matches between two images and estimate the fundamental matrix”, in: Interna- tional Journal of Computer Vision (IJCV) 57.3 (2004), pp. 201–218.
&
Pierre Moulon, Pascal Monasse, and Renaud Marlet, “Adaptive structure from mo- tion with a contrario model estimation”, in: Proceedings ot the Asian Conference of Computer Vision (ACCV), Springer, 2012, pp. 257–270.
```
```
Andrea Cohen and Christopher Zach, “The Likelihood-Ratio Test and Efficient Robust Estimation”, in: Proceedings of the IEEE International Conference on Computer Vision (ICCV), Dec. 2015, pp. 2282–2290, doi: 10.1109/ICCV.2015.263.
```
```
Daniel Barath, Jiri Matas, and Jana Noskova, “MAGSAC: marginalizing sample consensus”, in: Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2019, pp. 10197–10205.
```
```
Daniel Barath, Jana Noskova, Maksym Ivashechkin, and Jiri Matas, “MAGSAC++, a fast, reliable and accurate robust estimator”, in: Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR), 2020, pp. 1304-1312.
```
