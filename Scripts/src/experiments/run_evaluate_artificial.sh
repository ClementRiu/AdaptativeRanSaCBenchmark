#!/bin/bash

inputFileInlier=${inputFileInlier:--1} #Inliers will be read from this file.
inputFileOutlier=${inputFileOutlier:--1} #Outliers will be read from this file.
inputFileNumber=${inputFileNumber:-1} #Outliers will be read from this file.
outputFolder=${outputFolder:--1} #Datas will be written in this folder in a folder whose name is time of launch.

#Parameters of the ipol time experience:
expModel=${expModel:-0}
genNGen=${genNGen:-1}

expIterMax=${expIterMax:-10000}

expPrecision=${expPrecision:-3}
expPrecisionRANSACbis=${expPrecisionRANSACbis:-9}
expMaxPrecision=${expMaxPrecision:-16}
expCpIIT=${expCpIIT:-0.99}

expNModelMinRansac=${expNModelMinRansac:-5}

expCpILRT=${expCpILRT:-0.99}
expCpIIBLRT=${expCpIIBLRT:-0.95}

expForce=${expForce:-0}
expVerbose=${expVerbose:-0}
expSeed=${expSeed:--1}

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi

  shift
done

#Generation of the paths and creation of the necessary folders:

## NEEDS TO BE CHANGED
inFolderGlob="/home/riuclement/Documents/USAC/data/"
timeFolder=$(date +%Y_%m_%d_%H_%M_%S)

if [ "${inputFileInlier}" == "-1" ];then
    printf "\nInput file requiered for inliers (inputFileInlier):\n"
    exit 1
fi
if [ "${inputFileOutlier}" == "-1" ];then
    printf "\nInput file requiered for outliers (inputFileOutlier):\n"
    exit 1
fi
if [ "${outputFolder}" == "-1" ];then
    printf "\nOutput folder requiered (outputFolder):\n"
    exit 1
fi

outFolder="${outputFolder}${timeFolder}/"
printf "\nCreating data folder: ${outFolder}\n"
mkdir -p ${outFolder}

RansacPath="${outFolder}Ransac_"

RansacBisPath="${outFolder}RansacBis_"

ACRansacPath="${outFolder}ACRansac_"

LRTSacPath="${outFolder}LRTsac_"

expBegin=1
if [ "${expModel}" == "0" ];then
    modelPrefix="hom"
    inFolderGlob="${inFolderGlob}homog/"
    expEnd=10
fi
if [ "${expModel}" == "1" ];then
    modelPrefix="fun"
    inFolderGlob="${inFolderGlob}fundmatrix/"
    expEnd=11
fi
if [ "${expModel}" == "2" ];then
    modelPrefix="ess"
    inFolderGlob="${inFolderGlob}essential/"
    expEnd=6
fi

#Generation of the dataset:
expArgs="--model-analyse ${expModel} --num-gen-exp ${genNGen} --iterMax ${expIterMax} --cpIIT ${expCpIIT}"
if [ "${expForce}" -eq "1" ]; then
    expArgs="${expArgs} --force-compute"
fi
if [ "${expVerbose}" -eq "1" ]; then
    expArgs="${expArgs} --verbose"
fi
if [ "${expSeed}" -ge "0" ]; then
    expArgs="${expArgs} --time-seed ${expSeed}"
fi

expArgsRansac="${expArgs} --used-algo 0 --precision ${expPrecision} --nModelMinRansac ${expNModelMinRansac}"
expArgsRansacBis="${expArgs} --used-algo 0 --precision ${expPrecisionRANSACbis} --nModelMinRansac ${expNModelMinRansac}"
expArgsACRansac="${expArgs} --used-algo 1 --precision ${expMaxPrecision}"
expArgsLRTSac="${expArgs} --used-algo 2 --precision ${expMaxPrecision} --cpI ${expCpILRT} --cpIIB ${expCpIIBLRT}"


inFolder=${inFolderGlob}test${inputFileNumber}/
inImg1=${inFolder}im1.jpg
inImg2=${inFolder}im2.jpg
inCalib=${inFolder}calib_matrices.txt

outGoodMatches=demoIpol_${modelPrefix}_good_pts.txt

outInlierMatches=demoIpol_${modelPrefix}_inliersMatch.jpg
outOutlierMatches=demoIpol_${modelPrefix}_outliersMatch.jpg

outEpipolar=demoIpol_${modelPrefix}_epipolar.jpg

outMosaic0=demoIpol_${modelPrefix}_mosaic0.jpg
outMosaic1=demoIpol_${modelPrefix}_mosaic1.jpg
outMosaic2=demoIpol_${modelPrefix}_mosaic2.jpg

printf "\nReading ${inFolder} .\n"

execArgs="${inImg1} ${inImg2}"
if [ "${expModel}" == "2" ];then
    execArgs="${execArgs} ${inCalib}"
fi
execArgs="${execArgs} ${inputFileInlier} ${inputFileOutlier}"

execArgsRansac="${expArgsRansac} ${execArgs} ${RansacPath}${outGoodMatches} ${RansacPath}${outInlierMatches} ${RansacPath}${outOutlierMatches}"
if [ "${expModel}" == "0" ];then
    execArgsRansac="${execArgsRansac} ${RansacPath}${outMosaic0} ${RansacPath}${outMosaic1} ${RansacPath}${outMosaic2}"
else
    execArgsRansac="${execArgsRansac} ${RansacPath}${outEpipolar}"
fi

execArgsRansacBis="${expArgsRansacBis} ${execArgs} ${RansacBisPath}${outGoodMatches} ${RansacBisPath}${outInlierMatches} ${RansacBisPath}${outOutlierMatches}"
if [ "${expModel}" == "0" ];then
    execArgsRansacBis="${execArgsRansacBis} ${RansacBisPath}${outMosaic0} ${RansacBisPath}${outMosaic1} ${RansacBisPath}${outMosaic2}"
else
    execArgsRansacBis="${execArgsRansacBis} ${RansacBisPath}${outEpipolar}"
fi

execArgsACRansac="${expArgsACRansac} ${execArgs} ${ACRansacPath}${outGoodMatches} ${ACRansacPath}${outInlierMatches} ${ACRansacPath}${outOutlierMatches}"
if [ "${expModel}" == "0" ];then
    execArgsACRansac="${execArgsACRansac} ${ACRansacPath}${outMosaic0} ${ACRansacPath}${outMosaic1} ${ACRansacPath}${outMosaic2}"
else
    execArgsACRansac="${execArgsACRansac} ${ACRansacPath}${outEpipolar}"
fi

execArgsLRTSac="${expArgsLRTSac} ${execArgs} ${LRTSacPath}${outGoodMatches} ${LRTSacPath}${outInlierMatches} ${LRTSacPath}${outOutlierMatches}"
if [ "${expModel}" == "0" ];then
    execArgsLRTSac="${execArgsLRTSac} ${LRTSacPath}${outMosaic0} ${LRTSacPath}${outMosaic1} ${LRTSacPath}${outMosaic2}"
else
    execArgsLRTSac="${execArgsLRTSac} ${LRTSacPath}${outEpipolar}"
fi


~/Documents/RANSAC-benchmark/Scripts/Build/experiments/evaluate_artificial ${execArgsRansac}
printf "\n"
~/Documents/RANSAC-benchmark/Scripts/Build/experiments/evaluate_artificial ${execArgsRansacBis}
printf "\n"
~/Documents/RANSAC-benchmark/Scripts/Build/experiments/evaluate_artificial ${execArgsACRansac}
printf "\n"
~/Documents/RANSAC-benchmark/Scripts/Build/experiments/evaluate_artificial ${execArgsLRTSac}
printf "\n"

printf "Done.\n\n"
