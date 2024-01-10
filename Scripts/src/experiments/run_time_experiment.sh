#!/bin/bash

outputFolder=${outputFolder:--1} #Datas will be written in this folder in a folder whose name is time of launch.

#Parameters of the ipol time experience:
expModel=${expModel:-0}

expIterMax=${expIterMax:-20000}
expNRun=${expNRun:-100}


expMaxPrecision=${expMaxPrecision:-16}
expCpIIT=${expCpIIT:-0.99}

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

if [ "${outputFolder}" == "-1" ];then
    printf "\nOutput folder requiered (outputFolder):\n"
    exit 1
fi

printf "\nCreating data folder: ${outputFolder}${timeFolder}\n"
mkdir -p ${outputFolder}${timeFolder}

# expFolder="${outputFolder}${timeFolder}/time_exp_results/"
# printf "\nCreating experience folder: ${expFolder}"
# mkdir -p ${expFolder}

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
expCounter=${expBegin}

expArgs="--model-analyse ${expModel} --iterMax ${expIterMax} --num-run ${expNRun} --precision ${expMaxPrecision} --cpIIT ${expCpIIT} --cpI ${expCpILRT} --cpIIB ${expCpIIBLRT}"
if [ "${expForce}" -eq "1" ]; then
    expArgs="${expArgs} --force-compute"
fi
if [ "${expVerbose}" -eq "1" ]; then
    expArgs="${expArgs} --verbose"
fi
if [ "${expSeed}" -ge "0" ]; then
    expArgs="${expArgs} --time-seed ${expSeed}"
fi

while [ "${expCounter}" -ge "${expBegin}" ] && [ "${expCounter}" -le "${expEnd}" ]
do
    inFolder=${inFolderGlob}test${expCounter}/
    inMatches=${inFolder}orig_pts.txt
    inImg1=${inFolder}im1.jpg
    inImg2=${inFolder}im2.jpg
    inCalib=${inFolder}calib_matrices.txt

    outMetrics=${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_TimeMetrics.txt


    printf "\nReading ${inFolder} .\n"

    execArgs="${expArgs} ${inImg1} ${inImg2}"
    if [ "${expModel}" == "2" ];then
        execArgs="${execArgs} ${inCalib}"
    fi
    execArgs="${execArgs} ${inMatches} ${outMetrics}"

    ~/Documents/RANSAC-benchmark/Scripts/Build/experiments/time_experiment ${execArgs}

    printf "Done.\n\n"

     expCounter=$((${expCounter}+1))
done
