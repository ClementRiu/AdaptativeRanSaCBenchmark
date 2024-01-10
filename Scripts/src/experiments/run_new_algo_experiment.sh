#!/bin/bash

inFolderGlob=${inFolderGlob:--1}
outputFolder=${outputFolder:--1} #Datas will be written in this folder in a folder whose name is time of launch.
inputTimeFolder=${inputTimeFolder:--1} #Give a time folder here.

#Parameters of the experiments:

expModel=${expModel:-0}

expBegin=${expBegin:-1}
expEnd=${expEnd:--1}
expSeed=${expSeed:--1}

expUniform=${expUniform:-0}

expIterMax=${expIterMax:-10000}
expIterMaxStarSac=${expIterMaxStarSac:--1}

expMaxPrecision=${expMaxPrecision:-16}

expCpII=${expCpII:-0.95}

expNGen=${expNGen:-5}
expNRun=${expNRun:-5}

expForce=${expForce:--1}
expVerbose=${expVerbose:--1}

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi

  shift
done

#Generation of the paths and creation of the necessary folders:

## NEEDS TO BE CHANGED
if [ "${inFolderGlob}" == "-1" ];then
    inFolderGlob="/home/riuclement/Documents/USAC/data/"
    addInfo=1
fi

timeFolder=$(date +%Y_%m_%d_%H_%M_%S)

if [ "${outputFolder}" == "-1" ];then
    printf "\nOutput folder requiered (outputFolder):\n"
    exit 1
fi

if [ "${inputTimeFolder}" != "-1" ];then
    printf "\nReading data folder: ${outputFolder}${timeFolder}\n"
    timeFolder=${inputTimeFolder}
else
    printf "\nCreating data folder: ${outputFolder}${timeFolder}\n"
    mkdir -p ${outputFolder}${timeFolder}
fi

expFolder="${outputFolder}${timeFolder}/experiment_results/"
printf "\nCreating experiment folder: ${expFolder}"
mkdir -p ${expFolder}

if [ "${expModel}" == "0" ];then
    modelPrefix="hom"
    if [ "${addInfo}" == "1" ];then
        inFolderGlob="${inFolderGlob}homog/"
    fi
    # expEnd=10
fi
if [ "${expModel}" == "1" ];then
    modelPrefix="fun"
    if [ "${addInfo}" == "1" ];then
        inFolderGlob="${inFolderGlob}fundmatrix/"
    fi
    # expEnd=11
fi
if [ "${expModel}" == "2" ];then
    modelPrefix="ess"
    if [ "${addInfo}" == "1" ];then
        inFolderGlob="${inFolderGlob}essential/"
    fi
    # expEnd=6
fi


printf "\n\n---------- Generation done. Running experiment. ----------\n\n\n\n"

#Run of the experiment:

expStdNoiseMin=0.0
expStdNoiseMax=3.0
expStdNoiseStep=0.1

expStdNoiseCounterMin=0 #Used for the while loop.
expStdNoiseCounterMax=30 #Used for the while loop. = (expStdNoiseMax - expStdNoiseMin) / expStdNoiseStep

expStdNoise=${expStdNoiseMin}
expStdNoiseCounter=${expStdNoiseCounterMin} #Used for the while loop.

expOutlierRatioMin=.9
expOutlierRatioMax=0.9
expOutlierRatioStep=0.1

expOutlierRatioCounterMin=1 #Used for the while loop.
expOutlierRatioCounterMax=1 #Used for the while loop. = (expOutlierRatioMax - expOutlierRatioMin) / expOutlierRatioStep

expOutlierRatio=${expOutlierRatioMin}
expOutlierRatioCounter=${expOutlierRatioCounterMin} #Used for the while loop.

expCounter=${expBegin}

if [ "${expUniform}" == "0" ]; then
    outExpPathU="${expFolder}/uniform/"
else
    outExpPathU="${expFolder}/gaussian/"
fi
mkdir -p ${outExpPathU}

while [ "${expCounter}" -ge "${expBegin}" ] && [ "${expCounter}" -le "${expEnd}" ]
do
    while [ "${expOutlierRatioCounter}" -ge "${expOutlierRatioCounterMin}" ] && [ "${expOutlierRatioCounter}" -le "${expOutlierRatioCounterMax}" ]
    do
        while [ "${expStdNoiseCounter}" -ge "${expStdNoiseCounterMin}" ] && [ "${expStdNoiseCounter}" -le "${expStdNoiseCounterMax}" ]
        do
            expArgs="--iterMax ${expIterMax} --model-analyse ${expModel} --maximum-precision ${expMaxPrecision} --cpII ${expCpII} --num-gen-exp ${expNGen} --num-run ${expNRun}"

            if [ "${expForce}" -eq "1" ]; then
                expArgs="${expArgs} --force-compute"
            fi
            if [ "${expIterMaxStarSac}" -ge "1" ]; then
                expArgs="${expArgs} --iterMaxStarSac ${expIterMaxStarSac}"
            fi
            if [ "${expVerbose}" -eq "1" ]; then
                expArgs="${expArgs} --verbose"
            fi
            if [ "${expSeed}" -ge "0" ]; then
                expArgs="${expArgs} --time-seed ${expSeed}"
            fi
            if [ "${expOutlierRatioCounter}" == "0" ]; then
                expArgs="${expArgs} --ignore-outliers"
            fi

            inInfo=${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_info.txt
            inCalib=${inFolderGlob}test${expCounter}/calib_matrices.txt
            inNoisyInliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_NoisyIn.txt"
            inOutliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_Outliers.txt"

            outOutputInfo="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_facrOutput.txt"

            printf "Reading ${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_... .\n"
            printf "With parameters:\n\tNoise std: ${expStdNoise}\n\tOutlier ratio: ${expOutlierRatio}\n"

            execArgs="${expArgs} ${inInfo}"
            if [ "${expModel}" == "2" ];then
                execArgs="${execArgs} ${inCalib}"
            fi
            execArgs="${execArgs} ${inNoisyInliers} ${inOutliers} ${outOutputInfo}"

            ~/Documents/RANSAC-benchmark/Scripts/Build/experiments/new_algo_experiment ${execArgs}

            printf "Done.\n\n"

            expStdNoise="$(echo "${expStdNoise} + ${expStdNoiseStep}" | bc)"
            expStdNoiseCounter=$((${expStdNoiseCounter}+1))

        done # End of the std noise loop.

        expOutlierRatio="$(echo "${expOutlierRatio} + ${expOutlierRatioStep}" | bc)"
        expOutlierRatioCounter=$((${expOutlierRatioCounter}+1))

        expStdNoise=${expStdNoiseMin}
        expStdNoiseCounter=${expStdNoiseCounterMin}

    done # End of the outlier ratio loop.
    expCounter=$((${expCounter}+1))

    expOutlierRatio=${expOutlierRatioMin}
    expOutlierRatioCounter=${expOutlierRatioCounterMin}
done # End of the folder loop.
