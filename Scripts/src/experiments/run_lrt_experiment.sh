#!/bin/bash

outputFolder=${outputFolder:--1} #Datas will be written in this folder in a folder whose name is time of launch.
inputTimeFolder=${inputTimeFolder:--1} #To override generation of dataset, give a time folder here.

#Parameters of the generation:

genModel=${genModel:-0}

genPrecision=0
genIterMax=100000

genForce=${genForce:-0}
genVerbose=${genVerbose:-0}
genSeed=${genSeed:--1}


#Parameters of the experiments:

expIterMax=${expIterMax:-10000}


expMaxPrecision=${expMaxPrecision:-16}
expCpIIT=${expCpIIT:-0.99}

expCpILRT=${expCpILRT:-0.99}
expCpIIBLRT=${expCpIIBLRT:-0.95}

expNGen=${expNGen:-5}
expNRun=${expNRun:-5}
expMaxMatches=${expMaxMatches:-4000}
# expStdNoise=${expStdNoise:-1.0} # Those values are going to variate for the evaluation.
expUniform=${expUniform:-0}
expOutlierType=${expOutlierType:-2}
# expOutlierRatio=${expOutlierRatio:-0.0} # Those values are going to variate for the evaluation.

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
inFolderGlob="/home/riuclement/Documents/USAC/data/"
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

expFolder="${outputFolder}${timeFolder}/experience_results/"
printf "\nCreating experience folder: ${expFolder}"
mkdir -p ${expFolder}

genBegin=1
if [ "${genModel}" == "0" ];then
    modelPrefix="hom"
    inFolderGlob="${inFolderGlob}homog/"
    genEnd=10
fi
if [ "${genModel}" == "1" ];then
    modelPrefix="fun"
    inFolderGlob="${inFolderGlob}fundmatrix/"
    genEnd=11
fi
if [ "${genModel}" == "2" ];then
    modelPrefix="ess"
    inFolderGlob="${inFolderGlob}essential/"
    genEnd=6
fi

#Generation of the dataset:
genCounter=${genBegin}

if [ "${inputTimeFolder}" = "-1" ];then

    genArgs="--model-analyse ${genModel} --iterMax ${genIterMax} --precision ${genPrecision}"
    if [ "${genForce}" -eq "1" ]; then
        genArgs="${genArgs} --force-compute"
    fi
    if [ "${genVerbose}" -eq "1" ]; then
        genArgs="${genArgs} --verbose"
    fi
    if [ "${genSeed}" -ge "0" ]; then
        genArgs="${genArgs} --time-seed ${genSeed}"
    fi

    while [ "${genCounter}" -ge "${genBegin}" ] && [ "${genCounter}" -le "${genEnd}" ]
    do
        inFolder=${inFolderGlob}test${genCounter}/
        inMatches=${inFolder}orig_pts.txt
        inImg1=${inFolder}im1.jpg
        inImg2=${inFolder}im2.jpg
        inCalib=${inFolder}calib_matrices.txt

        outGM=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_gm.txt
        outInfo=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_info.txt

        outInlierImg=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_inlier.jpg
        outOutlierImg=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_outlier.jpg

        outWarpedImg=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_warped.jpg
        outEpiImg=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_epi.jpg


        printf "\nReading ${inFolder} .\n"

        execArgs="${genArgs} ${inMatches} ${inImg1} ${inImg2}"
        if [ "${genModel}" == "2" ];then
            execArgs="${execArgs} ${inCalib}"
        fi
        execArgs="${execArgs} ${outGM} ${outInfo}"
        if [ "${genModel}" == "0" ];then
            execArgs="${execArgs} ${outWarpedImg}"
        fi
        if [ "${genModel}" -ge "1" ];then
            execArgs="${execArgs} ${outInlierImg} ${outOutlierImg} ${outEpiImg}"
        fi

        ~/Documents/RANSAC-benchmark/Scripts/Build/experiments/generate_artificial_dataset ${execArgs}

        printf "Done.\n\n"

         genCounter=$((${genCounter}+1))
    done

fi

printf "\n\n---------- Generation done. Running experiment. ----------\n\n\n\n"

#Run of the experiment:
expModel=${genModel} # Forced equal to the generation parameter.
if [ "${expForce}" == "-1" ]; then
    expForce="${genForce}"
fi
if [ "${expVerbose}" == "-1" ]; then
    expVerbose="${genVerbose}"
fi
expSeed=${genSeed} # Forced equal to the generation parameter.

expStdNoiseMin=0.0
expStdNoiseMax=3.0
expStdNoiseStep=0.1

expStdNoiseCounterMin=0 #Used for the while loop.
expStdNoiseCounterMax=30 #Used for the while loop. = (expStdNoiseMax - expStdNoiseMin) / expStdNoiseStep

expStdNoise=${expStdNoiseMin}
expStdNoiseCounter=${expStdNoiseCounterMin} #Used for the while loop.

expOutlierRatioMin=0.0
expOutlierRatioMax=0.9
expOutlierRatioStep=0.1

expOutlierRatioCounterMin=0 #Used for the while loop.
expOutlierRatioCounterMax=9 #Used for the while loop. = (expOutlierRatioMax - expOutlierRatioMin) / expOutlierRatioStep

expOutlierRatio=${expOutlierRatioMin}
expOutlierRatioCounter=${expOutlierRatioCounterMin} #Used for the while loop.

expBegin=${genBegin}
expEnd=${genEnd}

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
            expArgs="--iterMax ${expIterMax} --model-analyse ${expModel} --maximum-precision ${expMaxPrecision} --cpIIT ${expCpIIT} --cpI ${expCpILRT} --cpIIB ${expCpIIBLRT} --num-gen-exp ${expNGen} --num-run ${expNRun} --max-match ${expMaxMatches} --noise-level ${expStdNoise} --outlier-type ${expOutlierType} --outlier-ratio ${expOutlierRatio}"

            if [ "${expForce}" -eq "1" ]; then
                expArgs="${expArgs} --force-compute"
            fi
            if [ "${expVerbose}" -eq "1" ]; then
                expArgs="${expArgs} --verbose"
            fi
            if [ "${expSeed}" -ge "0" ]; then
                expArgs="${expArgs} --time-seed ${expSeed}"
            fi
            if [ "${expUniform}" == "0" ]; then
                expArgs="${expArgs} --uniform-noise"
            fi

            inGM=${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_gm.txt
            inInfo=${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_info.txt
            inCalib=${inFolderGlob}test${expCounter}/calib_matrices.txt

            outOutputInfo="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_output_lrt.txt"
            outNoisyInliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_NoisyIn_lrt.txt"
            outOutliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_Outliers_lrt.txt"

            printf "Reading ${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_... .\n"
            printf "With parameters:\n\tNoise std: ${expStdNoise}\n\tOutlier ratio: ${expOutlierRatio}\n"

            execArgs="${expArgs} ${inGM} ${inInfo}"
            if [ "${genModel}" == "2" ];then
                execArgs="${execArgs} ${inCalib}"
            fi
            execArgs="${execArgs} ${outOutputInfo} ${outNoisyInliers} ${outOutliers}"

            ~/Documents/RANSAC-benchmark/Scripts/Build/experiments/lrt_experiment ${execArgs}

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
