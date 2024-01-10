#!/bin/bash

inFolderGlob=${inFolderGlob:--1} #If folder needs to be changed.
outputFolder=${outputFolder:--1} #Datas will be written in this folder in a folder whose name is time of launch.
inputTimeFolder=${inputTimeFolder:--1} #To override generation of dataset, give a time folder here.

#Parameters of the generation:

genModel=${genModel:-0}

genPrecision=0
genIterMax=100000

genEnd=${genEnd:--1}
genForce=${genForce:-0}
genVerbose=${genVerbose:-0}
genSeed=${genSeed:--1}


#Parameters of the experiments:

expIterMax=${expIterMax:-10000}

expPrecision=${expPrecision:-3}
expPrecisionRANSACbis=${expPrecisionRANSACbis:-9}
expMaxPrecision=${expMaxPrecision:-16}
expCpIITLRT=${expCpIITLRT:-0.99}

expNModelMinRansac=${expNModelMinRansac:-5}

expCpILRT=${expCpILRT:-0.99}
expCpIIBLRT=${expCpIIBLRT:-0.95}

expThresholdMagsac=${expThresholdMagsac:-10}
expCutoffMagsac=${expCutoffMagsac:-2}
expPartitionMagsac=${expPartitionMagsac:-10}
expTimeMagsac=${expTimeMagsac:-2}

expNGen=${expNGen:-5}
expNRun=${expNRun:-5}
expMaxMatches=${expMaxMatches:-4000}
# expStdNoise=${expStdNoise:-1.0} # Those values are going to variate for the evaluation.
expIncreaseRatio=${expIncreaseRatio:-1.1}
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
if [ "${genModel}" == "0" ];then
    modelPrefix="pnp"
fi
if [ "${inFolderGlob}" == "-1" ];then
    inFolderGlob="/home/riuclement/Documents/datasets/megadepth_sfm/pnp_ready/"
    # inFolderGlob="${inFolderGlob}${modelPrefix}/"
    genEnd=16
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

genBegin=1

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
        inMatches=${inFolder}2D3D_pts.txt
        inImg=${inFolder}im.jpg
        inCalib=${inFolder}calib_matrix.txt

        outGM=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_gm.txt
        outInfo=${outputFolder}${timeFolder}/${modelPrefix}_${genCounter}_info.txt

        printf "\nReading ${inFolder} .\n"

        execArgs="${genArgs} ${inMatches} ${inImg} ${inCalib} ${outGM} ${outInfo}"

        ~/Documents/RANSAC-benchmark/Scripts/Build/experiments/pnp_generate_artificial_dataset ${execArgs}

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
            expArgs="--iterMax ${expIterMax} --model-analyse ${expModel} --ransac-precision ${expPrecision} --ransac-precision-bis ${expPrecisionRANSACbis} --maximum-precision ${expMaxPrecision} --cpIIT ${expCpIITLRT} --nModelMinRansac ${expNModelMinRansac} --cpI ${expCpILRT} --cpIIB ${expCpIIBLRT} --threshold-magsac ${expThresholdMagsac} --cutoff-magsac ${expCutoffMagsac} --partition-number ${expPartitionMagsac} --maxTime ${expTimeMagsac} --num-gen-exp ${expNGen} --num-run ${expNRun} --max-match ${expMaxMatches} --noise-level ${expStdNoise} --outlier-type ${expOutlierType} --outlier-ratio ${expOutlierRatio} --increase-ratio ${expIncreaseRatio}"

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
            inCalib=${inFolderGlob}test${expCounter}/calib_matrix.txt

            outOutputInfo="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_output.txt"
            outNoisyInliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_NoisyIn.txt"
            outOutliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_Outliers.txt"

            outMagsacLabels="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacLabels.txt"
            outMagsacWeights="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacWeights.txt"
            outMagsacInliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacInliers.txt"
            outMagsacInErrors="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacInErrors.txt"
            outMagsacFullErrors="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacFullErrors.txt"
            outMagsacPPWeights="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacPPWeights.txt"
            outMagsacPPInliers="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacPPInliers.txt"
            outMagsacPPInErrors="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacPPInErrors.txt"
            outMagsacPPFullErrors="${outExpPathU}/${modelPrefix}_${expCounter}_std${expStdNoise}_ratio${expOutlierRatio}_outMagsacPPFullErrors.txt"

            printf "Reading ${outputFolder}${timeFolder}/${modelPrefix}_${expCounter}_... .\n"
            printf "With parameters:\n\tNoise std: ${expStdNoise}\n\tOutlier ratio: ${expOutlierRatio}\n"

            execArgs="${expArgs} ${inGM} ${inInfo} ${inCalib} ${outOutputInfo} ${outNoisyInliers} ${outOutliers}"
            execArgs="${execArgs} ${outMagsacLabels} ${outMagsacWeights} ${outMagsacInliers} ${outMagsacInErrors} ${outMagsacFullErrors} ${outMagsacPPWeights} ${outMagsacPPInliers} ${outMagsacPPInErrors} ${outMagsacPPFullErrors}"

            ~/Documents/RANSAC-benchmark/Scripts/Build/experiments/pnp_experiment ${execArgs}

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
