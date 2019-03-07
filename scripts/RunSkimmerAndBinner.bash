#!/bin/bash

#This runs the calls the RunSkimmerAndBinner.C macro which loads the necessary libraries
#and then runs skimmerAndBinner.cxx. The code is run for each root file in the data
# directory.

###########################################################
#SET THE DATA DIRECTORY HERE
dataDirectory=$1
#/scratch_menkar/FixedTargetData/AuAu_4_5GeV_2015/
#dataDirectory=/home/kinghorn/rootfiles/

#SET THE OUTPUT DIRECTORY HERE
outputDirectory=$2
#/home/sheppelmann/dev/davisdstanalysis/userfiles/AuAu_4_5GeV_2015/analysisCleanup/

#SET THE NUMBER OF EVENTS HERE (USE -1 FOR ALL)
nEvents=$3
#-1
#nEvents=1000

#SET THE STAR LIBRARY VERSION
starlib=$4
#SL16a

#SET THE COLLISION ENERGY HERE
energy=$5
#4.5

#SET THE EVENT CONFIGURATION HERE
eventConfig=$6
#FixedTarget2015

#tof file with gaussian fits
tofFile=$7
#../userfiles/AuAu_4_5GeV_2015/analysisCleanup/ztof_GausFits.root

#TOF Cut 
tofCut=$8
#1

########################################################### 

#Array containing all of the dataFiles
dataFiles=( $dataDirectory/*.root )
#dataFiles=( $dataDirectory/AlAu_4_9GeV_2015_0.root  )

processID=()
numberOfFiles=${#dataFiles[@]}
outFiles=()

for i in ${dataFiles[@]}
do
    echo "Running on dataFile: " $i

    outFile=$(basename $i .root)
    outFile=$outputDirectory/"$outFile"_Processed.root
    echo $outFile

    outFiles+=($outFile)

    root -l -q -b ../macros/RunSkimmerAndBinner.C\(\"$i\",\"$starlib\",$nEvents,\"$outFile\",\"$tofFile\",$tofCut,$energy,\"$eventConfig\"\) 

    processID+=($!)
    echo ${processID[@]}
done
wait ${processID[@]}

if [ $tofCut -eq 1 ]
then
    hadd $outputDirectory/YieldHistograms_tofCut.root ${outFiles[@]}
fi

if [ $tofCut -eq 0 ]
then
    hadd $outputDirectory/YieldHistograms.root ${outFiles[@]}
fi

wait

rm ${outFiles[@]}


exit
