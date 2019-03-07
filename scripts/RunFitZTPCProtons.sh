#!/bin/bash

#This calls the RunFitZTPCProtons.C file in the macros directory which loads the necessary
#libraries and then runs fitZTPCProtons.cxx in the src/analysis directory.

######################################################
#SET THE YIELD HISTOGRAM FILE NAME AND PATH (FULL PATH)
yieldHistoFile=$1
#yieldHistoFile=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/YieldHistograms.root

#SET THE SPECTRA FILE NAME AND PATH (FULL PATH)
spectraFile=$2
#spectraFile=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/zTPCProtonSpectra.root

#SET THE PROTONWIDTH FILE NAME AND PATH (FULL PATH)
protwidth=$3
#protwidth=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/zTPCProtonsTofCut.root

#SET THE DEUTERONWIDTH FILE NAME AND PATH (FULL PATH)
deutwidth=$4
#deutwidth=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/zTPCDeuteronsTofCut.root

#SET THE RATIO FILE
ampRatioFile=$5
#ampRatioFile=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/proton_Tof_spectra.root

#SET THE OUTPUT IMAGE DIRECTORY (LEAVE BLANK TO SAVE NO IMAGES)
imageDir=$6
#imageDir=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/zTPCImages

#SET THE RAPIDITY VALUE (USE -999 FOR ALL)
rapidityValue=-999

#SET THE CENTRALITY BIN INDEX (USE -1 FOR ALL)
#centralityIndex=-1
centralityIndex=-1

#SET START/END RAPIDITY BIN
start=7;
#start=0;
last=19;
#last=0;

#SET THE STARLIBRARY VERSION
starlib=$8
######################################################

if [ "$centralityIndex" = "-1" ]; then

outFiles=()
processID=()


for i in `seq 0 5`;
do
    outFile=$(basename $spectraFile .root)
    outFile="$outFile"_Cent"$i".root
    echo $outFile
    outFiles+=($outFile)
    
    root -l -b -q ../macros/RunFitZTPCProtons.C\(\"$yieldHistoFile\",\"$outFile\",\"$protwidth\",\"$deutwidth\",\"$ampRatioFile\",\"$starlib\",$i,$rapidityValue,\"$imageDir\"\,$start,$last\) #> /dev/null 2>&1 &

  processID+=($!)
    echo ${processID[@]}
done
wait ${processID[@]}

hadd -f $outputDirectory/$spectraFile ${outFiles[@]}

wait

#rm ${outFiles[@]}

else

root -l -b -q ../macros/RunFitZTPCProtons.C\(\"$yieldHistoFile\",\"$spectraFile\",\"$protwidth\",\"$deutwidth\",\"$ampRatioFile\",\"$starlib\",$centralityIndex,$rapidityValue,\"$imageDir\"\,$start,$last\) #> /dev/null 2>&1 &


fi

exit
