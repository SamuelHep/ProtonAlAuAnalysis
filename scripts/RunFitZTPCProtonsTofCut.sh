#!/bin/bash

#This calls the RunFitZTPCProtonsTofCut.C file in the macros directory which loads the necessary libraries
######################################################
#SET THE YIELD HISTOGRAM FILE NAME AND PATH (FULL PATH)
yieldHistoFile=$1
#yieldHistoFile=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/YieldHistograms.root

#SET THE OUTFILE FILE NAME AND PATH (FULL PATH)
outRoot=$2
#outRoot=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/zTPCProtonsTofCut.root

#SET THE OUTPUT IMAGE DIRECTORY (LEAVE BLANK TO SAVE NO IMAGES)
imageDir=$3
#imageDir=../userfiles/AuAu_4_5GeV_2015/analysisCleanup/tofCutImages/

#SET STARLIB
starlib=$4
######################################################

root -l -b -q ../macros/RunFitZTPCProtonsTofCut.C\(\"$yieldHistoFile\",\"$outRoot\",\"$starlib\",\"$imageDir\"\) #> /dev/null 2>&1 &

exit
