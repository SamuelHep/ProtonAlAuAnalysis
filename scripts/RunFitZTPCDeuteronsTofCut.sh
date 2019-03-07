#!/bin/bash

#This calls the RunFitZTPCDeuteronsTofCut.C file in the macros directory which loads the necessary libraries


######################################################
#SET THE YIELD HISTOGRAM FILE NAME AND PATH (FULL PATH)
yieldHistoFile=$1

#SET THE OUTFILE FILE NAME AND PATH (FULL PATH)
outRoot=$2

#SET THE OUTPUT IMAGE DIRECTORY (LEAVE BLANK TO SAVE NO IMAGES)
imageDir=$3

#SET THE STARLIBRARY
starlib=$4
######################################################

root -l -b -q ../macros/RunFitZTPCDeuteronsTofCut.C\(\"$yieldHistoFile\",\"$outRoot\",\"$starlib\",\"$imageDir\"\) #> /dev/null 2>&1 &

exit
