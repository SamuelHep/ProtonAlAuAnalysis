#!/bin/bash

#This class the RunFitZTOFPions.C file in the macros directory which loads the necessary
#libraries and then runs fitZTOFProtons.cxx in the src/analysis directory.

######################################################
#SET THE YIELD HISTOGRAM FILE NAME AND PATH (FULL PATH)
yieldHistoFile=$1

#SET THE SPECTRA FILE NAME AND PATH (FULL PATH)
spectraFile=$2

#SET THE OUTPUT IMAGE DIRECTORY (LEAVE BLANK TO SAVE NO IMAGES)
imageDir=$3

#SET GAUS FILE
gausFile=$4

#SET THE RAPIDITY VALUE (USE -999 FOR ALL)rapidityValue=-0.72
rapidityValue=-999

#SET THE STARLIBRARY VERSION
starlib=$5
######################################################

root -l -b -q ../macros/RunFitZTOFProtons.C\(\"$yieldHistoFile\",\"$spectraFile\",\"$starlib\",\"$gausFile\",$rapidityValue,\"$imageDir\"\) #> /dev/null 2>&1 &

exit
