#!/bin/bash

######################################################

spectraFile=$1
#spectraFile=./protonSpectra_7_8.root

imageDir=$2

outFile=$3

energy=$4

starLibrary=$5

eventConfig=$6
#resultsFile=../userfiles/AuAu_4_5GeV_2015/analysis/proton_AlAu_THERMAL.root

#pid=2
#SET THE CENTRALITY BIN INDEX (USE -1 FOR ALL)
#charge=1

#energy=4.9

#eventConfig=FixedTarget2015

#midY=-1.5
#SET THE STARLIBRARY VERSION
#starlib=SL16a
######################################################

root -l -q -b   ../macros/RunFitThermalProtonSpectra.C\(\"$spectraFile\",\"$outFile\",\"$imageDir\",1\)
exit
