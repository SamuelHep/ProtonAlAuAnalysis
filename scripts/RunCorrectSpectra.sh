#!/bin/bash

######################################################


spectraFile=$1
correctionFile=$2
#correctionFile=../userfiles/correctionFile.root
pid=2
#SET THE CENTRALITY BIN INDEX (USE -1 FOR ALL)
charge=1

energy=$3

eventConfig=$4

#SET THE STARLIBRARY VERSION
starlib=$5
######################################################

root -l -q -b ../macros/RunCorrectSpectra.C\(\"$spectraFile\",\"$correctionFile\",$pid,$charge,\"$starlib\",$energy,\"$eventConfig\"\)
exit
