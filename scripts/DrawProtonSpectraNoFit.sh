#!/bin/sh
spectraFile=$1
eventConfig=$2
systemLabel=$3
energy=$4
imageDir=$5

root -b -l -q  ../drawingmacros/DrawProtonSpectraNoFit.C\(\"$spectraFile\",\"$imageDir\",\"$eventConfig\",\"$systemLabel\",$energy,0\)
