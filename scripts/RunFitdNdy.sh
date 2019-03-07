#!/bin/bash


inputBlastWave=$1
outputFile=$2
imageDir=$3


root -l -b -q  ../macros/RunFitdNdy.C\(\"$inputBlastWave\",\"$outputFile\",\"$imageDir\"\)
