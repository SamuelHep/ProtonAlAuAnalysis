#!/bin/bash

#This script runs the Glauber Simulation

#USER DEFINED VARIABLE QUANTITIES
nEvents=100000           #Number of DESIRED events with at least one participating nucleon
nNucleonsA=197           #Number of nucleons in nucleus A
nNucleonsB=197           #Number of nucleons in nucleus B
nnCrossSection=35.5      #The nucleon-nucleon inelastic crosssection in mb
nucleonDistribution=1    #The method for distributing the nucleons in the nuclei
                         #0=Uniform-Hard Sphere, 1=Woods-Saxon
outDirectory=../data/       #path to output directory (the file name is automatically generated)

#Run the Glauber Model With the Above Options
root -l -b -q ../macros/RunGlauberSimulation.C\($nEvents\,$nNucleonsA\,$nNucleonsB\,$nnCrossSection\,$nucleonDistribution,\"$outDirectory\"\)

