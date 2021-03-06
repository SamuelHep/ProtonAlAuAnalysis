#!/bin/bash

#This script is used to find the centrality binning of data
#It needs a Glauber Simulation to have been completed first.

#USER DEFINED VARIABLES
dataFileName=../AuAu62_temp.root           #File containing histogram to be matched
dataHistoName=refMult_0_20_Corrected_RefMultConverted           #Name of the histogram in the file from above
glauberFileName=../data/Glauber_197_197_35.5mb_WoodsSaxon.root #Glauber file generated from runGlauberSimulation.sh
outputFileName=../data/outfile1.root                            #outputfile (will get re-written if it already exists!)
normStartBinCenter=80                                      #The bin center value to begin the chi^2 matching/optimization routine
normStopBinCenter=-1                                        #The bin center value to end the chi^2 matching/optimization routine. Use -1 to use the last bin of the data histogram. NOTE: This value must be larger than normStartBinCenter to make any sense.

root -l -b -q ../macros/RunCentralityDetermination.C\(\"$dataFileName\",\"$dataHistoName\",\"$glauberFileName\",\"$outputFileName\",$normStartBinCenter,$normStopBinCenter\)
