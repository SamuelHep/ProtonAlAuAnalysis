#!/bin/bash

########################################################
#                    USER VARIABLES                    #
########################################################
#Change for user enviroment

#DATA DIRECTORY
dataDir=/scratch_menkar/FixedTargetData/AlAu_4_9GeV_2015/

#IMAGE DIRECTORY, MUST EXIST BEFORE RUNNING
imageDir=../userfiles/AlAu_4_9GeV_2015/images_thermal/

#ANALYSIS OUT DIRECTORY
analysisDir=../userfiles/AlAu_4_9GeV_2015/analysis_thermal/

#CORRECTION FILE
correctionFile=../userfiles/CorCurves_2_14.root

#SET NUMBER OF EVENTS (PER DATA FILE)
nEvents=-1

#CM ENERGY
energy=4.9

#SET THE EVENT CONFIG (Only tested for FixedTarget2015)
eventConfig=FixedTarget2015

#LIBRARY CONFIG (Only tested on SL16a)
starlib=SL16a

#CENTRALITY INDEX
centIndex=-999

#PLOT LABEL
systemLabel=Al+Au

########################################################
#                     OUTPUTFILE NAMES                 #
########################################################
#Not necessary to change

#fit function for spectra
tofFitFile=zTOFFit.root

#Binned histogram of data 
yieldHistoFile=YieldHistograms.root

#TPC spectra file
tpcSpectraFile=zTPCProtonSpectra.root

#Clean proton TPC spectra. Used to parameterize the means and widths of protons 
protonWidthFile=zTPCProtonTofCut.root

#Clean deuteron TPC spectra. Used to parameterize the means and widths deuterons
deuteronWidthFile=zTPCDeuteronTofCut.root

#TOF spectra file
tofSpectraFile=zTOFSpectra.root

#dNdy Fit
dndyFile=protondNdyFit.root

#dNdy Fit File
dndyFitFile=protondNdyGausFit.root

#TOF Cut label
tofCutLabel=YieldHistograms_tofCut.root

########################################################

########################################################
echo ""
echo ":::::Files:::::"

yieldHistoFile=$analysisDir$yieldHistoFile
echo "   Yield Histogram file: $yieldHistoFile"

yieldHistoTofCutFile=$analysisDir$tofCutLabel
echo "   Yield Histogram file: $yieldHistoFile"

tpcSpectraFile=$analysisDir$tpcSpectraFile
echo "   TPC Spectra file: $tpcSpectraFile"

tofSpectraFile=$analysisDir$tofSpectraFile
echo "   TOF Spectra file: $tofSpectraFile"
 
tofFitFile=$analysisDir$tofFitFile
echo "   TOF fit file: $tofFitFile"

protonWidthFile=$analysisDir$protonWidthFile
echo "   protonWidth file: $protonWidthFile"

deuteronWidthFile=$analysisDir$deuteronWidthFile
echo "   deuteronWidth file: $deuteronWidthFile"

dndyFile=$analysisDir$dndyFile
echo "   dndy file: $dndyFile"

dndyFitFile=$analysisDir$dndyFitFile
echo "   dndy gaus fit file: $dndyFitFile"

########################################################
echo ""
echo ":::::Image Directories:::::"

tofImageDir=tofImages
tofImageDir=$imageDir$tofImageDir
echo "   TOF image Dir: $tofImageDir"

tpcImageDir=tpcImages
tpcImageDir=$imageDir$tpcImageDir
echo "   TPC image Dir: $tpcImageDir"

protonWidthImageDir=protonWidthImages
protonWidthImageDir=$imageDir$protonWidthImageDir
echo "   ProWidth image Dir: $protonWidthImageDir"

deuteronWidthImageDir=deuteronWidthImages
deuteronWidthImageDir=$imageDir$deuteronWidthImageDir
echo "   DeuWidth image Dir: $deuteronWidthImageDir"

#Make Image Sub-Directories#
mkdir $tofImageDir
mkdir $tpcImageDir
mkdir $protonWidthImageDir
mkdir $deuteronWidthImageDir


########################################################

echo ""

#./RunSkimmerAndBinner.bash $dataDir $analysisDir $nEvents $starlib $energy $eventConfig $tofFitFile 0

echo "YieldHisto file: $yieldHistoFile"
echo "Spectra file: $tofSpectraFile"
echo "tofImageDir file: $tofImageDir"
echo "tofFit file: $tofFitFile"
echo "Star Library: $starlib"
echo ""

#./RunFitZTOFProtons.bash $yieldHistoFile $tofSpectraFile $tofImageDir $tofFitFile $starlib

#./RunSkimmerAndBinner.bash $dataDir $analysisDir $nEvents $starlib $energy $eventConfig $tofFitFile 1

./RunFitZTPCProtonsTofCut.sh $yieldHistoTofCutFile $protonWidthFile $protonWidthImageDir $starlib

#./RunFitZTPCDeuteronsTofCut.sh $yieldHistoTofCutFile $deuteronWidthFile $deuteronWidthImageDir $starlib

#./RunFitZTPCProtons.sh $yieldHistoFile $tpcSpectraFile $protonWidthFile $deuteronWidthFile $tofSpectraFile $tpcImageDir $centIndex $starlib

#./RunCorrectSpectra.sh $tpcSpectraFile $correctionFile $energy $eventConfig $starlib

#./RunFitThermalProtonSpectra.sh $tpcSpectraFile $imageDir $dndyFile $energy $starlib $eventConfig

#./DrawProtonSpectraNoFit.sh $tpcSpectraFile $eventConfig $systemLabel $energy $imageDir

#./RunFitdNdy.sh $dndyFile $dndyFitFile $imageDir
