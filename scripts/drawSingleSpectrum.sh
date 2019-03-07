#!/bin/sh

spectraFile=../userfiles/AuAu_4_5GeV_2015/analysis/protondNdyFit.root
imageDir=../userfiles/AuAu_4_5GeV_2015/images/spectraFits/

root -b -l ../drawingmacros/DrawAllSpectra.C\(\"$spectraFile\",\"$imageDir\",\"FixedTarget2015\",4.5,2,1\)

