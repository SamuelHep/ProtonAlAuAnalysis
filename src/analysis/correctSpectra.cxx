#include <iostream>
#include <cstdio>
#include <vector>
#include <utility>
#include <algorithm>

#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TF3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TThread.h>

#include "globalDefinitions.h"
#include "utilityFunctions.h"
#include "StRefMultExtendedCorr.h"
#include "UserCuts.h"
#include "ParticleInfo.h"
#include "SpectraCorrectionCurves.h"
#include "SpectraCorrectionFunctions.h"
#include "SpectraClass.h"

bool draw = false;

//Structurs and functions for the multi-threaded systematic error algorithm
struct SystematicKernelArgs{
  TGraphErrors *rawSpectrum;
  TGraphErrors *correctedSpectrum;
  CorrectionCurves *corrCurves;
  double minTPCEfficiency;
};

void *SystematicKernel(void *ARGS){

  SystematicKernelArgs *args = (SystematicKernelArgs *)ARGS;
  TGraphErrors *rawSpectrum = args->rawSpectrum;
  CorrectionCurves *corrCurves = args->corrCurves;
  double minTPCEfficiency = args->minTPCEfficiency;
  
  args->correctedSpectrum = ApplyCorrectionChain(rawSpectrum,corrCurves,minTPCEfficiency,true);
   
}

std::pair<std::vector<TThread *>,std::vector<SystematicKernelArgs *> >
CreateThreads(TGraphErrors *rawSpectrum, CorrectionCurves *corrCurve,
	      double minTPCEfficiency,const int nThreads){
  std::vector<TThread *> threads(nThreads);
  std::vector<SystematicKernelArgs *> threadArgs(nThreads);
  for (int iThread=0; iThread<nThreads; iThread++){

    threadArgs.at(iThread) = new SystematicKernelArgs();
    threadArgs.at(iThread)->rawSpectrum = rawSpectrum;
    threadArgs.at(iThread)->correctedSpectrum = NULL;
    threadArgs.at(iThread)->corrCurves = corrCurve;
    threadArgs.at(iThread)->minTPCEfficiency = minTPCEfficiency;

    threads.at(iThread) = new TThread(Form("thread_%d",iThread),
				      SystematicKernel,(void *)threadArgs.at(iThread));

  }  
    
  return std::make_pair(threads,threadArgs);
}

void RunThreads(std::vector<TThread *> *threads, const int nThreads, const int nConcurrentThreads){

  int finishedThreads(0);
  std::vector<TThread *> runningThread;
  
  while (finishedThreads < nThreads){

    for (int iThread=0; iThread<nConcurrentThreads; iThread++){

      runningThread.push_back(threads->at(finishedThreads));
      runningThread.back()->Run();
      finishedThreads++;      
      
    }//End Loop Over Submit Thread

    //Join the current batch
    for (unsigned int iThread=0; iThread<runningThread.size(); iThread++){
      runningThread.at(iThread)->Join();
      runningThread.at(iThread)->Delete();
    }
    runningThread.clear();
    
  }//End Finish Threads
  

}

//____MAIN_______________________________________________________________
void CorrectSpectra(TString spectraFileName, TString correctionFileName, Int_t pid, Int_t charge){

  ParticleInfo *particleInfo = new ParticleInfo();
  const int nCentBins = GetNCentralityBins();
  const int pidIndex = pid;
  const int pidCharge = charge;
  Double_t minTPCEfficiency=0.3;
  Double_t pidMass = particleInfo->GetParticleMass(pid);
  
  //Open the SpectraFile
  TFile *spectraFile = new TFile(spectraFileName,"UPDATE");
  spectraFile->cd();
  spectraFile->mkdir(Form("CorrectedSpectra_%s",particleInfo->GetParticleName(pid,charge).Data()));

  //Open the Correction File
  TFile *corrFile = new TFile(correctionFileName,"READ");

  //Spectra Pointers
  TGraphErrors *rawSpectrum[2] = {NULL,NULL};
  TGraphErrors *fullyCorrectedSpectrum[2] = {NULL,NULL};
  TGraphErrors *systematicErrorSpectrum[2] = {NULL,NULL};
  TGraphErrors *totalErrorSpectrum[2] = {NULL,NULL};
  TGraphErrors *widthSpectrum[2] =  {NULL,NULL};


  //For the Systematic Error Computation
  const int nSpectra(1000);
  const int nConcurrentThreads(100);
  std::vector<TGraphErrors *> tpcSpectraEnsemble(nSpectra);
  std::vector<TGraphErrors *> tofSpectraEnsemble(nSpectra);

  //Loop Over the Centrality Bins and Rapidity Bins and apply the Corrections
  for (int iCentBin=0; iCentBin<nCentBins; iCentBin++){
    for (int yIndex=GetMinRapidityIndexOfInterest(pid); yIndex<=GetMaxRapidityIndexOfInterest(pid); yIndex++){
      
      cout <<  "iCentBin: " <<  iCentBin << " yIndex: " << yIndex << endl;
      
      //------------------------------------------------------
      //----Get the Correction Curves     --------------------
      //----See SpectraCorrectionCurves.cxx   ----------------
      //------------------------------------------------------      
      
      CorrectionCurves corrCurve(corrFile,pid,charge,iCentBin,yIndex);
      
      //Make Sure we have everything we need
      bool isGoodToGo = true;
      //  if (!corrCurve.GetInterBinMeanCurve()){
      //	cout <<Form("INFO - correctSpectra.cxx - iCentBin: %02d yIndex: %02d InterBinMeanCurve does not exist. SKIPPING!\n",iCentBin,yIndex);
      //	isGoodToGo = false;
      //      }
      //      if (pid == PION && !corrCurve.GetMuonBackgroundCurvePowerLaw()){
      //	cout <<Form("INFO - correctSpectra.cxx - iCentBin: %02d yIndex: %02d MuonBackgroundCurve does not exist. SKIPPING!\n",iCentBin,yIndex);
      //	isGoodToGo = false;
      //      }
      if (pid != KAON && !corrCurve.GetFeedDownBackgroundCurveExpo()){
   	cout <<Form("INFO - correctSpectra.cxx - iCentBin: %02d yIndex: %02d FeedDownCurve does not exist. SKIPPING!\n",iCentBin,yIndex);
    	isGoodToGo = false;
      }
      if (!corrCurve.GetEnergyLossCurve()){
	cout <<Form("INFO - correctSpectra.cxx - iCentBin: %02d yIndex: %02d EnergyLossCurve does not exist. SKIPPING!\n",iCentBin,yIndex);
	isGoodToGo = false;
      }
      if (!corrCurve.GetEfficiencyCurve()){
	cout <<Form("INFO - correctSpectra.cxx - iCentBin: %02d yIndex: %02d TPCEfficiencyCurve does not exist. SKIPPING!\n",iCentBin,yIndex);
	isGoodToGo = false;
      }
      
      if (!isGoodToGo){
	continue;
      }
      
      TString detector[2] = {"","TOF"};
      //Two Loops one for the TPC and the Other for the TOF
      for (int iSpectrum=0; iSpectrum<2; iSpectrum++){

	
	TString rawSpectraName = TString::Format("rawSpectra%s_%s_Cent%02d_yIndex%02d",
						 detector[iSpectrum].Data(),
						 particleInfo->GetParticleName(pid,charge).Data(),
						 iCentBin,yIndex);
	TString corrSpectraName = TString::Format("correctedSpectra%s_%s_Cent%02d_yIndex%02d",
						  detector[iSpectrum].Data(),
						  particleInfo->GetParticleName(pid,charge).Data(),
						  iCentBin,yIndex);
	
	cout <<"Correcting Spectrum: " <<rawSpectraName.Data() <<endl;
	
	
	//Get the Spectrum to be Corrected
	spectraFile->cd();
	spectraFile->cd(Form("RawSpectra_%s",particleInfo->GetParticleName(pid,charge).Data()));
	rawSpectrum[iSpectrum] = (TGraphErrors *)gDirectory->Get(Form("%s",rawSpectraName.Data()));
	
	//Skip if the Spectrum Does not exist
	if (!rawSpectrum[iSpectrum]){
	  cout <<"    Spectrum Does not Exist. SKIPPING!" <<endl;
	  delete rawSpectrum[iSpectrum];
	  rawSpectrum[iSpectrum] = NULL;
	  continue;
	}
	
	//Skip if the Spectrum has no points
	if (rawSpectrum[iSpectrum]->GetN() == 0){
	  cout <<"    Spectrum does not have any points. SKIPPING" <<endl;
	  delete rawSpectrum[iSpectrum];
	  rawSpectrum[iSpectrum] = NULL;
	  continue;
	}
	
	//------------------------------------------------------
	//----Apply the Nominal Corrections --------------------
	//----See SpectraCorrectionFunctions.cxx----------------
	//------------------------------------------------------
	cout <<"     Nominal Corrections...." <<endl;
	
	fullyCorrectedSpectrum[iSpectrum] = ApplyCorrectionChain(rawSpectrum[iSpectrum],&corrCurve,minTPCEfficiency,false);

	cout << "N points : " << rawSpectrum[iSpectrum]->GetN()  <<endl;
	
	//Set Names and Titles
	fullyCorrectedSpectrum[iSpectrum]->SetName(corrSpectraName.Data());
	fullyCorrectedSpectrum[iSpectrum]->
	  SetTitle(Form("%s Spectrum %s Cent=[%d,%d]%% y_{%s}=[%.02f,%.02f];m_{T}-m_{%s};",
			iSpectrum == 0 ? "TPC" : "TOF",
			particleInfo->GetParticleSymbol(pidIndex,charge).Data(),
			//	iCentBin!=nCentBins-1 ? (int)GetCentralityPercents().at(iCentBin+1):0,
			0,
			(int)GetCentralityPercents().at(iCentBin),
			particleInfo->GetParticleSymbol(pidIndex).Data(),
			GetRapidityRangeLow(yIndex),GetRapidityRangeHigh(yIndex),
			particleInfo->GetParticleSymbol(pidIndex).Data()));
	fullyCorrectedSpectrum[iSpectrum]->SetMarkerStyle(particleInfo->GetParticleMarker(pid,charge));
	fullyCorrectedSpectrum[iSpectrum]->SetMarkerColor(particleInfo->GetParticleColor(pid));
	
	//------------------------------------------------------
	//----Compute Systematic Errors ------------------------
	//------------------------------------------------------
	//Create an ensemble of n corrected spectra where each spectrum
	//has been corrected in systematic error mode.
	
	
	cout <<"     Systematic Corrections...." <<endl;
	
	if (nSpectra < nConcurrentThreads){
	  cout <<"ERROR - CorrectSpectra.cxx -- The number of spectra used in the systematic error\n calculation must be larger than the number of concurrent threads. EXITING!\n";
	  exit (EXIT_FAILURE);
	  
	}
	
	//Create the Threads
	std::pair<std::vector<TThread *>,std::vector<SystematicKernelArgs *> > threads
	  = CreateThreads(rawSpectrum[iSpectrum],&corrCurve,minTPCEfficiency,nSpectra);
	
	//Run The Threads
	RunThreads(&threads.first,threads.first.size(),nConcurrentThreads);

	std::vector<TGraphErrors *> *spectraEnsemble = NULL;
	if (iSpectrum == 0)
	  spectraEnsemble = &tpcSpectraEnsemble;
	else
	  spectraEnsemble = &tofSpectraEnsemble;
	for (unsigned int iSpectra=0; iSpectra<nSpectra; iSpectra++){
	  spectraEnsemble->at(iSpectra) = threads.second.at(iSpectra)->correctedSpectrum;
	}
	//Create the Systematic Error Spectrum
	systematicErrorSpectrum[iSpectrum] =
	MakeSystematicErrorSpectrum(fullyCorrectedSpectrum[iSpectrum],nSpectra,spectraEnsemble);

	//Create the Total Error Spectrum
	totalErrorSpectrum[iSpectrum] =
	  MakeTotalErrorSpectrum(fullyCorrectedSpectrum[iSpectrum],systematicErrorSpectrum[iSpectrum]);

	//Set Names so they have StatErr at the end before saving
	fullyCorrectedSpectrum[iSpectrum]->SetName(Form("%s_StatErr",corrSpectraName.Data()));

	//Set the Names so that the Spectrum with the total error does not have any extension
	totalErrorSpectrum[iSpectrum]->SetName(Form("%s",corrSpectraName.Data()));
	totalErrorSpectrum[iSpectrum]->SetTitle(fullyCorrectedSpectrum[iSpectrum]->GetTitle());

	//Clean Up Some Vectors
	for (int i=0; i<nSpectra; i++){
	  if (threads.second.at(i))
	    delete threads.second.at(i);
	}
	
      }
      //-----------------------------------------------
      //End Loop Over TPC and TOF Spectra
      //-----------------------------------------------

      //Check to Make sure at least one of either the TPC or TOF Spectra Exist
      if (!fullyCorrectedSpectrum[0] && !fullyCorrectedSpectrum[1]){
	cout <<"Neither the TPC nor the TOF Spectrum Exists for this centrality and Rapidity Bin. SKIPPING Rapidity Bin!" <<endl;
	continue;
      }

      if (fullyCorrectedSpectrum[0]->GetN() <= 0 ) continue;
      
      //Clean the Spectra
      cout <<"Cleaning Spectra by removing outling points." <<endl;
      std::vector<int> removedPointsTPC;
      std::vector<int> removedPointsTOF;
      
      cout <<"Combining TPC and TOF Spectra into SpectraClass" <<endl;
      
      //Store the Spectra in the SpectraClass Object
      SpectraClass spectraTotal(pid,charge,iCentBin,yIndex,
				totalErrorSpectrum[0],fullyCorrectedSpectrum[1]);
      
      spectraTotal.SetSpectrumName(Form("CorrectedSpectra_%s_Cent%02d_yIndex%02d_Total",
					particleInfo->GetParticleName(pid,charge).Data(),
					iCentBin,yIndex));

      spectraTotal.SetSpectraFile(spectraFile);

      spectraTotal.SetSpectraDir(Form("CorrectedSpectra_%s",
				      particleInfo->GetParticleName(pid,charge).Data()));

      spectraTotal.SetSpectraClassDir(Form("SpectraClass_%s",
					   particleInfo->GetParticleName(pid,charge).Data()));

      //Require that the full spectrum have some minimum number of points
      if (spectraTotal.GetN() > 4){
	
	//Save - Overwrite any spectrum or spectra class object that exists
	spectraTotal.WriteSpectrum(true,true);
	cout <<"Spectrum Saved" <<endl;
      
	
	//OPTIONAL DRAW
	if (draw && iCentBin==8 && yIndex == 20 && pid == PION && charge == -1 ){
	  TH1F *yieldSpreadHisto = new TH1F("yieldSpreadHisto","yieldSpreadHisto",200,0,30);
	  TH1F *mTm0SpreadHisto  = new TH1F("mTm0SpreadHisto","mTm0SpreadHisto",100,0,.0002);
	  TCanvas *ensembleCanvas = new TCanvas("ensembleCanvas","ensembleCanvas",20,20,800,600);
	  TH1F *frame = ensembleCanvas->DrawFrame(0,0.1,2,500);
	  for (int iSpectra=0; iSpectra<nSpectra; iSpectra++){
	    tpcSpectraEnsemble.at(iSpectra)->SetMarkerStyle(7);
	    tofSpectraEnsemble.at(iSpectra)->SetMarkerStyle(7);
	    tpcSpectraEnsemble.at(iSpectra)->Draw("PZ");
	    tofSpectraEnsemble.at(iSpectra)->Draw("PZ");
	    int point=TGraphFindPoint(tpcSpectraEnsemble.at(iSpectra),.1125);
	    double yield = tpcSpectraEnsemble.at(iSpectra)->GetY()[point];
	    double mtm0  = tpcSpectraEnsemble.at(iSpectra)->GetX()[point];
	    yieldSpreadHisto->Fill(fabs(yield - fullyCorrectedSpectrum[0]->GetY()[point]));
	    mTm0SpreadHisto->Fill(fabs(mtm0 - fullyCorrectedSpectrum[0]->GetX()[point]));
	  }
	  fullyCorrectedSpectrum[0]->Draw("PZ");
	  fullyCorrectedSpectrum[1]->Draw("PZ");
	  ensembleCanvas->SaveAs(Form("ensembleCanvas_%s_Cent%02d_yIndex%02d.C",
				      particleInfo->GetParticleName(pid,charge).Data(),
				      iCentBin,yIndex));
	  TCanvas *canvas1 = new TCanvas("canvas1","canvas1",20,20,800,600);
	  yieldSpreadHisto->Draw();
	  canvas1->SaveAs("YieldSpreadHisto.C");
	  mTm0SpreadHisto->Draw();
	  canvas1->SaveAs("mTm0SpreadHisto.C");
	  
	}

      }//End Require Minimum number of points
	
      //Do Some Clean Up
      for (int i=0; i<2; i++){
	if (rawSpectrum[i]){
	  delete rawSpectrum[i];
	  rawSpectrum[i] = NULL;
	}
	if (fullyCorrectedSpectrum[i]){
	  delete fullyCorrectedSpectrum[i];
	  fullyCorrectedSpectrum[i] = NULL;
	}
	if (systematicErrorSpectrum[i]){
	  delete systematicErrorSpectrum[i];
	  systematicErrorSpectrum[i] = NULL;
	}
	if (totalErrorSpectrum[i]){
	  delete totalErrorSpectrum[i];
	  totalErrorSpectrum[i] = NULL;
	}
      }
      
      tpcSpectraEnsemble.clear();
      tofSpectraEnsemble.clear();

      tpcSpectraEnsemble.resize(nSpectra);
      tofSpectraEnsemble.resize(nSpectra);
      
    }//End Loop Over Rapidity Bins

  }//End Loop Over Centrality Bins
  

}//END MAIN


