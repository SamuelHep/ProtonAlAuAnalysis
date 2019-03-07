//Extracts the yield of protons

#include <iostream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TVirtualFitter.h>

#include "globalDefinitions.h"
#include "ParticleInfo.h"
#include "TrackInfo.h"
#include "PrimaryVertexInfo.h"
#include "EventInfo.h"

#include "utilityFunctions.h"
#include "fitZTPCUtilities.h"


//-------FUNCTIONS DEFINED BELOW MAIN ------
void FitProtonMinus();
void FitProtonPlus();
void DrawPredictionLines(TCanvas *canvas);
void ComputeBinKinematics(Int_t mTm0Index, Double_t y, Int_t PID);
Double_t GetYFromTGraph(TGraphErrors* graph, Int_t mtm0Index);
Bool_t AddPointToSpectra(TH1D *yieldHisto, TF1 *yieldFit, TGraphErrors *rawSpectra,TGraphErrors *deutRawSpectra,TGraphErrors *deuteronRatioGraph, Int_t yIndex,
			 Int_t mTm0Index, Int_t PID, Int_t CHARGE,bool badfit, TString IMAGEDIR="");

//Canvases
TCanvas *fittingCanvas, *fitParCanvas, *spectraCanvas,*deutSpectraCanvas, *deuteronRatioCanvas;

//Particle Parameterizations and predictions
ParticleInfo *particleInfo;
ParticlePars *pion, *kaon, *proton, *electron, *deuteron;
TLine *pionLine, *kaonLine,*protonLine, *electronLine, *deuteronLine;
Double_t rapidity, mTm0, mT, pZ, pT, pTotal;
Double_t pionPredict, kaonPredict, protonPredict, electronPredict, deuteronPredict;
Double_t lineHeight;

//Pion Ratios
std::vector<Double_t> pionRatios;

//ZTOF Variable Histograms
std::vector<std::vector<TH1D *> > tofProtonMinusHisto(nRapidityBins,std::vector<TH1D *>(nmTm0Bins,(TH1D *)NULL));
std::vector<std::vector<TH1D *> > tofProtonPlusHisto(nRapidityBins,std::vector<TH1D *>(nmTm0Bins,(TH1D *)NULL));

//Fit Functions
std::vector<std::vector<TF1 *> > tofProtonMinusFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));
std::vector<std::vector<TF1 *> > tofProtonPlusFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));

//Fit Functions for zTPC fitting
std::vector<std::vector<TF1 *> > tofProtonOnlyFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));
std::vector<std::vector<TF1 *> > tofDeuteronOnlyFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));
std::vector<std::vector<TF1 *> > tofKaonOnlyFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));

//Raw Spectra
std::vector<TGraphErrors *> rawSpectraTOFProPlus(nRapidityBins,(TGraphErrors *)NULL);
std::vector<TGraphErrors *> deutRawSpectraTOFProPlus(nRapidityBins,(TGraphErrors *)NULL);
std::vector<TGraphErrors *> deuteronSpectraRatio(nRapidityBins,(TGraphErrors *)NULL);

std::vector <TF1 * > deuteronRatioTOFFit (nRapidityBins, (TF1 *)NULL);

//Number of Events Histogram
TH1D *nEventsHisto;

//User Specified Indices
Int_t rapidityIndex;

bool badFit [nRapidityBins][nmTm0Bins] = {0};


//-------MAIN-------
void fitZTOFProtons(TString ZTOFFILE, TString SPECTRAFILE, TString STARLIB, TString ZTOFGAUSFILE, Int_t RAPIDITYINDEX=-1, TString IMAGEDIR=""){



  Bool_t SAVEOBJ = true;
  Bool_t SAVEIMG = true;

  //If the ImageDir argument is empty then turn off saving images                                        
  if (!IMAGEDIR.CompareTo("") && SAVEIMG == true)
    SAVEIMG = false;


  //If we are going to save images then create subdirectories                           
  if (SAVEIMG && !gSystem->OpenDirectory(IMAGEDIR+"/ProtonPlus/")){
    gSystem->mkdir(IMAGEDIR+"/ProtonPlus/");
  }
  if (SAVEIMG && !gSystem->OpenDirectory(IMAGEDIR+"ProtonMinus/")){
    gSystem->mkdir(IMAGEDIR+"/ProtonMinus/");
  }
  
  //Pass the Value of the user's input indices to the globals
  rapidityIndex = RAPIDITYINDEX;

  //Create the Particle Info Object which has all the Bichesl Curves
  particleInfo = new ParticleInfo(STARLIB,true);
  
  //Create Tof Functions File for ZTPC fitting 
  TFile *tofFuncFile = new TFile(ZTOFGAUSFILE,"RECREATE");
  
  //Create the Output file and its directory structure
  TFile *outFile = new TFile(SPECTRAFILE,"RECREATE");
  
  if (SAVEOBJ){
    outFile->mkdir(Form("RawSpectra_%s",particleInfo->GetParticleName(PROTON,1).Data()));
    outFile->mkdir("FitParameterizations");
  }
  
  //Create the Particle Parameterization Objects
  pion = new ParticlePars();
  kaon = new ParticlePars();
  proton = new ParticlePars();
  electron = new ParticlePars();
  deuteron = new ParticlePars();
  
  DefineParticlePars(pion,"Pion");
  DefineParticlePars(kaon,"Kaon");
  DefineParticlePars(proton,"Proton");
  DefineParticlePars(electron,"Electron");
  DefineParticlePars(deuteron,"Deuteron");

  //Particle ZTOF Prediction Lines
  pionLine = new TLine();
  kaonLine = new TLine();
  protonLine = new TLine();
  electronLine = new TLine();
  deuteronLine = new TLine();

  pionLine->SetLineColor(particleInfo->GetParticleColor(PION));
  kaonLine->SetLineColor(particleInfo->GetParticleColor(KAON));
  protonLine->SetLineColor(particleInfo->GetParticleColor(PROTON));
  electronLine->SetLineColor(particleInfo->GetParticleColor(ELECTRON));
  deuteronLine->SetLineColor(particleInfo->GetParticleColor(DEUTERON));

  pionLine->SetLineWidth(2);
  kaonLine->SetLineWidth(2);
  protonLine->SetLineWidth(2);
  electronLine->SetLineWidth(2);
  deuteronLine->SetLineWidth(2);

  //Get the Yield Histogram file and Get the 3D Z Histograms
  //Because the ZTOF is currently being used to create a binned nSigma variable, only the All centrality histogram is loaded
  
  TFile *inFile = new TFile(ZTOFFILE,"READ");
  TH3D *protonPlusZTOF3D = (TH3D *)inFile->Get("ProtonPlus/zTOF_ProtonPlus_AllCent_tofCut");
  TH3D *protonMinusZTOF3D = (TH3D *)inFile->Get("ProtonMinus/zTOF_ProtonMinus_AllCent_tofCut");
  LoadYieldHistograms(protonPlusZTOF3D,&tofProtonPlusHisto,"TOF",rapidityIndex);
  LoadYieldHistograms(protonMinusZTOF3D,&tofProtonMinusHisto,"TOF",rapidityIndex);
  nEventsHisto = (TH1D *)inFile->Get("nEvents");
  
  //Set Style Options
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(112);
  
  //Create Global Canvases
  fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas",20,20,1000,1000);
  fittingCanvas->SetLogy();
  
  fitParCanvas = new TCanvas("fitParCanvas","Fit Parameter Canvas",20,20,1000,1000);
  fitParCanvas->Divide(3,4);
  spectraCanvas = new TCanvas("spectraCanvas","Spectra Canvas",20,850,800,600);
  spectraCanvas->SetLogy();
  deutSpectraCanvas = new TCanvas("deut_spectraCanvas","Deuteron Spectra Canvas",20,850,800,600);
  deutSpectraCanvas->SetLogy();
  deuteronRatioCanvas = new TCanvas("deuteronRatioCanvas","Ratio Canvas",20,850,800,600);
  
  
  //Fit the Proton Minuses
  //  FitProtonMinus();
  
  //Fit the Proton Pluses
  FitProtonPlus();

  //Make the Proton Plus Spectra
  for (Int_t yIndex=0; yIndex<nRapidityBins; yIndex++){
    
    //If the user has specified a rapidityIndex then only do that index
    if (rapidityIndex >= 0 && yIndex != rapidityIndex)
      continue;
    
    for (Int_t mTm0Index=12; mTm0Index<nmTm0Bins; mTm0Index++){
      
      tofProtonOnlyFit.at(yIndex).at(mTm0Index) = NULL;
      tofDeuteronOnlyFit.at(yIndex).at(mTm0Index) = NULL;
      tofKaonOnlyFit.at(yIndex).at(mTm0Index) = NULL;
      
      //If the spectra doesn't exist already then make it
      if (!rawSpectraTOFProPlus.at(yIndex)){
	
	rawSpectraTOFProPlus.at(yIndex) = new TGraphErrors();
	rawSpectraTOFProPlus.at(yIndex)->SetMarkerStyle(20);
	rawSpectraTOFProPlus.at(yIndex)->SetName(Form("rawSpectraTOF_%s_yIndex%02d",
						      particleInfo->GetParticleName(PROTON,1).Data(),						
						      yIndex));
	
	deutRawSpectraTOFProPlus.at(yIndex) = new TGraphErrors();
	deutRawSpectraTOFProPlus.at(yIndex)->SetMarkerStyle(20);
	deutRawSpectraTOFProPlus.at(yIndex)->SetName(Form("deut_rawSpectraTOF_%s_yIndex%02d",
							   particleInfo->GetParticleName(PROTON,1).Data(),						
							   yIndex));
	
	deuteronSpectraRatio.at(yIndex) = new TGraphErrors();
	deuteronSpectraRatio.at(yIndex)->SetMarkerStyle(20);
	deuteronSpectraRatio.at(yIndex)->SetName(Form("deuteronSpectraRatioTOF_%s_yIndex%02d",
						      particleInfo->GetParticleName(PROTON,1).Data(),
						      yIndex));
	
	deuteronRatioTOFFit.at(yIndex) = new TF1(Form("deuteronRatioFit_%02d",yIndex),"[0]*exp([1]+[2]*x)",0.4,2);
	
	
      }
      
      Bool_t success = AddPointToSpectra(tofProtonPlusHisto.at(yIndex).at(mTm0Index),
					 tofProtonPlusFit.at(yIndex).at(mTm0Index),
					 rawSpectraTOFProPlus.at(yIndex),deutRawSpectraTOFProPlus.at(yIndex),deuteronSpectraRatio.at(yIndex),
					 yIndex,mTm0Index,PROTON,1,badFit[yIndex][mTm0Index],IMAGEDIR);
      
      if (!success)
	continue;
      

      //Save the Gaussian Fits for use in the zTPC fitting       
      if( tofProtonOnlyFit.at(yIndex).at(mTm0Index) ){
	
	tofFuncFile->cd();
	tofProtonOnlyFit.at(yIndex).at(mTm0Index)->Write();
	
      }

      if( tofDeuteronOnlyFit.at(yIndex).at(mTm0Index) ){

        tofFuncFile->cd();
        tofDeuteronOnlyFit.at(yIndex).at(mTm0Index)->Write();

      }
      
      if( tofKaonOnlyFit.at(yIndex).at(mTm0Index) ){
	
        tofFuncFile->cd();
        tofKaonOnlyFit.at(yIndex).at(mTm0Index)->Write();
	
      }
      
    }//End loop over mTm0Index
    
    //Save the Spectra
    if (rawSpectraTOFProPlus.at(yIndex)->GetN() > 3 && SAVEOBJ){
      outFile->cd();
      outFile->cd(Form("RawSpectra_%s",particleInfo->GetParticleName(PROTON,1).Data()));
      rawSpectraTOFProPlus.at(yIndex)->Write();
      
      double fitMin = TMath::MinElement(deuteronSpectraRatio.at(yIndex)->GetN(), deuteronSpectraRatio.at(yIndex)->GetX() );
      double fitMax = TMath::MaxElement(deuteronSpectraRatio.at(yIndex)->GetN(), deuteronSpectraRatio.at(yIndex)->GetX() );

      deuteronRatioTOFFit.at(yIndex)->SetRange(fitMin,fitMax);
      deuteronSpectraRatio.at(yIndex)->Fit(deuteronRatioTOFFit.at(yIndex),"RQM");
      deuteronRatioTOFFit.at(yIndex)->Write();

      deutRawSpectraTOFProPlus.at(yIndex)->Write();
      deuteronSpectraRatio.at(yIndex)->Write();
    }

  }//End loop over yIndex
}//END OF FNXN

void FitProtonMinus(){


  //Local Variables
  TH1D *yieldHisto = NULL;
  TF1 *yieldFit   = NULL;

  for (Int_t yIndex=0; yIndex<nRapidityBins; yIndex++){

    //If the user has specified a rapidityIndex then only do that index
    if (rapidityIndex >= 0 && yIndex != rapidityIndex)
      continue;

    //Compute the Rapidity
    rapidity = GetRapidityRangeCenter(yIndex);

    cout <<"Fit ProtonMinus --- Working on Rapidity Index: " <<yIndex <<" Rapidity: " <<rapidity <<endl;

    //----- ROUND 1 - 
    for (Int_t mTm0Index=0; mTm0Index<nmTm0Bins; mTm0Index++){

      //Make sure the histogram exists
      if (!tofProtonMinusHisto.at(yIndex).at(mTm0Index))
  continue;

      //Pass the histogram to a local pointer
      yieldHisto = tofProtonMinusHisto.at(yIndex).at(mTm0Index);
      yieldHisto->SetAxisRange(-0.5,0.5);

      //Make sure it has enough entires
      if (yieldHisto->GetEntries() < 100) 
	continue;
      
      //Compute the kinematics for this bin
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      
      //Make the Predictions for the particle species
      pionPredict = particleInfo->PredictZTOF(pTotal,PROTON,PION);
      kaonPredict = particleInfo->PredictZTOF(pTotal,PROTON,KAON);
      protonPredict = particleInfo->PredictZTOF(pTotal,PROTON,PROTON);
      electronPredict = particleInfo->PredictZTOF(pTotal,PROTON,ELECTRON);
      
      
      //Stop at this mTm0
      if (mTm0 > 1.0)
	break;
      
      //Draw the Yield Histogram
      yieldHisto->SetAxisRange(-0.5,0.5);
      yieldHisto->SetMarkerStyle(kOpenCircle);
      fittingCanvas->cd();
      yieldHisto->Draw("E");
      
      Double_t piAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict));
      Double_t KAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(kaonPredict));
      
      //Create the Fit and pass it to a local pointer
      if (tofProtonMinusFit.at(yIndex).at(mTm0Index))
	tofProtonMinusFit.at(yIndex).at(mTm0Index)->Delete();
      
      Bool_t kaonOverlap = false;
      if (pionPredict-0.015 < kaonPredict && kaonPredict < pionPredict+0.015) kaonOverlap=true;
      if (KAmpSeed < 5 || kaonOverlap) 
	tofProtonMinusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)",-0.5,0.);
      else tofProtonMinusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)+gaus(3)",-0.5,0.);
      
      yieldFit = tofProtonMinusFit.at(yIndex).at(mTm0Index);
      
      
      //Set the Parameters of the Fit Function
      yieldFit->SetParameter(0, yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict)));
      yieldFit->SetParameter(1, particleInfo->PredictZTOF(pTotal, PROTON, PION));
      yieldFit->SetParameter(2, 0.005);
      yieldFit->SetParLimits(2, 0,.1);
      
      
      //Set Kaon/Electron parameters
      if (KAmpSeed >=5 && !kaonOverlap) {
	yieldFit->SetParameter(3, yieldHisto->GetBinContent(yieldHisto->FindBin(kaonPredict)));
	yieldFit->SetParameter(4, particleInfo->PredictZTOF(pTotal, PROTON, KAON));
	yieldFit->SetParameter(5, 0.005);
	yieldFit->SetParLimits(5, 0,0.1);
	yieldFit->SetParLimits(4,kaonPredict-0.02,kaonPredict+0.02);
	yieldFit->SetParLimits(5,0.001,0.1);
      }
      
      //Fit the Histogram
      yieldHisto->Fit(yieldFit,"RQ");
      
      //Draw the Prediction Lines
      lineHeight = yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict));
      DrawPredictionLines(fittingCanvas);
      fittingCanvas->Update();
      
      //Add the Mean and Width to their respective graphs
      AddPointToParGraph(pion->amp[yIndex],yieldFit,0,mTm0Index);
      AddPointToParGraph(pion->mean[yIndex],yieldFit,1,mTm0Index);
      AddPointToParGraph(pion->width[yIndex],yieldFit,2,mTm0Index);
      fitParCanvas->cd(1);
      pion->amp[yIndex]->Draw("APZ");
      fitParCanvas->cd(2);
      pion->mean[yIndex]->Draw("APZ");
      fitParCanvas->cd(3);
      pion->width[yIndex]->Draw("APZ");
      fitParCanvas->Update();
      
      if (KAmpSeed >=5 && !kaonOverlap) {
	AddPointToParGraph(kaon->amp[yIndex],yieldFit,3,mTm0Index);
	AddPointToParGraph(kaon->mean[yIndex],yieldFit,4,mTm0Index);
	AddPointToParGraph(kaon->width[yIndex],yieldFit,5,mTm0Index);
	fitParCanvas->cd(4);
	kaon->amp[yIndex]->Draw("APZ");
	fitParCanvas->cd(5);
	kaon->mean[yIndex]->Draw("APZ");
	fitParCanvas->cd(6);
	kaon->width[yIndex]->Draw("APZ");
	fitParCanvas->Update();
      }
      
    }//End of mtm0 loop
  }//End of rapidity loop
}//END FIT PI MINUS

//--------------------------------------------------------------------
void FitProtonPlus(){
  
  //Local Variables
  TH1D *yieldHisto = NULL;
  TF1 *yieldFit = NULL;
  
  for (Int_t yIndex=0; yIndex<nRapidityBins; yIndex++){
    
    //If the user has specified a rapidityIndex then only do that index
    if (rapidityIndex >= 0 && yIndex != rapidityIndex)
      continue;
    
    //Compute the Rapidity
    rapidity = GetRapidityRangeCenter(yIndex);
    
    cout <<"Fit ProtonPlus --- Working on Rapidity Index: " <<yIndex <<" Rapidity: " <<rapidity <<endl;
    
    //-----ROUND 1 
    for (Int_t mTm0Index=12; mTm0Index<nmTm0Bins; mTm0Index++){
      
      //Make sure the histogram exists
      if (!tofProtonPlusHisto.at(yIndex).at(mTm0Index))
	continue;
      
      //Pass the histogram to a local pointer
      yieldHisto = tofProtonPlusHisto.at(yIndex).at(mTm0Index);
      yieldHisto->SetAxisRange(-0.5,0.5);
      
      //Make sure it has enough entries
      if (yieldHisto->GetEntries() < 100)
	continue;
      
      //Compute the kinematics for this bin
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      
      //Make the Predictions for the particle species
      pionPredict = particleInfo->PredictZTOF(pTotal,PROTON,PION);
      kaonPredict = particleInfo->PredictZTOF(pTotal,PROTON,KAON);
      protonPredict = particleInfo->PredictZTOF(pTotal,PROTON,PROTON);
      electronPredict = particleInfo->PredictZTOF(pTotal,PROTON,ELECTRON);
      deuteronPredict = particleInfo->PredictZTOF(pTotal,PROTON,DEUTERON);
      
      //Stop at this mTm0
      if (mTm0 > 1.0)
	break;
      
      //Draw the Histogram
      yieldHisto->SetMarkerStyle(kOpenCircle);
      fittingCanvas->cd();
      yieldHisto->Draw("E");
      
      Double_t KAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(kaonPredict));
      Double_t proAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(protonPredict));
      Double_t deuAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(deuteronPredict));
      Double_t piPlusAmp     = GetYFromTGraph(pion->amp[yIndex], mTm0Index)*1.0;
      Double_t pionMeanPrev  = pionPredict;
      Double_t pionWidthPrev = .005;
      Double_t kaonMeanPrev  = kaonPredict;
      Double_t kaonWidthPrev = .005;
      
      if (proAmpSeed < 5) continue;
      
      //Create the Fit and pass it to a local pointer
      if (tofProtonPlusFit.at(yIndex).at(mTm0Index))
	tofProtonPlusFit.at(yIndex).at(mTm0Index)->Delete();
      
      tofProtonPlusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)+gaus(3)+gaus(6)+gaus(9)",-0.5,0.5);
      yieldFit = tofProtonPlusFit.at(yIndex).at(mTm0Index);
      
      //Set Parameter Names
      yieldFit->SetParNames(Form("%s_{Amp.}",particleInfo->GetParticleSymbol(PION).Data()),
			    "#mu_{#pi}","#sigma_{#pi}",
			    "p_{Amp.}","#mu_{p}","#sigma_{p}","K_{Amp.}","#mu_{K}","#sigma_{K}");
      
      //Set the Fit Parameters
      yieldFit->SetParameter(0,piPlusAmp);
      yieldFit->SetParameter(1,pionMeanPrev);
      yieldFit->SetParameter(2,pionWidthPrev);
      yieldFit->SetParLimits(2, 0,.01);
      
      yieldFit->SetParameter(3,proAmpSeed);
      yieldFit->SetParameter(4,protonPredict);
      yieldFit->SetParameter(5,0.01);
      yieldFit->SetParLimits(5, 0,0.1);
      yieldFit->SetParLimits(3,proAmpSeed*0.8,proAmpSeed*1.5);
      
      yieldFit->SetParameter(9,deuAmpSeed);
      yieldFit->SetParameter(10,deuteronPredict);
      yieldFit->SetParameter(11,0.01);
      yieldFit->SetParLimits(11,0.0001,0.2);
      yieldFit->SetParLimits(9,deuAmpSeed*0.8,deuAmpSeed*1.5);
      
      yieldFit->SetParameter(6,KAmpSeed);
      yieldFit->FixParameter(7,kaonMeanPrev);
      yieldFit->FixParameter(8,kaonWidthPrev);
      yieldFit->SetParLimits(6,KAmpSeed*0.8,KAmpSeed*1.5);
      
      //Fit the Histogram
      yieldHisto->Fit(yieldFit,"R");
      
      //Draw the Prediction Lines
      lineHeight = yieldHisto->GetBinContent(yieldHisto->FindBin(protonPredict));
      DrawPredictionLines(fittingCanvas);
      fittingCanvas->Update();
      
      AddPointToParGraph(deuteron->amp[yIndex],yieldFit,9,mTm0Index);
      AddPointToParGraph(deuteron->mean[yIndex],yieldFit,10,mTm0Index);
      AddPointToParGraph(deuteron->width[yIndex],yieldFit,11,mTm0Index);
      
      AddPointToParGraph(proton->amp[yIndex],yieldFit,3,mTm0Index);
      AddPointToParGraph(proton->mean[yIndex],yieldFit,4,mTm0Index);
      AddPointToParGraph(proton->width[yIndex],yieldFit,5,mTm0Index);
      
      //Check for empty bins
      badFit[yIndex][mTm0Index] = false;
           
      if (true){
	fitParCanvas->cd(10);
	gPad->SetLogy();
	deuteron->amp[yIndex]->Draw("APZ");
	fitParCanvas->cd(11);
	deuteron->mean[yIndex]->Draw("APZ");
	fitParCanvas->cd(12);
	deuteron->width[yIndex]->Draw("APZ");
	fitParCanvas->Update();
	
	fitParCanvas->cd(7);
	gPad->SetLogy();
	proton->amp[yIndex]->Draw("APZ");
	fitParCanvas->cd(8);
	proton->mean[yIndex]->Draw("APZ");
	fitParCanvas->cd(9);
	proton->width[yIndex]->Draw("APZ");
	fitParCanvas->Update();
	
      }
      
    }//End loop over mTm0 Index
    
  }//End loop over yIndex
  
}//END FIT PI PLUS

//--------------------------------------------------------------------
void DrawPredictionLines(TCanvas *canvas){
  
  canvas->cd();
  
  //Make prediction lines 
  PredictionLine(pionLine,pionPredict,lineHeight)->Draw("SAME");
  PredictionLine(kaonLine,kaonPredict,lineHeight)->Draw("SAME");
  PredictionLine(protonLine,protonPredict,lineHeight)->Draw("SAME");
  PredictionLine(electronLine,electronPredict,lineHeight)->Draw("SAME");
  PredictionLine(deuteronLine,deuteronPredict,lineHeight)->Draw("SAME");

}

//--------------------------------------------------------------------
void ComputeBinKinematics(Int_t mTm0Index, Double_t y, Int_t PID){

  mTm0 = GetmTm0RangeCenter(mTm0Index);
  mT   = mTm0 + particleInfo->GetParticleMass(PID);
  pZ   = ComputepZ(mT,y);
  pT   = ConvertmTm0ToPt(mTm0,particleInfo->GetParticleMass(PID));
  pTotal = ComputepTotal(pT,pZ);

}

//--------------------------------------------------------------------
Double_t GetYFromTGraph(TGraphErrors* graph, Int_t mtm0Index){

  Bool_t badmTm0Index = true;

  Int_t nPoints = graph->GetN();
  if (!graph || nPoints <= 0){ 
    cout<<"WARNING: GetYFromTGraph -- TGraph has no points"<<endl;
    return 0.0; 
  }
  Double_t *xPointsArr = graph->GetX();
  Double_t *yPointsArr = graph->GetY();
  for (Int_t iPoint = 0; iPoint<nPoints; iPoint++){

    Int_t xIndex = GetmTm0Index(xPointsArr[iPoint]);
    if (xIndex == mtm0Index){
      badmTm0Index == false;
      return yPointsArr[iPoint];
    }        

    else if (iPoint == nPoints-1 && badmTm0Index == true){
      cout<<"WARNING: Amplitude point for mTm0 Index does not exist"<<endl;
      return 0.0; 
    }        
  }

}

//--------------------------------------------------------------------
Bool_t AddPointToSpectra(TH1D *yieldHisto, TF1 *yieldFit, TGraphErrors *rawSpectra,TGraphErrors *deutRawSpectra,TGraphErrors *deuteronRatioGraph,
			 Int_t yIndex, Int_t mTm0Index, Int_t PID, Int_t CHARGE,bool badfit, TString IMAGEDIR){

  //Adds a point to the spectra
  //returns false on fail, true on success

  Bool_t SAVEIMG = true;

  //Make sure both the histogram and fit exist
  if (!yieldHisto || !yieldFit){ 
    return false;
	}

  if(badfit) return false;

  rapidity = GetRapidityRangeCenter(yIndex);
  ComputeBinKinematics(mTm0Index,rapidity,PROTON);

  //Make the Predictions for the particle species
  pionPredict = particleInfo->PredictZTOF(pTotal,PROTON,PION);
  kaonPredict = particleInfo->PredictZTOF(pTotal,PROTON,KAON);
  protonPredict = particleInfo->PredictZTOF(pTotal,PROTON,PROTON);
  electronPredict = particleInfo->PredictZTOF(pTotal,PROTON,ELECTRON);
  deuteronPredict = particleInfo->PredictZTOF(pTotal,PROTON,DEUTERON);

  //Draw the Final Version of the Fit
  fittingCanvas->Clear();
  fittingCanvas->cd();
  fittingCanvas->SetTopMargin(.05);
  fittingCanvas->SetRightMargin(.05);
  yieldHisto->SetMarkerStyle(kOpenCircle);
  yieldHisto->SetTitle("");
  yieldHisto->SetAxisRange(-0.5,0.5);
  yieldHisto->Draw("E");
  yieldHisto->SetStats(0);

  //Draw the Fit
  yieldFit->SetLineWidth(3);
  yieldFit->SetLineColor(kBlack);
  yieldFit->SetNpx(10000);
  yieldFit->Draw("SAME");

  //Create a Label
  TPaveText label(.12,.8,.35,.92,"BL NDC");
  label.SetTextAlign(11);
  label.AddText(Form("Z(%s) Yield Extraction Fit",particleInfo->GetParticleSymbol(PID,CHARGE).Data()));
  label.AddText(Form("  %.2f #leq y_{%s} < %.2f | %.2f #leq m_{T}-m_{%s} < %.2f",
                     GetRapidityRangeLow(yIndex),particleInfo->GetParticleSymbol(PID,CHARGE).Data(),
                     GetRapidityRangeHigh(yIndex),GetmTm0RangeLow(mTm0Index),
                     particleInfo->GetParticleSymbol(PID,CHARGE).Data(),GetmTm0RangeHigh(mTm0Index)));
  label.GetLine(0)->SetTextSize(.03);
  label.GetLine(1)->SetTextSize(.02);
  label.SetFillColor(kWhite);
  label.SetBorderSize(0);
  label.Draw("SAME");

  //Make A Function for Each of the Particles
  TF1 pionGaus("pion","gaus(0)",-20,20);
  pionGaus.FixParameter(0,yieldFit->GetParameter(0));
  pionGaus.FixParameter(1,yieldFit->GetParameter(1));
  pionGaus.FixParameter(2,yieldFit->GetParameter(2));
  pionGaus.SetParError(0,yieldFit->GetParError(0));
  pionGaus.SetParError(1,yieldFit->GetParError(1));
  pionGaus.SetParError(2,yieldFit->GetParError(2));
  pionGaus.SetLineColor(particleInfo->GetParticleColor(PION));
  pionGaus.SetLineWidth(3);
  pionGaus.SetLineStyle(7);
  pionGaus.SetNpx(10000);
  pionGaus.Draw("SAME");

  TF1 kaonGaus("kaon","gaus(0)",-20,20);
  kaonGaus.FixParameter(0,yieldFit->GetParameter(6));
  kaonGaus.FixParameter(1,yieldFit->GetParameter(7));
  kaonGaus.FixParameter(2,yieldFit->GetParameter(8));
  kaonGaus.SetParError(0,yieldFit->GetParError(6));
  kaonGaus.SetParError(1,yieldFit->GetParError(7));
  kaonGaus.SetParError(2,yieldFit->GetParError(8));
  kaonGaus.SetLineColor(particleInfo->GetParticleColor(KAON));
  kaonGaus.SetLineWidth(3);
  kaonGaus.SetLineStyle(7);
  kaonGaus.SetNpx(10000);
  kaonGaus.Draw("SAME");
    
  TF1 protonGaus("proton","gaus(0)",-20,20);
  protonGaus.FixParameter(0,yieldFit->GetParameter(3));
  protonGaus.FixParameter(1,yieldFit->GetParameter(4));
  protonGaus.FixParameter(2,yieldFit->GetParameter(5));
  protonGaus.SetParError(0,yieldFit->GetParError(3));
  protonGaus.SetParError(1,yieldFit->GetParError(4));
  protonGaus.SetParError(2,yieldFit->GetParError(5));
  protonGaus.SetLineColor(particleInfo->GetParticleColor(PROTON));
  protonGaus.SetLineWidth(3);
  protonGaus.SetLineStyle(7);
  protonGaus.SetNpx(10000);
  protonGaus.Draw("SAME");

  TF1 deuteronGaus("deuteron","gaus(0)",-20,20);
  deuteronGaus.FixParameter(0,yieldFit->GetParameter(9));
  deuteronGaus.FixParameter(1,yieldFit->GetParameter(10));
  deuteronGaus.FixParameter(2,yieldFit->GetParameter(11));
  deuteronGaus.SetParError(0,yieldFit->GetParError(9));
  deuteronGaus.SetParError(1,yieldFit->GetParError(10));
  deuteronGaus.SetParError(2,yieldFit->GetParError(11));
  deuteronGaus.SetLineColor(particleInfo->GetParticleColor(DEUTERON));
  deuteronGaus.SetLineWidth(3);
  deuteronGaus.SetLineStyle(7);
  deuteronGaus.SetNpx(10000);
  deuteronGaus.Draw("SAME");

  
  //Make single Gaussian Functions for zTPC fitting
  tofProtonOnlyFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)",-.5,.5);
  tofProtonOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(0,protonGaus.GetParameter(0));
  tofProtonOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(1,protonGaus.GetParameter(1));
  tofProtonOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(2,protonGaus.GetParameter(2));
  tofProtonOnlyFit.at(yIndex).at(mTm0Index)->SetName(Form("ProtonPlusTofFit_yIndex%02d_mTm0Index%02d",yIndex,mTm0Index));

  tofDeuteronOnlyFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)",-.5,.5);
  tofDeuteronOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(0,deuteronGaus.GetParameter(0));
  tofDeuteronOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(1,deuteronGaus.GetParameter(1));
  tofDeuteronOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(2,deuteronGaus.GetParameter(2));
  tofDeuteronOnlyFit.at(yIndex).at(mTm0Index)->SetName(Form("DeuteronPlusTofFit_yIndex%02d_mTm0Index%02d",yIndex,mTm0Index));

  tofKaonOnlyFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)",-.5,.5);
  tofKaonOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(0,kaonGaus.GetParameter(0));
  tofKaonOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(1,kaonGaus.GetParameter(1));
  tofKaonOnlyFit.at(yIndex).at(mTm0Index)->FixParameter(2,kaonGaus.GetParameter(2));
  tofKaonOnlyFit.at(yIndex).at(mTm0Index)->SetName(Form("KaonPlusTofFit_yIndex%02d_mTm0Index%02d",yIndex,mTm0Index));

  fittingCanvas->Update();

  gPad->Update();

  //Make A Legend
  TLegend leg;
  leg.SetOption("brNDC");

  leg.SetX1NDC(.65);
  leg.SetY1NDC(.89);
  leg.SetX2NDC(.92);
  leg.SetY2NDC(.94);

  leg.SetTextSize(.035);
  leg.SetNColumns(5);
  leg.AddEntry(yieldFit,"Fit","L");
  leg.AddEntry("pion","#pi","L");
  leg.AddEntry("kaon","K","L");
  leg.AddEntry("proton","p","L");
  leg.SetFillColor(kWhite);
  leg.SetBorderSize(0);
  leg.Draw("SAME");

  DrawPredictionLines(fittingCanvas);
  fittingCanvas->Update();

  //Extract the Yield
  Double_t amp = protonGaus.GetParameter(0);
  Double_t ampError = protonGaus.GetParError(0);
  Double_t width = protonGaus.GetParameter(2);
  Double_t widthError = protonGaus.GetParError(2);
  Double_t pi = 3.14159;

  Double_t deut_amp = deuteronGaus.GetParameter(0);
  Double_t deut_ampError = deuteronGaus.GetParError(0);
  Double_t deut_width = deuteronGaus.GetParameter(2);
  Double_t deut_widthError = deuteronGaus.GetParError(2);

  Double_t zBinWidth = (Double_t)yieldHisto->GetBinWidth(2);
  Double_t mTm0BinCenter = GetmTm0RangeCenter(mTm0Index);
  
  Double_t nEvents = 0;

  //Sum All Events
  for (int iCentBin=0; iCentBin<=5;iCentBin++)
    {
      nEvents = nEvents + (Double_t)nEventsHisto->GetBinContent(nEventsHisto->FindBin(iCentBin));
    }

  //Set the Normalization Factor
  Double_t normFactor = (1.0/(mTm0BinCenter+particleInfo->GetParticleMass(PID))) * (1.0/mTm0BinWidth) * (1.0/rapidityBinWidth) *
    (1.0/nEvents) * (1.0/(2.0*pi));

  //Compute the Raw Yield
  Double_t rawYield, deut_rawYield, rawYieldError, deut_rawYieldError;;
  rawYield = ( (amp * width * sqrt(2.0*pi)) / zBinWidth ) * normFactor;
  rawYieldError = ( (amp * width * sqrt(pow(ampError/amp,2)+pow(widthError/width,2)) *
		     sqrt(2.0*pi)) / zBinWidth) * normFactor;
  
  
  //Compute the Raw Yield for Deuteron
  deut_rawYield = ( (deut_amp * deut_width * sqrt(2.0*pi)) / zBinWidth ) * normFactor;
  deut_rawYieldError = ( (deut_amp * width * sqrt(pow(deut_ampError/deut_amp,2)+pow(widthError/width,2)) *
			  sqrt(2.0*pi)) / zBinWidth) * normFactor;
  
  double deut_fracErr = deut_rawYieldError/deut_rawYield;
  double fracErr = rawYieldError/rawYield;
  
  
  //Calculate ratio of deuteron to proton yield
  if (rawYield > 0)
    { 
      Double_t ratioDeuteron = deut_rawYield/rawYield;
      Double_t totDeuteron_err = sqrt(deut_fracErr*deut_fracErr + fracErr*fracErr);  
      
      //Add point to deuteron to proton yield graph
      deuteronRatioGraph->SetPoint(deuteronRatioGraph->GetN(),GetmTm0RangeCenter(mTm0Index),ratioDeuteron);
      deuteronRatioGraph->SetPointError(deuteronRatioGraph->GetN()-1,0,totDeuteron_err);
      
    }
  
  //Add Point to Spectrum
  rawSpectra->SetPoint(rawSpectra->GetN(),GetmTm0RangeCenter(mTm0Index),rawYield);
  rawSpectra->SetPointError(rawSpectra->GetN()-1,0,rawYieldError);
  
  deutRawSpectra->SetPoint(deutRawSpectra->GetN(),GetmTm0RangeCenter(mTm0Index),deut_rawYield);
  deutRawSpectra->SetPointError(deutRawSpectra->GetN()-1,0,deut_rawYieldError);
  
  
  //Draw the Spectra
  spectraCanvas->cd();
  rawSpectra->Draw("APZ");
  spectraCanvas->Update();
  
  deutSpectraCanvas->cd();
  deutRawSpectra->Draw("APZ");
  deutSpectraCanvas->Update();
  
  //draw ratio graph
  deuteronRatioCanvas->cd();
  deuteronRatioGraph->Draw("APZ");
  deuteronRatioCanvas->Update();
  
  //Save the Fit Image
  TString fullImageDir = IMAGEDIR+"/"+particleInfo->GetParticleName(PID,CHARGE);
  TString imageName = TString::Format("zTOFFit_%s_%02d_%02d",
				      particleInfo->GetParticleName(PROTON,CHARGE).Data(),yIndex,mTm0Index);
  
  fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
  
  return true;
  
}
