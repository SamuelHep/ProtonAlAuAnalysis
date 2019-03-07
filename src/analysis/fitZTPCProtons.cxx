//Extracts the yield of protons
//Currently assumes amplitude of pi+ yield to be the same as the amplitude of the pi- yield
//Currently does not account for deuteron contamination

#include <iostream>
#include <vector>
#include <algorithm>

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
#include <THistPainter.h>

#include "globalDefinitions.h"
#include "ParticleInfo.h"
#include "TrackInfo.h"
#include "PrimaryVertexInfo.h"
#include "EventInfo.h"
#include "utilityFunctions.h"
#include "fitZTPCUtilities.h"

//-------FUNCTIONS DEFINED BELOW MAIN ------
void FitProtonMinus(TString imageDirect);
void FitProtonPlus(TString imageDirect);
Bool_t DrawFits(TH1D *yieldHisto, TF1 *yieldFit, Int_t yIndex,			 
	      Int_t mTm0Index, Int_t PID, Int_t CHARGE, TString IMAGEDIR="");
void DrawPredictionLines(TCanvas *canvas, Int_t charge);
void ComputeBinKinematics(Int_t mTm0Index, Double_t y, Int_t PID);
Bool_t AddPointToSpectra(TH1D *yieldHisto, TF1 *yieldFit, TGraphErrors *rawSpectra,TGraphErrors * ratioTPCdeut,TGraphErrors * ratioTPCpion, Int_t yIndex,
			 Int_t mTm0Index, Int_t PID, Int_t CHARGE, TString IMAGEDIR="");
Double_t GetYFromTGraph(TGraphErrors* graph, Int_t mtm0Index);
Double_t GetYFromMeanTGraph(TGraphErrors* graph, Int_t mtm0Index);
void loadWidths(TString protonWidth, TString deuteronWidth);
Float_t LowStatAmpSeed(TH1D * yieldHisto);

//-------GLOBALS--------

//Save?
Bool_t SAVEOBJ = true;
Bool_t SAVEIMG = true;
Bool_t DRAWIMG = true;

//Canvases
TCanvas *fittingCanvas, *fitParCanvas, *spectraCanvas,*ratioCanvas, *confCanvas;

//Particle Parameterizations and predictions
ParticleInfo *particleInfo;
ParticlePars *pion, *kaon, *proton, *electron, *deuteron;
TLine *pionLine, *kaonLine,*protonLine, *electronLine, *deuteronLine;
Double_t rapidity, mTm0, mT, pZ, pT, pTotal;
Double_t pionPredict, kaonPredict, protonPredict, electronPredict, deuteronPredict;
Double_t lineHeight;

//ZTPC Variable Histograms
std::vector<std::vector<TH1D *> > tpcProtonMinusHisto(nRapidityBins,std::vector<TH1D *>(nmTm0Bins,(TH1D *)NULL));
std::vector<std::vector<TH1D *> > tpcProtonPlusHisto(nRapidityBins,std::vector<TH1D *>(nmTm0Bins,(TH1D *)NULL));

//Fit Functions
std::vector<std::vector<TF1 *> > tpcProtonMinusFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));
std::vector<std::vector<TF1 *> > tpcProtonPlusFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));

//Raw Spectra
std::vector<TGraphErrors *> rawSpectraProtonPlus(nRapidityBins,(TGraphErrors *)NULL);
std::vector<TGraphErrors *> ratioTPC_deut_prot(nRapidityBins,(TGraphErrors *)NULL);
std::vector<TGraphErrors *> ratioTPC_pion_prot(nRapidityBins,(TGraphErrors *)NULL);
std::vector<TF1 *> ratioFit(nRapidityBins,(TF1 *)NULL);

std::vector<TF1 *> protonWidthFit(nRapidityBins,(TF1 *)NULL);
std::vector<TGraphErrors *> protonMeanGraph(nRapidityBins,(TGraphErrors *)NULL);
std::vector<TF1 *> deuteronWidthFit(nRapidityBins,(TF1 *)NULL);
std::vector<TGraphErrors *> deuteronMeanGraph(nRapidityBins,(TGraphErrors *)NULL);
//NOTE: ABOVE IS TEMPORARILY COMMENTED OUT TO PREVENT COMPILING ERROR

//Number of Events Histogram
TH1D *nEventsHisto;

//Species Controls
Bool_t killElectron[nRapidityBins] = {0}, killElectronMean[nRapidityBins] = {0}, 
	killFirstProtonAttempt[nRapidityBins] = {0}, skipPiMinus[nRapidityBins] = {0};

//User Specified Indices
Int_t rapidityIndex;
Int_t centralityIndex;
Int_t firstRapidityIndex;
Int_t lastRapidityIndex;

Double_t mTm0Cutoff; //cutoff for electrons

//-------MAIN-------
void fitZTPCProtons(TString ZTPCFILE, TString SPECTRAFILE,TString PROTONWIDTH, TString DEUTERONWIDTH, TString AMPRATIOFILE, TString STARLIB, Int_t CENTRALITYINDEX,
		    Int_t RAPIDITYINDEX=-1, TString IMAGEDIR="",Int_t RAPIDITYINDEX_START = 0, Int_t RAPIDITYINDEX_END = 0){
  
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
  centralityIndex = CENTRALITYINDEX;
  
  cout << " ~~~~~~~~~~~~~~~~~~ " << "CENTRALITY INDEX:  " <<centralityIndex << "  ~~~~~~~~~~~~~~~~~~" << endl;
  
  //Options for running a few rapidity indexes
  if (RAPIDITYINDEX_END == 0 ) lastRapidityIndex = nRapidityBins;
  else lastRapidityIndex = RAPIDITYINDEX_END;
  firstRapidityIndex = RAPIDITYINDEX_START;
  
  //Create the Particle Info Object which has all the Bichsel Curves
  particleInfo = new ParticleInfo(STARLIB,true);
  
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
  
  //Particle ZTPC Prediction Lines
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
  
  pionLine->SetLineWidth(3);
  kaonLine->SetLineWidth(3);
  protonLine->SetLineWidth(3);
  electronLine->SetLineWidth(3);
  deuteronLine->SetLineWidth(3);
  
  //Load tofFile which contains the ratio of deuterons to protons in the TOF
  TFile * tofFile = new TFile(AMPRATIOFILE,"READ");
  tofFile->cd(); 

  //Load each ratioFit from the tofFile
  for (int iRapidityBin =7; iRapidityBin < 19; iRapidityBin++){

   ratioFit.at(iRapidityBin) = (TF1 *)tofFile->Get(Form("RawSpectra_ProtonPlus/deuteronRatioFit_%02d",iRapidityBin));
   TF1 * ratioTest = (TF1 *)tofFile->Get(Form("RawSpectra_ProtonPlus/deuteronRatioFit_%02d",iRapidityBin));
    
    if ( ratioFit.at(iRapidityBin) ){
    }
    else { 

      cout << "Quiting on: "<< Form("RawSpectra_ProtonPlus/deuteronRatioFit_%02d",iRapidityBin) << endl;
      cout << "TOF NOT FOUND... EXITING..." << endl;
      return;
    }
  }
  
  //Load the tofCut ZTPC proton and deuteron parameter graphs
  loadWidths(PROTONWIDTH,DEUTERONWIDTH);
  
  //Get the Yield Histogram file and Get the 3D Z Histograms 
  cout << "DEBUG:::0" << endl;

  TFile *inFile = new TFile(ZTPCFILE,"READ");
  TH3D *protonPlusZTPC3D = (TH3D *)inFile->Get(Form("ProtonPlus/zTPC_ProtonPlus_Cent%d",centralityIndex));
  TH3D *protonMinusZTPC3D = (TH3D *)inFile->Get(Form("ProtonMinus/zTPC_ProtonMinus_Cent%d",centralityIndex));
  LoadYieldHistograms(protonPlusZTPC3D,&tpcProtonPlusHisto,"TPC",rapidityIndex);
  LoadYieldHistograms(protonMinusZTPC3D,&tpcProtonMinusHisto,"TPC",rapidityIndex);
  nEventsHisto = (TH1D *)inFile->Get("nEvents");
  
  cout << "DEBUG:::1" << endl;

  //Set Style Options
  gStyle->SetOptFit(112);
  gStyle->SetTitleFontSize(0.1);
  
  cout << "DEBUG:::2" << endl;

  //Create Global Canvases
  fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas",20,20,800,600);
  fittingCanvas->SetLogy();
  fitParCanvas = new TCanvas("fitParCanvas","Fit Parameter Canvas",20,20,1100,1100);
  fitParCanvas->Divide(3,4);
  spectraCanvas = new TCanvas("spectraCanvas","Spectra Canvas",20,850,800,600);
  spectraCanvas->SetLogy();
  ratioCanvas = new TCanvas("ratioCanvas","Ratio Canvas",20,850,800,600);
  ratioCanvas->Divide(2,1);
  confCanvas = new TCanvas("confCanvas","Width Confidence Interval Canvas",20,20,800,600);
  confCanvas ->Divide(3,1);
  
  cout << "DEBUG:::3" << endl;

  //Fit the Proton Minuses (to parameterize the pions)
  FitProtonMinus(IMAGEDIR);
  
  //Fit the Proton Pluses 
  FitProtonPlus(IMAGEDIR);
  
  //Make the Proton Plus Spectra
  for (Int_t yIndex=firstRapidityIndex; yIndex<lastRapidityIndex; yIndex++){//nRapidityBins
    
    //If the user has specified a rapidityIndex then only do that index
    if (rapidityIndex >= 0 && yIndex != rapidityIndex)
      continue;
    
    //Kill if too few points to parameterize proton
    if (proton->mean[yIndex]->GetN() < 5) continue;
    
    
    //If the spectra doesn't exist already then make it
    if (!rawSpectraProtonPlus.at(yIndex)){
      rawSpectraProtonPlus.at(yIndex) = new TGraphErrors();
      rawSpectraProtonPlus.at(yIndex)->SetMarkerStyle(20);
      rawSpectraProtonPlus.at(yIndex)->SetName(Form("rawSpectra_%s_Cent%02d_yIndex%02d",
						    particleInfo->GetParticleName(PROTON,1).Data(),
						    centralityIndex,yIndex));
      ratioTPC_deut_prot.at(yIndex) = new TGraphErrors();
      ratioTPC_deut_prot.at(yIndex)->SetMarkerStyle(20);
      ratioTPC_deut_prot.at(yIndex)->SetName(Form("ratioTPC_%s_Cent%02d_yIndex%02d",
						  particleInfo->GetParticleName(DEUTERON,1).Data(),
						  centralityIndex,yIndex));
      ratioTPC_pion_prot.at(yIndex) = new TGraphErrors();
      ratioTPC_pion_prot.at(yIndex)->SetMarkerStyle(20);
      ratioTPC_pion_prot.at(yIndex)->SetName(Form("ratioTPC_%s_Cent%02d_yIndex%02d",
						  particleInfo->GetParticleName(PION,1).Data(),
						  centralityIndex,yIndex));
    


  }
    
    
    //Adds point to spectra
    for (Int_t mTm0Index=2; mTm0Index<nmTm0Bins; mTm0Index++){
      
      Bool_t success = AddPointToSpectra(tpcProtonPlusHisto.at(yIndex).at(mTm0Index),
					 tpcProtonPlusFit.at(yIndex).at(mTm0Index),rawSpectraProtonPlus.at(yIndex),ratioTPC_deut_prot.at(yIndex),ratioTPC_pion_prot.at(yIndex),
					 yIndex,mTm0Index,PROTON,1,IMAGEDIR);
            
      //Draws QA fits if fit success
      if (success) {
	DrawFits(tpcProtonMinusHisto.at(yIndex).at(mTm0Index),tpcProtonMinusFit.at(yIndex).at(mTm0Index),yIndex,mTm0Index,PROTON,-1,IMAGEDIR);
      }
    }//End Loop Over mTm0Index
    
    //Save the Spectra
    if (SAVEOBJ){
      outFile->cd();
      outFile->cd(Form("RawSpectra_%s",particleInfo->GetParticleName(PROTON,1).Data()));
      rawSpectraProtonPlus.at(yIndex)->Write();
      ratioTPC_deut_prot.at(yIndex)->Write();
      ratioTPC_pion_prot.at(yIndex)->Write();
    }    
    
    //Save the Parameterizations
    if (SAVEOBJ){
      outFile->cd();
      outFile->cd("FitParameterizations");
      pion->amp[yIndex]->Write();
      pion->mean[yIndex]->Write();
      pion->width[yIndex]->Write();
      pion->widthConf[yIndex]->Write();
      electron->amp[yIndex]->Write();
      electron->mean[yIndex]->Write();
      electron->width[yIndex]->Write();
      electron->widthConf[yIndex]->Write();
      proton->mean[yIndex]->Write();
      proton->width[yIndex]->Write();     
      proton->widthConf[yIndex]->Write();
    }    
  }//End Loop Over yIndex
}
//--------------------------------------------------------------------

void FitProtonMinus(TString imageDirect){

  //Make QA directory in /ProtonMinus/
  if (!gSystem->OpenDirectory(imageDirect+"/ProtonMinus/QA/")){
    gSystem->mkdir(imageDirect+"/ProtonMinus/QA/");
  }

  //Local Variables
  Int_t chargePID = -1;
  TH1D *yieldHisto = NULL;
  TF1 *yieldFit   = NULL;

  for (Int_t yIndex=firstRapidityIndex; yIndex<lastRapidityIndex; yIndex++){

    //Check if width parameter fit is properly loaded
    if (!protonWidthFit.at(yIndex)) continue;

    killElectron[yIndex] = false;
    killElectronMean[yIndex] = false;
    killFirstProtonAttempt[yIndex] = false;
    skipPiMinus[yIndex] = false;
    
    //If the user has specified a rapidityIndex then only do that index
    if (rapidityIndex >= 0 && yIndex != rapidityIndex)
      continue;
    
    //Compute the Rapidity
    rapidity = GetRapidityRangeCenter(yIndex);
    
    cout <<"Working on Rapidity Index: " <<yIndex <<" Rapidity: " <<rapidity <<endl;
    
    //----- ROUND 1 - PARAMETERIZE THE PION AND ELECTRON MEAN AND WIDTH
    for (Int_t mTm0Index=2; mTm0Index<nmTm0Bins; mTm0Index++){
      
      //Make sure the histogram exists
      if (!tpcProtonMinusHisto.at(yIndex).at(mTm0Index))
	continue;
      
      //Pass the histogram to a local pointer
      yieldHisto = tpcProtonMinusHisto.at(yIndex).at(mTm0Index);
      yieldHisto->SetAxisRange(-2.,2.5);
      
      //Make sure it has enough entries
      if (yieldHisto->GetEntries() < 100)
	continue;
      
      //Compute the kinematics for this bin
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      
      //Make the Predictions for the particle species
      pionPredict = particleInfo->PredictZTPC(pTotal,PROTON,PION);
      kaonPredict = particleInfo->PredictZTPC(pTotal,PROTON,KAON);
      protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);
      electronPredict = particleInfo->PredictZTPC(pTotal,PROTON,ELECTRON);
      deuteronPredict = particleInfo->PredictZTPC(pTotal,PROTON,DEUTERON);
      
      //Stop at this mTm0
      if (mTm0 > 1)
	break;
      
      //Determine if electron is underthreshold to be fit
      if (yieldHisto->GetBinContent(yieldHisto->FindBin(electronPredict)) < 5) killElectronMean[yIndex] = true;
      else killElectronMean[yIndex] = false;
      
      if (killElectronMean[yIndex] || (!killElectronMean[yIndex] && mTm0 > 0.5)){
	
	//Create the Fit and pass it to a local pointer
	if (!tpcProtonMinusFit.at(yIndex).at(mTm0Index)){
	  tpcProtonMinusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)",pionPredict-0.2,pionPredict+0.2);
	  yieldFit = tpcProtonMinusFit.at(yIndex).at(mTm0Index);
	}
	
	//Draw the Histogram
	fittingCanvas->cd();
	tpcProtonMinusHisto.at(yIndex).at(mTm0Index)->Draw("E");
	
	
	//Set the Parameters of the Fit Function
	yieldFit->SetParName(0,"#pi_{Amp}");
	yieldFit->SetParName(1,"#mu_{#pi}");
	yieldFit->SetParName(2,"{#sigma}_{#pi}");
	yieldFit->SetParameter(0,yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict)));
	yieldFit->SetParameter(1,pionPredict);
	yieldFit->SetParameter(2,.08);
	
	//Fit the Histogram
	yieldHisto->Fit(yieldFit,"RQ");
	
	//Draw the Prediction Lines
	lineHeight = yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict));
	if (DRAWIMG) DrawPredictionLines(fittingCanvas, chargePID);
	fittingCanvas->Update();
	
	//Save QA plots 
	if (SAVEIMG){
	  
	  gPad->Update();
	  
	  TString fullImageDir = imageDirect+"/ProtonMinus/QA/";
	  TString imageName = TString::Format("Fit_OnlyPionMinus_%s_%02d_%02d_%02d",
					      particleInfo->GetParticleName(2,chargePID).Data(),
					      centralityIndex,yIndex,mTm0Index);
	  fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	}
	
	//Add the Mean, Width and Amplitude to their respective graphs
	AddPointToParGraph(pion->amp[yIndex],yieldFit,0,mTm0Index);
	AddPointToParGraph(pion->mean[yIndex],yieldFit,1,mTm0Index);
	AddPointToParGraph(pion->width[yIndex],yieldFit,2,mTm0Index);
	if (DRAWIMG){
	  fitParCanvas->cd(1);
	  gPad->SetLogy();
	  pion->amp[yIndex]->Draw("APZ");
	  fitParCanvas->cd(2);
	  pion->mean[yIndex]->Draw("APZ");
	  fitParCanvas->cd(3);
	  pion->width[yIndex]->Draw("APZ");
	  fitParCanvas->Update();
	}
      }
      
      //Fit if electron is not killed
      if (!killElectronMean[yIndex]){
	
        if (mTm0 > 0.5) continue;
	
	//Create the Fit and pass it to a local pointer
        if (!tpcProtonMinusFit.at(yIndex).at(mTm0Index)){
	  tpcProtonMinusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0) + gaus(3)",pionPredict-0.3,electronPredict+0.3);
	  yieldFit = tpcProtonMinusFit.at(yIndex).at(mTm0Index);
        }
	
	yieldFit->SetParName(0,"#pi_{Amp}");
	yieldFit->SetParName(1,"#mu_{#pi}");
	yieldFit->SetParName(2,"{#sigma}_{#pi}");
	yieldFit->SetParName(3,"e_{Amp}");
	yieldFit->SetParName(4,"#mu_{e}");
	yieldFit->SetParName(5,"{#sigma}_{e}");

        //Draw the Histogram
        fittingCanvas->cd();
        tpcProtonMinusHisto.at(yIndex).at(mTm0Index)->Draw("E");
        	
        //Set the Parameters of the Fit Function
        yieldFit->SetParameter(0,yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict)));
        yieldFit->SetParameter(1,pionPredict);
        yieldFit->SetParameter(2,.08);
	
        yieldFit->SetParameter(3,yieldHisto->GetBinContent(yieldHisto->FindBin(electronPredict)));
	yieldFit->SetParameter(4,electronPredict);
	yieldFit->SetParameter(5,0.08);
	yieldFit->SetParLimits(3,yieldFit->GetParameter(3)*.8, yieldFit->GetParameter(3)*3);
	yieldFit->SetParLimits(5,.07,.12);
	
        //Fit the Histogram
        yieldHisto->Fit(yieldFit,"RQ");
	
	if (SAVEIMG){
	  
	  gPad->Update();
	  
	  TString fullImageDir = imageDirect+"/ProtonMinus/QA/";
	  TString imageName = TString::Format("Fit_PionElectronMinus_%s_%02d_%02d_%02d",
					      particleInfo->GetParticleName(2,chargePID).Data(),
					      centralityIndex,yIndex,mTm0Index);

	  fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));

	}


        //Draw the Prediction Lines
        lineHeight = yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict));
        if (DRAWIMG) DrawPredictionLines(fittingCanvas, chargePID);
        fittingCanvas->Update();
	
        //Add the Mean, Width and Amplitude to their respective graphs
        AddPointToParGraph(pion->amp[yIndex],yieldFit,0,mTm0Index);
        AddPointToParGraph(pion->mean[yIndex],yieldFit,1,mTm0Index);
        AddPointToParGraph(pion->width[yIndex],yieldFit,2,mTm0Index);
	AddPointToParGraph(electron->mean[yIndex],yieldFit,4,mTm0Index);
	AddPointToParGraph(electron->width[yIndex],yieldFit,5,mTm0Index);
        
	if (DRAWIMG){
	  fitParCanvas->cd(1);
	  pion->amp[yIndex]->Draw("APZ");
	  fitParCanvas->cd(2);
	  pion->mean[yIndex]->Draw("APZ");
	  fitParCanvas->cd(3);
	  pion->width[yIndex]->Draw("APZ");
	  fitParCanvas->cd(5);
	  electron->mean[yIndex]->Draw("APZ");
	  fitParCanvas->cd(6);
	  electron->width[yIndex]->Draw("APZ");
	  fitParCanvas->Update();
	}
      }
    }//End Loop Over mTm0Index
    
    //Get mtm0 cutoff value for electrons
    Int_t nPoints = electron->mean[yIndex]->GetN();
    Double_t x[1], y[1];
    electron->mean[yIndex]->GetPoint(nPoints-1, x[0], y[0]);
    mTm0Cutoff = ceil((x[0]+0.025)*100.)/100.;
    
    //Check if the pion parameters can be parameterized
    if (pion->mean[yIndex]->GetN() < 5){
      cout <<"fitZTPCProtons.cxx:FitProtonMinus() WARNING: No successful parameterizations for rapidity index: "
	   <<yIndex <<", rapidity = " <<rapidity <<endl;
      killElectron[yIndex] = true;
      skipPiMinus[yIndex] = true;
      continue;
    }
    
    //Fit the Pion Width as a function of mTm0
    pion->widthFit[yIndex] = new TF1(Form("pionMinus_WidthFit_%02d",yIndex),
				     SwooshFunc,.05,1.3,4);
    pion->widthFit[yIndex]->SetParameters(.2,.5,.3,.1);
    //pion->widthFit[yIndex]->SetParLimits(1,-100,0.);
    pion->widthFit[yIndex]->SetParLimits(3,.0,3);
    pion->widthFit[yIndex]->SetNpx(10000);
    
    //Get the Maximum mTm0 from the graph and use it as a bound on the parameter
    Double_t max = pion->width[yIndex]->GetX()[pion->width[yIndex]->GetN()-1]; //TGraphErrors->GetX() returns an array of X points
    pion->widthFit[yIndex]->SetParLimits(2,.1,max);
    
    fitParCanvas->cd(3);
    pion->width[yIndex]->Fit(pion->widthFit[yIndex],"RQ");
    
    fitParCanvas->Update();
    
    
    if (SAVEIMG){

      gPad->Update();
      
      TString fullImageDir = imageDirect+"/ProtonMinus/QA/";
      TString imageName = TString::Format("Fit_ParameterizedPionMinus_%s_%02d_%02d",
					  particleInfo->GetParticleName(2,chargePID).Data(),
					  centralityIndex,yIndex);
      
      fitParCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
    }
    
    //Determine Confidence Interval 
    Double_t minRange = pion->width[yIndex]->GetX()[0];;
    Double_t maxRange = max;
    
    pion->widthConf[yIndex]->Set(pion->widthFit[yIndex]->GetNpx());
    Double_t step = (maxRange-minRange)/(double)pion->widthConf[yIndex]->GetN();               
    for (int i=0; i<pion->widthConf[yIndex]->GetN(); i++){                                     
      pion->widthConf[yIndex]->SetPoint(i,minRange + i*step,0);                                
    }                                                                                          
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(pion->widthConf[yIndex],.95);//.68);
    pion->widthConf[yIndex]->SetFillColor(kOrange);                                            
    pion->widthConf[yIndex]->SetFillStyle(3001);
    
    //Graph Width Confidence Interval
    if (DRAWIMG){
      confCanvas->cd(1);
      pion->widthConf[yIndex]->Draw(); //A=>without axis. Not sure what Z is, P not applicable 
      confCanvas->Update();
    }
    
    cout <<"Finished Parameterizing the Pion Mean and Width!" <<endl;

    //Get the number of points in the electron mean graph and decide if it
    //should be killed or not
    if (electron->mean[yIndex]->GetN() < 6){
      killElectron[yIndex] = true;
      cout <<"INFO: fitZTPCProtons.cxx:FitProtonMinus(): The electron has been killed due to too few points in the mean parameterization in rapidity index: " <<yIndex <<" rapidity = " <<rapidity <<endl;
    }
    else killElectron[yIndex] = false;
    
    //If the electron is not killed then fit its mean and width
    if (!killElectron[yIndex]){
      
      //Electron Width
      electron->widthFit[yIndex] = new TF1(Form("electronMinus_WidthFit_%02d",yIndex),
					   SwooshFunc,.05,1.0,4);
      electron->widthFit[yIndex]->SetParameters(.2,.5,.3,.1);
      electron->widthFit[yIndex]->SetParLimits(3,.0,3);
      electron->widthFit[yIndex]->SetNpx(10000);
      
      //Get the Maximum mTm0 from the graph and use it as a bound on the parameter
      max = electron->width[yIndex]->GetX()[electron->width[yIndex]->GetN()-1];
      electron->widthFit[yIndex]->SetParLimits(2,.1,max);
      
      fitParCanvas->cd(6);
      electron->width[yIndex]->Fit(electron->widthFit[yIndex],"R");
      fitParCanvas->Update();
      
      if (SAVEIMG){
	
	gPad->Update();
	
	TString fullImageDir = imageDirect+"/ProtonMinus/QA/";
	TString imageName = TString::Format("parameterizedElectronMinus_%s_%02d_%02d",
					    particleInfo->GetParticleName(2,chargePID).Data(),
					    centralityIndex,yIndex);
	
	fitParCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
      }
      
      //Determine Confidence Interval
      minRange = electron->width[yIndex]->GetX()[0];
      maxRange = max;
      
      electron->widthConf[yIndex]->Set(electron->widthFit[yIndex]->GetNpx());
      step = (maxRange-minRange)/(double)electron->widthConf[yIndex]->GetN();
      
      
      for (int i=0; i<electron->widthConf[yIndex]->GetN(); i++){
        electron->widthConf[yIndex]->SetPoint(i,minRange + i*step,0);
      }
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(electron->widthConf[yIndex],.68);
      electron->widthConf[yIndex]->SetFillColor(kOrange);
      electron->widthConf[yIndex]->SetFillStyle(3001);
      
      //Graph Width Confidence Interval
      if (DRAWIMG) {
	confCanvas->cd(2);
	electron->widthConf[yIndex]->Draw(); //A=>without axis. Not sure what Z is, P not applicable
	confCanvas->Update();
      }
    }
    
    else{
      electron->widthFit[yIndex] = NULL;
    }
    
    cout <<"Finished Parameterizing the Electron mean and width!" <<endl;
    
    //----- ROUND 2: PARAMETERIZE THE ELECTRON AMPLITUDE
    for (Int_t mTm0Index=2; mTm0Index<nmTm0Bins; mTm0Index++){
      
      //Make sure that the electron hasn't already been killed
      if (killElectron[yIndex])
	continue;
      
      //Make sure the histogram exists
      if (!tpcProtonMinusHisto.at(yIndex).at(mTm0Index))
	continue;
      
      //Pass the histogram to a local pointer
      yieldHisto = tpcProtonMinusHisto.at(yIndex).at(mTm0Index);
      yieldHisto->SetAxisRange(-2.,2.5);
      
      //Make sure it has enough entries
      if (yieldHisto->GetEntries() < 100)
	continue;
      
      //Compute the kinematics for this bin
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      
      //Make the Predictions for the particle species
      pionPredict = particleInfo->PredictZTPC(pTotal,PROTON,PION);
      kaonPredict = particleInfo->PredictZTPC(pTotal,PROTON,KAON);
      protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);
      electronPredict = particleInfo->PredictZTPC(pTotal,PROTON,ELECTRON);
      
      //Stop at this mTm0
      if (mTm0 > mTm0Cutoff)
	break;
      
      //Get the Previous Pion Mean
      Double_t pionMeanPrev = pion->mean[yIndex]->Eval(mTm0);
      
      //Create the Fit and pass it to a local pointer
      if (tpcProtonMinusFit.at(yIndex).at(mTm0Index))
	tpcProtonMinusFit.at(yIndex).at(mTm0Index)->Delete();
      
      tpcProtonMinusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)+gaus(3)",pionPredict-0.3,electronPredict+0.3);
      yieldFit = tpcProtonMinusFit.at(yIndex).at(mTm0Index);
      
      //Pion Parameters
      yieldFit->SetParName(0,"#pi_{Amp}");
      yieldFit->SetParName(1,"#mu_{#pi}");
      yieldFit->SetParName(2,"{#sigma}_{#pi}");
      yieldFit->SetParName(3,"e_{Amp}");
      yieldFit->SetParName(4,"#mu_{e}");
      yieldFit->SetParName(5,"{#sigma}_{e}");
      
      yieldFit->SetParameter(0,yieldHisto->GetBinContent(yieldHisto->FindBin(pionMeanPrev)));
      yieldFit->FixParameter(1,pionMeanPrev);
      yieldFit->FixParameter(2,pion->widthFit[yIndex]->Eval(mTm0));
      
      //Electron Parameters
      yieldFit->SetParameter(3,yieldHisto->GetBinContent(yieldHisto->FindBin(electronPredict)));
      yieldFit->FixParameter(4,electron->mean[yIndex]->Eval(mTm0));
      yieldFit->FixParameter(5,electron->widthFit[yIndex]->Eval(mTm0));
      
      fittingCanvas->cd();
      yieldHisto->Draw("E");
      //if (DRAWIMG)  yieldHisto->Draw("E");
      yieldHisto->Fit(yieldFit,"RQ");
      lineHeight = yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict));
      if (DRAWIMG) DrawPredictionLines(fittingCanvas, chargePID);
      fittingCanvas->Update();
      
      if (SAVEIMG){
	
	gPad->Update();

	TString fullImageDir = imageDirect+"/ProtonMinus/QA/";
	TString imageName = TString::Format("FinalElectronPionMinus_%s_%02d_%02d_%02d",
					    particleInfo->GetParticleName(2,chargePID).Data(),
					    centralityIndex,yIndex,mTm0Index);
	
	fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	
      }
      
      //Add Point to the Electron Parameter Graph
      AddPointToParGraph(electron->amp[yIndex],yieldFit,3,mTm0Index);
      
      if (DRAWIMG){
	fitParCanvas->cd(4);
	electron->amp[yIndex]->Draw("APZ");
	fitParCanvas->Update();
      }
      
    }//End Loop Over mTm0Index
    
    //----- ROUND 3: Final Fitting Loop
    for (Int_t mTm0Index=2; mTm0Index<nmTm0Bins; mTm0Index++){
      
      //Compute the kinematics for this bin
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      
      //Make sure that the electron hasn't already been killed
      if (killElectron[yIndex] || mTm0 > mTm0Cutoff)
	continue;
      
      //Make sure the histogram exists
      if (!tpcProtonMinusHisto.at(yIndex).at(mTm0Index))
	continue;
      
      //Pass the histogram to a local pointer
      yieldHisto = tpcProtonMinusHisto.at(yIndex).at(mTm0Index);
      
      //Make sure it has enough entires
      if (yieldHisto->GetEntries() < 100)
	continue;
      
      //Make the Predictions for the particle species
      pionPredict = particleInfo->PredictZTPC(pTotal,PROTON,PION);
      kaonPredict = particleInfo->PredictZTPC(pTotal,PROTON,KAON);
      protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);
      electronPredict = particleInfo->PredictZTPC(pTotal,PROTON,ELECTRON);
      
      //Get the Previous Pion Mean
      Double_t pionMeanPrev = pion->mean[yIndex]->Eval(mTm0);
      
      //Draw the Yield Histogram
      yieldHisto->SetAxisRange(-2.,2.5);
      yieldHisto->SetMarkerStyle(kOpenCircle);
      fittingCanvas->cd();
      yieldHisto->Draw("E");
      lineHeight = yieldHisto->GetBinContent(yieldHisto->FindBin(pionPredict));
      fittingCanvas->Update();
      
      //Create the Fit and pass it to a local pointer
      if (tpcProtonMinusFit.at(yIndex).at(mTm0Index))
	tpcProtonMinusFit.at(yIndex).at(mTm0Index)->Delete();
      
      tpcProtonMinusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)+gaus(3)",pionPredict-0.3,electronPredict+0.3);
      yieldFit = tpcProtonMinusFit.at(yIndex).at(mTm0Index);
      
      //Set Parameter Names
      yieldFit->SetParNames(Form("%s_{Amp.}",particleInfo->GetParticleSymbol(PION).Data()),
			    "#mu_{#pi}","#sigma_{#pi}","e_{Amp.}",
			    "#mu_{e}","#sigma_{e}");
      
      //Pion Parameters
      yieldFit->SetParameter(0,yieldHisto->GetBinContent(yieldHisto->FindBin(pionMeanPrev)));
      yieldFit->FixParameter(1,pionMeanPrev);
      yieldFit->FixParameter(2,pion->widthFit[yIndex]->Eval(mTm0));
      
      //Electron Parameters
      Double_t eAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(electron->mean[yIndex]->Eval(mTm0)));
      if (mTm0 >= electron->amp[yIndex]->GetX()[0] &&
	  mTm0 <= electron->amp[yIndex]->GetX()[electron->amp[yIndex]->GetN()-1])
	eAmpSeed = electron->amp[yIndex]->Eval(mTm0);
      yieldFit->FixParameter(3,eAmpSeed);
      yieldFit->FixParameter(4,electron->mean[yIndex]->Eval(mTm0));
      yieldFit->FixParameter(5,electron->widthFit[yIndex]->Eval(mTm0));
      
      //Fit the Yield Histogram
      yieldHisto->Fit(yieldFit,"RQ");
      fittingCanvas->Update();
      
      
      if (SAVEIMG){
        TString fullImageDir = imageDirect+"/ProtonMinus/QA/";
        TString imageName = TString::Format("finalfitMinus_%s_%02d_%02d_%02d",
                                            particleInfo->GetParticleName(2,chargePID).Data(),
                                            centralityIndex,yIndex,mTm0Index);
	
        fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	
      }
    }       //End Loop Over mTm0Index
    
    cout <<"Finished Final Loop of Pion Minus Fitting Procedure" <<endl;
    
  }      //End Loop Over yIndex
} // End Fxn

void FitProtonPlus(TString imageDirect){

  //Make QA directory in /ProtonPlus/
  if (!gSystem->OpenDirectory(imageDirect+"/ProtonPlus/QA/")){
    gSystem->mkdir(imageDirect+"/ProtonPlus/QA/");
  }
   
  //Local Variables
  Int_t chargePID = 1;
  TH1D *yieldHisto = NULL;
  TF1 *yieldFit = NULL;
  
  for (Int_t yIndex=firstRapidityIndex; yIndex<lastRapidityIndex; yIndex++){
    
    if (!protonWidthFit.at(yIndex)) continue;
    
    //If the user has specified a rapidityIndex then only do that index
    if (rapidityIndex >= 0 && yIndex != rapidityIndex)
      continue;
    
    //Compute the Rapidity
    rapidity = GetRapidityRangeCenter(yIndex);
    
    cout <<"Working on Rapidity Index: " <<yIndex <<" Rapidity: " <<rapidity <<endl;
    
    //-----ROUND 1 - Parameterize the Mean of the Proton
    for (Int_t mTm0Index=2; mTm0Index<nmTm0Bins; mTm0Index++){
      
      cout<<"mTm0Index: "<<mTm0Index<<endl;
      
      //Make sure that proton hasn't already been killed
      if (killFirstProtonAttempt[yIndex]) continue;
      
      //Make sure the histogram exists
      if (!tpcProtonPlusHisto.at(yIndex).at(mTm0Index))
	continue;
      
      //Pass the histogram to a local pointer
      yieldHisto = tpcProtonPlusHisto.at(yIndex).at(mTm0Index);
      yieldHisto->SetAxisRange(-2.,2.5);
      
      
      //Make sure it has enough entries
      if (yieldHisto->GetEntries() < 100)
	continue;
      
      //Compute the kinematics for this bin
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      
      //Make the Predictions for the particle species
      pionPredict = particleInfo->PredictZTPC(pTotal,PROTON,PION);
      kaonPredict = particleInfo->PredictZTPC(pTotal,PROTON,KAON);
      protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);
      electronPredict = particleInfo->PredictZTPC(pTotal,PROTON,ELECTRON);
      deuteronPredict = particleInfo->PredictZTPC(pTotal,PROTON,DEUTERON);
      
      if (mTm0 > 1.0)
	break;
      
      
      //Draw the Histogram
      fittingCanvas->cd();
      
      yieldHisto->Draw("E");
      
      //Delete the Previous fit
      if (tpcProtonPlusFit.at(yIndex).at(mTm0Index))
	tpcProtonPlusFit.at(yIndex).at(mTm0Index)->Delete();
      
      //Fit range
      Double_t min      = TMath::Min(pionPredict, protonPredict);
      Double_t max      = TMath::Max(pionPredict, protonPredict);
      Double_t minBelow = FindMinBinCenter(yieldHisto,min-.3,min);
      Double_t minAbove = FindMinBinCenter(yieldHisto,1,1.3);
      
      //Fit all four particles.......
      if (!killElectron[yIndex] && !skipPiMinus[yIndex] && mTm0 < mTm0Cutoff){
	
	cout<<"Electron not killed"<<endl;

	//Pull width and mean from zTPCProtonTofCut fit 
	Double_t protonWidthSeed = protonWidthFit.at(yIndex)->Eval(mTm0);
	Double_t protonMeanSeed = GetYFromMeanTGraph(protonMeanGraph.at(yIndex),mTm0Index);

	//if the mean seed does not exist, use the bischel prediction
	if (protonMeanSeed == -999) {
	  cout << "INFO:  Proton Mean Seed set to ProtonPredict value" << endl;
	  protonMeanSeed = protonPredict;
	}
	
	//Amplitude set to histo value in bin of protonMeanSeed
	Double_t protonAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(protonMeanSeed));
	
	//If the statistics are low, it is better to find the highest value in the histogram and set proton amp to the value, the fit tends to underpredict amplitude
	if ( LowStatAmpSeed(yieldHisto) < 40 ){
	  protonAmpSeed = LowStatAmpSeed(yieldHisto);
	  cout << "Low Stats... Setting ProtonAmp to highest bin" << endl;
	  cout << "AMP SEED = " << protonAmpSeed << endl;
	}
	
	//Pull width and mean from zTPCDeuteronTofCut fit 
	Double_t deuteronWidthSeed = deuteronWidthFit.at(yIndex)->Eval(mTm0);
	Double_t deuteronMeanSeed = GetYFromMeanTGraph(deuteronMeanGraph.at(yIndex),mTm0Index);

	if (deuteronMeanSeed == -999){
	  cout << "INFO:  Deuteron Mean Seed set to DeuteronPredict value" << endl;
	  deuteronMeanSeed = deuteronPredict;
        }

	//Amplitude set to histo value in bin of deuteronMeanSeed
	Double_t deuteronAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(deuteronMeanSeed));

	tpcProtonPlusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)+gaus(6) + gaus(3) + gaus(9)",minBelow,minAbove);
	yieldFit = tpcProtonPlusFit.at(yIndex).at(mTm0Index);
	
	yieldFit->SetParName(0,"#pi_{Amp}");
        yieldFit->SetParName(1,"#mu_{#pi}");
        yieldFit->SetParName(2,"{#sigma}_{#pi}");
        yieldFit->SetParName(6,"e_{Amp}");
        yieldFit->SetParName(7,"#mu_{e}");
        yieldFit->SetParName(8,"{#sigma}_{e}");
        yieldFit->SetParName(3,"p_{Amp}");
        yieldFit->SetParName(4,"#mu_{p}");
        yieldFit->SetParName(5,"{#sigma}_{p}");
        yieldFit->SetParName(9,"d_{Amp} ratio");
        yieldFit->SetParName(10,"#mu_{d}");
        yieldFit->SetParName(11,"{#sigma}_{d}");

	//Get the previous pion parameters
	Double_t piPlusAmp     = GetYFromTGraph(pion->amp[yIndex], mTm0Index);
	Double_t pionMeanPrev  = pion->mean[yIndex]->Eval(mTm0);
	Double_t pionWidthPrev = pion->widthFit[yIndex]->Eval(mTm0);
	Bool_t PiOverlap = (abs(pionMeanPrev - protonPredict) < (2*protonWidthSeed) )  ||  (abs(pionMeanPrev - deuteronPredict) < (2*deuteronWidthSeed) );  
	
	if (PiOverlap) yieldFit->FixParameter(0, piPlusAmp);
	else yieldFit->SetParameter(0, piPlusAmp);

	yieldFit->FixParameter(1, pionMeanPrev);
	yieldFit->FixParameter(2, pionWidthPrev);
	
	//Electron Parameters
	if (fabs(electron->mean[yIndex]->Eval(mTm0)) <
	    (electron->widthFit[yIndex]->Eval(mTm0)+protonWidthSeed) ||
	    fabs(electron->mean[yIndex]->Eval(mTm0) - pion->mean[yIndex]->Eval(mTm0)) <
	    (electron->widthFit[yIndex]->Eval(mTm0)+pion->widthFit[yIndex]->Eval(mTm0)))
	  
	  yieldFit->FixParameter(6,electron->amp[yIndex]->Eval(mTm0));
	
	else 
	  yieldFit->SetParameter(6,
				 yieldHisto->GetBinContent(yieldHisto->
							   FindBin(electron->
								   mean[yIndex]->Eval(mTm0))));
	
	yieldFit->FixParameter(7,electron->mean[yIndex]->Eval(mTm0));
	yieldFit->FixParameter(8,electron->widthFit[yIndex]->Eval(mTm0));
	
	yieldFit->SetParameter(9,deuteronAmpSeed);
        yieldFit->SetParameter(10, deuteronMeanSeed);
	yieldFit->SetParLimits(10,(deuteronMeanSeed-deuteronWidthSeed),(deuteronMeanSeed+deuteronWidthSeed));
        yieldFit->FixParameter(11,deuteronWidthSeed);
	
	//Proton Parameters
	yieldFit->SetParameter(3,protonAmpSeed);
	yieldFit->SetParameter(4,protonMeanSeed);
	yieldFit->SetParLimits(4,(protonMeanSeed-protonWidthSeed),(protonMeanSeed+protonWidthSeed));
        yieldFit->FixParameter(5,protonWidthSeed);
		
        yieldFit->SetName(Form("TPCPionPlusFit_%02d_%02d",yIndex,mTm0Index));

          //Check to see if proton peak was correctly fit
          bool protonGoodFit =false;

          if ( abs( abs(protonMeanSeed) - abs(deuteronMeanSeed) ) < protonMeanSeed ){

	    deuteronAmpSeed = protonAmpSeed*ratioFit.at(yIndex)->Eval(mTm0);

            //Deuteron Parameters                                                                   
            yieldFit->SetParameter(9,deuteronAmpSeed);
            yieldFit->FixParameter(10, deuteronMeanSeed);
	    yieldFit->FixParameter(11,deuteronWidthSeed);

            //Proton Parameters 
            yieldFit->SetParameter(3,protonAmpSeed);
            yieldFit->SetParameter(4,protonMeanSeed);
            yieldFit->SetParLimits(4,(protonMeanSeed-protonWidthSeed),(protonMeanSeed+protonWidthSeed));
            yieldFit->FixParameter(5,protonWidthSeed);
          }

          yieldHisto->Fit(yieldFit,"RQ");
	AddPointToParGraph(proton->mean[yIndex],yieldFit,4,mTm0Index);	

	//refit if bad fit
          if (yieldFit->GetParameter(3) > yieldFit->GetParameter(9) ) protonGoodFit = true;
          if(protonGoodFit == false){
	    cout << "WARNING: BAD PROTON FIT, HIGH DEUT" << endl;
	    cout << "Attempting better fit!" << endl;
	    float nudge = .02;
	    protonMeanSeed = protonMeanSeed += nudge;

	    yieldFit->SetParameter(9,deuteronAmpSeed);
	    yieldFit->SetParLimits(9,deuteronAmpSeed*.6,deuteronAmpSeed*.12);
            yieldFit->FixParameter(10, deuteronMeanSeed);
            yieldFit->FixParameter(11,deuteronWidthSeed);

	    yieldFit->SetParameter(3,protonAmpSeed);
            yieldFit->SetParameter(4,protonMeanSeed);
            yieldFit->SetParLimits(4,2*(protonMeanSeed-protonWidthSeed),2*(protonMeanSeed+protonWidthSeed));
            yieldFit->FixParameter(5,protonWidthSeed);

	    yieldHisto->Fit(yieldFit,"RQ");

	  }

	DrawPredictionLines(fittingCanvas, chargePID);
        fittingCanvas->Update();
	
	if (SAVEIMG){

	  gPad->Update();
	  TString fullImageDir = imageDirect+"/ProtonPlus/QA/";
	  TString imageName = TString::Format("Fit_9_%s_%02d_%02d_%02d",
					      particleInfo->GetParticleName(2,chargePID).Data(),
					      centralityIndex,yIndex,mTm0Index);
	  
	  fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	  
	}
	
	
	
      }//End if electron is not killed
      
     
//Second case.... fit 3 particles but not the electron

      else if (!skipPiMinus[yIndex] && ( killElectron[yIndex] || (!killElectron[yIndex] && mTm0 > mTm0Cutoff) )){
        cout<<"Electron killed but not proton or pion"<<endl;

	//Pull width and mean from zTPCProtonTofCut fit 	
	Double_t protonWidthSeed = protonWidthFit.at(yIndex)->Eval(mTm0);
	Double_t protonMeanSeed = GetYFromMeanTGraph(protonMeanGraph.at(yIndex),mTm0Index);

	if (protonMeanSeed == -999){
	  cout << "INFO:  Proton Mean Seed set to ProtonPredict value" << endl;  
	  protonMeanSeed = protonPredict;
	}

	//Amplitude set to histo value in bin of protonMeanSeed
	Double_t protonAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(protonMeanSeed));

	//check if low statistics 
	if ( LowStatAmpSeed(yieldHisto) < 40 ){
	  protonAmpSeed = LowStatAmpSeed(yieldHisto);
	  cout << "Low Stats... Setting ProtonAmp to highest bin" << endl;
	  cout << "AMP SEED = " << protonAmpSeed << endl;
	}

	//Pull width and mean from zTPCDeuteronTofCut fit 	
        Double_t deuteronWidthSeed = deuteronWidthFit.at(yIndex)->Eval(mTm0);

	Double_t deuteronMeanSeed = GetYFromMeanTGraph(deuteronMeanGraph.at(yIndex),mTm0Index);

	if (deuteronMeanSeed == -999){
	  cout << "INFO:  Deuteron Mean Seed set to DeuteronPredict value" << endl;
	  deuteronMeanSeed = deuteronPredict;
        }

	//Amplitude set to histo value in bin of deuteronMeanSeed
	Double_t deuteronAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(deuteronMeanSeed));

	tpcProtonPlusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0) + gaus(3) + gaus(9)",minBelow,minAbove);
	yieldFit = tpcProtonPlusFit.at(yIndex).at(mTm0Index);
	
	yieldFit->SetParName(0,"#pi_{Amp}");
        yieldFit->SetParName(1,"#mu_{#pi}");
        yieldFit->SetParName(2,"{#sigma}_{#pi}");
	yieldFit->SetParName(3,"p_{Amp}");
        yieldFit->SetParName(4,"#mu_{p}");
        yieldFit->SetParName(5,"{#sigma}_{p}");
        yieldFit->SetParName(9,"d_{Amp} ratio");
        yieldFit->SetParName(10,"#mu_{d}");
        yieldFit->SetParName(11,"{#sigma}_{d}");

	//Pion Parameters
	Double_t piPlusAmp   = GetYFromTGraph(pion->amp[yIndex], mTm0Index);
	Double_t pionMeanPrev  = pion->mean[yIndex]->Eval(mTm0);
	Double_t pionWidthPrev = pion->widthFit[yIndex]->Eval(mTm0);
	Bool_t PiOverlap = (abs(pionMeanPrev - protonPredict) < (2*protonWidthSeed) )  || (abs(pionMeanPrev - deuteronPredict) < (2*deuteronWidthSeed) );  
	
        if (PiOverlap) yieldFit->FixParameter(0, piPlusAmp);
        else yieldFit->SetParameter(0, piPlusAmp);
	
	yieldFit->FixParameter(1, pionMeanPrev);
	yieldFit->FixParameter(2, pionWidthPrev);
	
	//Deuteron Parameters
	yieldFit->SetParameter(9,deuteronAmpSeed);
        yieldFit->SetParameter(10, deuteronMeanSeed);
	yieldFit->SetParLimits(10,(deuteronMeanSeed-deuteronWidthSeed),(deuteronMeanSeed+deuteronWidthSeed));
        yieldFit->FixParameter(11,deuteronWidthSeed);
	
        //Proton Parameters
        yieldFit->SetParameter(3,protonAmpSeed);
	yieldFit->SetParameter(4,protonMeanSeed);
        yieldFit->SetParLimits(4,(protonMeanSeed-protonWidthSeed),(protonMeanSeed+protonWidthSeed));
	yieldFit->FixParameter(5,protonWidthSeed);
	yieldFit->SetName(Form("TPCPionPlusFit_%02d_%02d",yIndex,mTm0Index));
	
	//check if proton and deuteron are overlapping 
	bool protonGoodFit =false;
	
	if ( abs( abs(protonMeanSeed) - abs(deuteronMeanSeed) ) < protonMeanSeed ){
	  
	  deuteronAmpSeed = protonAmpSeed*ratioFit.at(yIndex)->Eval(mTm0);
	  
	  //Deuteron Parameters
	  yieldFit->SetParameter(9,deuteronAmpSeed);
	  yieldFit->FixParameter(10, deuteronMeanSeed);
	  yieldFit->FixParameter(11,deuteronWidthSeed);
	  
	  //Proton Parameters
	  yieldFit->SetParameter(3,protonAmpSeed);
	  yieldFit->SetParameter(4,protonMeanSeed);
	  yieldFit->SetParLimits(4,(protonMeanSeed-protonWidthSeed),(protonMeanSeed+protonWidthSeed));
	  yieldFit->FixParameter(5,protonWidthSeed);
	}
	
	//Fit the Yield Histogram
	yieldHisto->Fit(yieldFit,"RQ");
	AddPointToParGraph(proton->mean[yIndex],yieldFit,4,mTm0Index);	
	
	
	//check quality of fit
	if (yieldFit->GetParameter(3) > yieldFit->GetParameter(9) ) protonGoodFit = true;
	if(protonGoodFit == false){
	  cout << "WARNING: BAD PROTON FIT, HIGH DEUT" << endl;
	  cout << "Attempting better fit!" <<endl;
	  float nudge = .02;
	  protonMeanSeed = protonMeanSeed += nudge;
	  
	  yieldFit->SetParameter(9,deuteronAmpSeed);
	  yieldFit->SetParLimits(9,deuteronAmpSeed*.6,deuteronAmpSeed*.12);
	  yieldFit->FixParameter(10, deuteronMeanSeed);
	  yieldFit->FixParameter(11,deuteronWidthSeed);
	  
	  yieldFit->SetParameter(3,protonAmpSeed);
	  yieldFit->SetParameter(4,protonMeanSeed);
	  yieldFit->SetParLimits(4,2*(protonMeanSeed-protonWidthSeed),2*(protonMeanSeed+protonWidthSeed));
	  yieldFit->FixParameter(5,protonWidthSeed);
	  
	  yieldHisto->Fit(yieldFit,"RQ");
	  
	}
	
	DrawPredictionLines(fittingCanvas, chargePID);
	fittingCanvas->Update();
	
	if (SAVEIMG){
	  
	  gPad->Update();
	  
	  TString fullImageDir = imageDirect+"/ProtonPlus/QA/";
	  TString imageName = TString::Format("Fit_10_%s_%02d_%02d_%02d",
					      particleInfo->GetParticleName(2,chargePID).Data(),
					      centralityIndex,yIndex,mTm0Index);
	  
	  fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	  
	}
      }
      
      //Fit 2 particles... proton and deuteron
      
      if (skipPiMinus[yIndex]){
        cout<<"Pion killed but not proton"<<endl;
	
	//Pull width and mean from zTPCProtonTofCut fit 		
        Double_t protonWidthSeed = protonWidthFit.at(yIndex)->Eval(mTm0);
	Double_t protonMeanSeed = GetYFromMeanTGraph(protonMeanGraph.at(yIndex),mTm0Index);

	//check if fit exists, if not fit with bischel value
	if (protonMeanSeed == -999){
	  cout << "INFO:  Proton Mean Seed set to ProtonPredict value" << endl;
	  protonMeanSeed = protonPredict;
	}

	//Amplitude set to histo value in bin of protonMeanSeed
	Double_t protonAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(protonMeanSeed));

	//Check if low statistics 
	if ( LowStatAmpSeed(yieldHisto) < 40 ){
	  protonAmpSeed = LowStatAmpSeed(yieldHisto);
	  cout << "Low Stats... Setting ProtonAmp to highest bin" << endl;
	  cout << "AMP SEED = " << protonAmpSeed << endl;
	}

	//Pull width and mean from zTPCDeuteronTofCut fit 		
        Double_t deuteronWidthSeed = deuteronWidthFit.at(yIndex)->Eval(mTm0);
	Double_t deuteronMeanSeed = GetYFromMeanTGraph(deuteronMeanGraph.at(yIndex),mTm0Index);

	if (deuteronMeanSeed == -999){
	  cout << "INFO:  Deuteron Mean Seed set to DeuteronPredict value" << endl;
	  deuteronMeanSeed = deuteronPredict;
	}

	//Amplitude set to histo value in bin of deuteronMeanSeed
        Double_t deuteronAmpSeed = yieldHisto->GetBinContent(yieldHisto->FindBin(deuteronMeanSeed));

        tpcProtonPlusFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)+gaus(3)",minBelow,minAbove);
        yieldFit = tpcProtonPlusFit.at(yIndex).at(mTm0Index);

        yieldFit->SetParName(0, "p_{Amp.}");
        yieldFit->SetParName(1, "#mu_{p}");
        yieldFit->SetParName(2, "#sigma_{p}");

	yieldFit->SetParName(3, "d_{Amp.}");
        yieldFit->SetParName(4, "#mu_{deut}");
        yieldFit->SetParName(5, "#sigma_{deut}");

	//Deuteron Parameters
        yieldFit->SetParameter(3,deuteronAmpSeed);
        yieldFit->SetParameter(4, deuteronMeanSeed);
	yieldFit->SetParLimits(4,(deuteronMeanSeed-deuteronWidthSeed),(deuteronMeanSeed+deuteronWidthSeed));
	yieldFit->FixParameter(5,deuteronWidthSeed);

        //Proton Parameters
        yieldFit->SetParameter(0,protonAmpSeed);
	yieldFit->SetParameter(1,protonMeanSeed);
        yieldFit->SetParLimits(1,(protonMeanSeed-protonWidthSeed),(protonMeanSeed+protonWidthSeed));
	yieldFit->FixParameter(2,protonWidthSeed);
        yieldFit->SetName(Form("TPCPionPlusFit_%02d_%02d",yIndex,mTm0Index));


	//check if proton and deuteron seeds overlap
	bool protonGoodFit =false;

	if ( abs( abs(protonMeanSeed) - abs(deuteronMeanSeed) ) < protonMeanSeed ){

	  deuteronAmpSeed = protonAmpSeed*ratioFit.at(yIndex)->Eval(mTm0);

	  //Deuteron Parameters                                                                                                                                                                           

	  yieldFit->SetParameter(3,deuteronAmpSeed);
	  yieldFit->FixParameter(4, deuteronMeanSeed);
	  yieldFit->FixParameter(5,deuteronWidthSeed);

	  //Proton Parameters                                                                                                                                                                             
	  yieldFit->SetParameter(0,protonAmpSeed);
	  yieldFit->SetParameter(1,protonMeanSeed);
	  yieldFit->SetParLimits(1,(protonMeanSeed-protonWidthSeed),(protonMeanSeed+protonWidthSeed));
	  yieldFit->FixParameter(2,protonWidthSeed);
	}

        //Fit the Yield Histogram
	yieldHisto->Fit(yieldFit,"RQ");
	AddPointToParGraph(proton->mean[yIndex],yieldFit,1,mTm0Index);	

	//Refit if needed
	if (yieldFit->GetParameter(3) > yieldFit->GetParameter(9) ) protonGoodFit = true;
	if(protonGoodFit == false){
	  cout << "WARNING: BAD PROTON FIT, HIGH DEUT" << endl;
	  cout << "Attempting better fit!" <<endl;
	  float nudge = .02;
	  protonMeanSeed = protonMeanSeed += nudge;

	  yieldFit->SetParameter(3,deuteronAmpSeed);
	  yieldFit->SetParLimits(3,deuteronAmpSeed*.6,deuteronAmpSeed*.12);
	  yieldFit->FixParameter(4, deuteronMeanSeed);
	  yieldFit->FixParameter(5,deuteronWidthSeed);

	  yieldFit->SetParameter(0,protonAmpSeed);
	  yieldFit->SetParameter(1,protonMeanSeed);
	  yieldFit->SetParLimits(1,2*(protonMeanSeed-protonWidthSeed),2*(protonMeanSeed+protonWidthSeed));
	  yieldFit->FixParameter(2,protonWidthSeed);

	  yieldHisto->Fit(yieldFit,"RQ");

	}
      }

      DrawPredictionLines(fittingCanvas, chargePID);
        fittingCanvas->Update();
	
	if (SAVEIMG){

	  //	  yieldHisto->SetStats(0);
	  gPad->Update();

          TString fullImageDir = imageDirect+"/ProtonPlus/QA/";
          TString imageName = TString::Format("Fit_11_%s_%02d_%02d_%02d",
                                              particleInfo->GetParticleName(2,chargePID).Data(),
                                              centralityIndex,yIndex,mTm0Index);
	  
          fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	  
        }
    }//End Loop Over mTm0Index
  }//End Loop Over yIndex
}

//--------------------------------------------------------------------
void DrawPredictionLines(TCanvas *canvas, Int_t charge){

  canvas->cd();
  
  PredictionLine(pionLine,pionPredict,lineHeight)->Draw("SAME");
  PredictionLine(electronLine,electronPredict,lineHeight)->Draw("SAME");
  PredictionLine(kaonLine,kaonPredict,lineHeight)->Draw("SAME");
  PredictionLine(electronLine,electronPredict,lineHeight)->Draw("SAME");
	PredictionLine(protonLine,protonPredict,lineHeight)->Draw("SAME");
	if (charge > 0){
		PredictionLine(deuteronLine,deuteronPredict,lineHeight)->Draw("SAME");
  }
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
Bool_t AddPointToSpectra(TH1D *yieldHisto, TF1 *yieldFit, TGraphErrors *rawSpectra,TGraphErrors *ratioTPCdeut, TGraphErrors * ratioTPCpion, 
			 Int_t yIndex, Int_t mTm0Index, Int_t PID, Int_t CHARGE, TString IMAGEDIR){
  
  //Adds a point to the spectra
  //returns false on fail, true on success
  
  //Make sure both the histogram and fit exist
  if (!yieldHisto || !yieldFit) return false;
  
  
  rapidity = GetRapidityRangeCenter(yIndex);
  
  ComputeBinKinematics(mTm0Index,rapidity,PROTON);
  
  //Make the Predictions for the particle species
  pionPredict = particleInfo->PredictZTPC(pTotal,PROTON,PION);
  kaonPredict = particleInfo->PredictZTPC(pTotal,PROTON,KAON);
  protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);
  electronPredict = particleInfo->PredictZTPC(pTotal,PROTON,ELECTRON);
  deuteronPredict = particleInfo->PredictZTPC(pTotal,PROTON,DEUTERON);
  
  
  //Draw the Final Version of the Fit
  fittingCanvas->Clear();
  fittingCanvas->cd();
  fittingCanvas->SetTopMargin(.05);
  fittingCanvas->SetRightMargin(.05);
  yieldHisto->SetMarkerStyle(kOpenCircle);
  yieldHisto->SetTitle("");
  yieldHisto->SetAxisRange(-2.0,2.0);
  yieldHisto->Draw("E");
  
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
  TF1 protonGaus("proton","gaus(0)",-20,20);
  TF1 deuteronGaus("deuteron","gaus(0)",-20,20);
  
  if (!skipPiMinus[yIndex]){
        
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
    
    Double_t deuteronAmp =  (yieldFit->GetParameter(9) );
    deuteronGaus.FixParameter(0,deuteronAmp);
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
    
  }
  
  if (skipPiMinus[yIndex]){
    
    protonGaus.FixParameter(0,yieldFit->GetParameter(0));
    protonGaus.FixParameter(1,yieldFit->GetParameter(1));
    protonGaus.FixParameter(2,yieldFit->GetParameter(2));
    protonGaus.SetParError(0,yieldFit->GetParError(0));
    protonGaus.SetParError(1,yieldFit->GetParError(1));
    protonGaus.SetParError(2,yieldFit->GetParError(2));
    protonGaus.SetLineColor(particleInfo->GetParticleColor(PROTON));
    protonGaus.SetLineWidth(3);
    protonGaus.SetLineStyle(7);
    protonGaus.SetNpx(10000);
    protonGaus.Draw("SAME");
    
    Double_t deuteronAmp =  (yieldFit->GetParameter(3) );
    deuteronGaus.FixParameter(0,deuteronAmp);
    deuteronGaus.FixParameter(1,yieldFit->GetParameter(4));
    deuteronGaus.FixParameter(2,yieldFit->GetParameter(5));
    deuteronGaus.SetParError(0,yieldFit->GetParError(3));
    deuteronGaus.SetParError(1,yieldFit->GetParError(4));
    deuteronGaus.SetParError(2,yieldFit->GetParError(5));
    
    deuteronGaus.SetLineColor(particleInfo->GetParticleColor(DEUTERON));
    deuteronGaus.SetLineWidth(3);
    deuteronGaus.SetLineStyle(7);
    deuteronGaus.SetNpx(10000);
    deuteronGaus.Draw("SAME");
      
  }

  TF1 electronGaus("electron","gaus(0)",-20,20);
  if (!killElectron[yIndex]){
    
    electronGaus.FixParameter(0,yieldFit->GetParameter(6));
    electronGaus.FixParameter(1,yieldFit->GetParameter(7));
    electronGaus.FixParameter(2,yieldFit->GetParameter(8));
    electronGaus.SetParError(0,yieldFit->GetParError(6));
    electronGaus.SetParError(1,yieldFit->GetParError(7));
    electronGaus.SetParError(2,yieldFit->GetParError(8));
    electronGaus.SetLineColor(particleInfo->GetParticleColor(ELECTRON));
    electronGaus.SetLineWidth(3);
    electronGaus.SetLineStyle(7);
    electronGaus.SetNpx(10000);
    electronGaus.Draw("SAME");
    
  }

  //Get the Stats Box and Set its attributes
  
  TPaveStats *statsBox = (TPaveStats*)yieldHisto->FindObject("stats");
  statsBox->SetX1NDC(.14);
  statsBox->SetY1NDC(.12);
  statsBox->SetX2NDC(.37);
  statsBox->SetY2NDC(.78);
  
  //yieldHisto->SetStats(0);
  
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
  if (!killElectron[yIndex]) leg.AddEntry("electron","e","L");
  leg.AddEntry("proton","p","L");
  leg.SetFillColor(kWhite);
  leg.SetBorderSize(0);
  leg.Draw("SAME");
  
  DrawPredictionLines(fittingCanvas,1);
  
  fittingCanvas->Update();
  
  //Extract the Yield
  Double_t amp = protonGaus.GetParameter(0);
  Double_t ampError = protonGaus.GetParError(0);
  Double_t width = protonGaus.GetParameter(2);
  Double_t widthError = protonGaus.GetParError(2);
  Double_t pi = TMath::Pi();

  Double_t zBinWidth = (Double_t)yieldHisto->GetBinWidth(2);
  Double_t mTm0BinCenter = GetmTm0RangeCenter(mTm0Index);
  Double_t nEvents = (Double_t)nEventsHisto->GetBinContent(nEventsHisto->FindBin(centralityIndex));

  //Set the Normalization Factor
  Double_t normFactor = (1.0/(mTm0BinCenter+particleInfo->GetParticleMass(PID))) * (1.0/mTm0BinWidth) * (1.0/rapidityBinWidth) * (1.0/nEvents) * (1.0/(2.0*pi));
 
  //Compute the Raw Yield
  Double_t rawYield, rawYieldError;
  rawYield = ( (amp * width * sqrt(2.0*pi)) / zBinWidth ) * normFactor;
  rawYieldError = ( (amp * width * sqrt(pow(ampError/amp,2)+pow(widthError/width,2)) *
         sqrt(2.0*pi)) / zBinWidth) * normFactor;

  Double_t deut_amp = deuteronGaus.GetParameter(0);
  Double_t deut_ampError = deuteronGaus.GetParError(0);
  Double_t deut_width = deuteronGaus.GetParameter(2);
  Double_t deut_widthError = deuteronGaus.GetParError(2);

  Double_t pion_amp = pionGaus.GetParameter(0);
  Double_t pion_ampError = pionGaus.GetParError(0);
  Double_t pion_width = pionGaus.GetParameter(2);
  Double_t pion_widthError = pionGaus.GetParError(2);
  
  //Set the Normalization Factor
  Double_t deut_normFactor = (1.0/(mTm0BinCenter+particleInfo->GetParticleMass(PID))) * (1.0/mTm0BinWidth) * (1.0/rapidityBinWidth) * (1.0/nEvents) * (1.0/(2.0*pi));
  Double_t pion_normFactor = (1.0/(mTm0BinCenter+particleInfo->GetParticleMass(PID))) * (1.0/mTm0BinWidth) * (1.0/rapidityBinWidth) * (1.0/nEvents) * (1.0/(2.0*pi));
  
  //Compute Deuteron Raw Yield
  Double_t deut_rawYield, deut_rawYieldError;
  deut_rawYield = ( (deut_amp * deut_width * sqrt(2.0*pi)) / zBinWidth ) * deut_normFactor;
  deut_rawYieldError = ( (deut_amp * deut_width * sqrt(pow(deut_ampError/deut_amp,2)+pow(deut_widthError/deut_width,2)) *
		     sqrt(2.0*pi)) / zBinWidth) * deut_normFactor;

  Double_t pion_rawYield, pion_rawYieldError;
  pion_rawYield = ( (pion_amp * pion_width * sqrt(2.0*pi)) / zBinWidth ) * pion_normFactor;
  pion_rawYieldError = ( (pion_amp * pion_width * sqrt(pow(pion_ampError/pion_amp,2)+pow(pion_widthError/pion_width,2)) *
		     sqrt(2.0*pi)) / zBinWidth) * pion_normFactor;
  
  if (rawYield > 0)
    {
      Double_t ratioDeut = deut_rawYield/rawYield;
      Double_t ratioErrorDeut = sqrt(deut_rawYieldError*deut_rawYieldError + rawYieldError*rawYieldError);

      Double_t ratioPion = pion_rawYield/rawYield;
      Double_t ratioErrorPion = sqrt(pion_rawYieldError*pion_rawYieldError + rawYieldError*rawYieldError);
      
      ratioTPCdeut->SetPoint(ratioTPCdeut->GetN(),GetmTm0RangeCenter(mTm0Index),deut_rawYield/rawYield);
      ratioTPCdeut->SetPointError(ratioTPCdeut->GetN()-1,0,ratioErrorDeut);

      ratioTPCpion->SetPoint(ratioTPCpion->GetN(),GetmTm0RangeCenter(mTm0Index),pion_rawYield/rawYield);
      ratioTPCpion->SetPointError(ratioTPCpion->GetN()-1,0,ratioErrorPion);

    }
  
  //Add Point to Spectrum 
  rawSpectra->SetPoint(rawSpectra->GetN(),GetmTm0RangeCenter(mTm0Index),rawYield);
  rawSpectra->SetPointError(rawSpectra->GetN()-1,0,rawYieldError);
  
  //Draw the Spectra
  spectraCanvas->cd();
  rawSpectra->Draw("APZ");
  spectraCanvas->Update();
  ratioCanvas->cd(1);
  ratioTPCdeut->Draw("APZ");
  ratioCanvas->cd(2);
  ratioTPCpion->Draw("APZ");
  ratioCanvas->Update();
  
  //Save the Fit Image
  TString fullImageDir = IMAGEDIR+"/"+particleInfo->GetParticleName(PID,CHARGE);
  TString imageName = TString::Format("zTPCFit_%s_%02d_%02d_%02d",
				      particleInfo->GetParticleName(PID,CHARGE).Data(),
				      centralityIndex,yIndex,mTm0Index);
  
  fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
  
  TString imageRatioName = TString::Format("zTPC_Ratio_%s_%02d_%02d_%02d",
				      particleInfo->GetParticleName(PID,CHARGE).Data(),
				      centralityIndex,yIndex,mTm0Index);

  ratioCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageRatioName.Data()));

  return true;
}


//USE TO DRAW QA PLOTS
Bool_t DrawFits(TH1D *yieldHisto, TF1 *yieldFit, Int_t yIndex, 
	      Int_t mTm0Index, Int_t PID, Int_t CHARGE, TString IMAGEDIR){

  if (!yieldHisto || !yieldFit) return false;
  //Draw the Final Version of the Fit

  rapidity = GetRapidityRangeCenter(yIndex);
  
  ComputeBinKinematics(mTm0Index,rapidity,PROTON);
  
  //Make the Predictions for the particle species
  pionPredict = particleInfo->PredictZTPC(pTotal,PROTON,PION);
  kaonPredict = particleInfo->PredictZTPC(pTotal,PROTON,KAON);
  protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);
  electronPredict = particleInfo->PredictZTPC(pTotal,PROTON,ELECTRON);
  deuteronPredict = particleInfo->PredictZTPC(pTotal,PROTON,DEUTERON);
  
  fittingCanvas->Clear();
  fittingCanvas->cd();
  fittingCanvas->SetTopMargin(.05);
  fittingCanvas->SetRightMargin(.05);
  yieldHisto->SetMarkerStyle(kOpenCircle);
  yieldHisto->SetTitle("");
  yieldHisto->SetAxisRange(-2.0,2.0);
  yieldHisto->Draw("E");
  
  //Draw the Fit
  yieldFit->SetLineWidth(3);
  yieldFit->SetLineColor(kBlack);
  yieldFit->SetNpx(10000);
  yieldFit->Draw("SAME");
  
  //Create a Label
  TPaveText label(.12,.8,.35,.92,"BL NDC");
  label.SetTextAlign(11);
  label.AddText(Form("Z(%s) Yield Extraction",particleInfo->GetParticleSymbol(PID,CHARGE).Data()));
  label.AddText(Form("  %.2f #leq y_{%s} < %.2f | %.2f #leq m_{T}-m_{%s} < %.2f",
                     GetRapidityRangeLow(yIndex),particleInfo->GetParticleSymbol(PID,CHARGE).Data(),
                     GetRapidityRangeHigh(yIndex),GetmTm0RangeLow(mTm0Index),
                     particleInfo->GetParticleSymbol(PID,CHARGE).Data(),GetmTm0RangeHigh(mTm0Index)));
  label.GetLine(0)->SetTextSize(.05);
  label.GetLine(1)->SetTextSize(.02);
  label.SetFillColor(kWhite);
  label.SetBorderSize(0);
  label.Draw("SAME");

  //Make A Function for Each of the Particles
  TF1 pionGaus("pion","gaus(0)",-20,20);
  TF1 electronGaus("electron","gaus(0)",-20,20);
  TF1 protonGaus("proton","gaus(0)",-20,20);
  TF1 deuteronGaus("deuteron","gaus(0)",-20,20);

  if ( CHARGE ==-1 )
{
  if (!skipPiMinus[yIndex]){
    
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
  }
  
  if (!killElectron[yIndex]){
    
    electronGaus.FixParameter(0,yieldFit->GetParameter(3));
    electronGaus.FixParameter(1,yieldFit->GetParameter(4));
    electronGaus.FixParameter(2,yieldFit->GetParameter(5));
    electronGaus.SetParError(0,yieldFit->GetParError(3));
    electronGaus.SetParError(1,yieldFit->GetParError(4));
    electronGaus.SetParError(2,yieldFit->GetParError(5));
    electronGaus.SetLineColor(particleInfo->GetParticleColor(ELECTRON));
    electronGaus.SetLineWidth(3);
    electronGaus.SetLineStyle(7);
    electronGaus.SetNpx(10000);
    electronGaus.Draw("SAME");
    
  }
 }
  
  if ( CHARGE == 1 )
    {
      
      if (!skipPiMinus[yIndex]){
	
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

	Double_t deuteronAmp =  ( (yieldFit->GetParameter(9)) );
	deuteronGaus.FixParameter(0,deuteronAmp);
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
	
      }
      
      if (skipPiMinus[yIndex]){
	
	protonGaus.FixParameter(0,yieldFit->GetParameter(0));
	protonGaus.FixParameter(1,yieldFit->GetParameter(1));
	protonGaus.FixParameter(2,yieldFit->GetParameter(2));
	protonGaus.SetParError(0,yieldFit->GetParError(0));
	protonGaus.SetParError(1,yieldFit->GetParError(1));
	protonGaus.SetParError(2,yieldFit->GetParError(2));
	protonGaus.SetLineColor(particleInfo->GetParticleColor(PROTON));
	protonGaus.SetLineWidth(3);
	protonGaus.SetLineStyle(7);
	protonGaus.SetNpx(10000);
	protonGaus.Draw("SAME");

	Double_t deuteronAmp =  ( (yieldFit->GetParameter(0))*(yieldFit->GetParameter(3)) );
	deuteronGaus.FixParameter(0,deuteronAmp);
	deuteronGaus.FixParameter(1,yieldFit->GetParameter(4));
	deuteronGaus.FixParameter(2,yieldFit->GetParameter(5));
	deuteronGaus.SetParError(0,yieldFit->GetParError(3));
	deuteronGaus.SetParError(1,yieldFit->GetParError(4));
	deuteronGaus.SetParError(2,yieldFit->GetParError(5));

	deuteronGaus.SetLineColor(particleInfo->GetParticleColor(DEUTERON));
	deuteronGaus.SetLineWidth(3);
	deuteronGaus.SetLineStyle(7);
	deuteronGaus.SetNpx(10000);
	deuteronGaus.Draw("SAME");

      }

      if (!killElectron[yIndex]){

	electronGaus.FixParameter(0,yieldFit->GetParameter(6));
	electronGaus.FixParameter(1,yieldFit->GetParameter(7));
	electronGaus.FixParameter(2,yieldFit->GetParameter(8));
	electronGaus.SetParError(0,yieldFit->GetParError(6));
	electronGaus.SetParError(1,yieldFit->GetParError(7));
	electronGaus.SetParError(2,yieldFit->GetParError(8));
	electronGaus.SetLineColor(particleInfo->GetParticleColor(ELECTRON));
	electronGaus.SetLineWidth(3);
	electronGaus.SetLineStyle(7);
	electronGaus.SetNpx(10000);
	electronGaus.Draw("SAME");

      }
    }

  //Get the Stats Box and Set its attributes
  TPaveStats *statsBox = (TPaveStats*)yieldHisto->FindObject("stats");
  statsBox->SetX1NDC(.84);
  statsBox->SetY1NDC(.18);
  statsBox->SetX2NDC(.97);
  statsBox->SetY2NDC(.78);
  
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
  if (CHARGE < 0 ) leg.AddEntry("pion","#pi","L");
  if (CHARGE > 0 && (!skipPiMinus[yIndex])) leg.AddEntry("pion","#pi","L");
  if (CHARGE > 0 )  leg.AddEntry("deuteron","d","L");
  if (CHARGE > 0 )  leg.AddEntry("proton","p","L");


  if (!killElectron[yIndex]) leg.AddEntry("electron","e","L");
    leg.SetFillColor(kWhite);
  leg.SetBorderSize(0);
  leg.Draw("SAME");

  fittingCanvas->Update();
  
  gPad->Update();
 

    TString fullImageDir = IMAGEDIR+"/"+particleInfo->GetParticleName(PID,CHARGE);
    TString imageName = TString::Format("zTPCFit_%s_%02d_%02d_%02d",
					particleInfo->GetParticleName(PID,CHARGE).Data(),
					centralityIndex,yIndex,mTm0Index);

    fittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));

    return true;

} //End Fxn


//--------------------------------------------------------------------
Double_t GetYFromTGraph(TGraphErrors* graph, Int_t mtm0Index){

	Bool_t badmTm0Index = true;

	Int_t nPoints = graph->GetN();
	if (nPoints <= 0){
		cout<<"WARNING: GetYFromTGraph -- TGraph has no points"<<endl;
		return 0.;
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
Double_t GetYFromMeanTGraph(TGraphErrors* graph, Int_t mtm0Index){

  Bool_t badmTm0Index = true;

  Int_t nPoints = graph->GetN();
  if (nPoints <= 0){
    cout<<"WARNING: GetYFromTGraph -- TGraph has no points"<<endl;
    return -999;
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
      return -999;
    }
  }

}

//--------------------------------------------------------------------
void loadWidths(TString protonWidth, TString deuteronWidth){

  TFile *proFile = new TFile(protonWidth,"READ");


  proFile->cd();
  
  for (int iRapidityBin =firstRapidityIndex; iRapidityBin < lastRapidityIndex; iRapidityBin++){
    if(gDirectory->GetListOfKeys()->Contains(Form("matchedWidthFit_yIndex_%02d",iRapidityBin)))
      {
	cout << Form("Width file loading... Proton:::matchedWidthFit_yIndex_%02d",iRapidityBin) << endl;
	protonWidthFit.at(iRapidityBin) = (TF1 *)proFile->Get(Form("matchedWidthFit_yIndex_%02d",iRapidityBin));
	protonMeanGraph.at(iRapidityBin) = (TGraphErrors *)proFile->Get(Form("ProtonMean_%02d",iRapidityBin));

      }
    if ( protonWidthFit.at(iRapidityBin)) cout << "Rapidity bin " << iRapidityBin << " loaded... " << endl;

  }

  proFile->Close();


  TFile *deuFile = new TFile(deuteronWidth,"READ");

  deuFile->cd();

  for (int iRapidityBin =firstRapidityIndex; iRapidityBin < lastRapidityIndex; iRapidityBin++){
    if(gDirectory->GetListOfKeys()->Contains(Form("matchedWidthFit_yIndex_%02d",iRapidityBin)))
      {
	cout << Form("Width file loading... Deuteron:::matchedWidthFit_yIndex_%02d",iRapidityBin) << endl;
        deuteronWidthFit.at(iRapidityBin) = (TF1 *)deuFile->Get(Form("matchedWidthFit_yIndex_%02d",iRapidityBin));
	deuteronMeanGraph.at(iRapidityBin) = (TGraphErrors *)deuFile->Get(Form("DeuteronMean_%02d",iRapidityBin));
      }
  }

  deuFile->Close();


  cout << "Files loaded" << endl;

}

//--------------------------------------------------------------------
Float_t LowStatAmpSeed(TH1D * yieldHisto){

 
  int NBins = yieldHisto->GetNbinsX();
  float highBin =0;

  for (int i = 0; i < NBins; i++){

    float binCounts = (float) yieldHisto->GetBinContent(i);
    if (binCounts > highBin) highBin = binCounts;
  }

  return highBin;

}






