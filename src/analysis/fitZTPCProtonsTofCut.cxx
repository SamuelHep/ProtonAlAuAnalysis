//Extracts the yield of protons
//Currently assumes amplitude of pi+ yield to be the same as the amplitude of the pi- yield
//Currently does not account for deuteron contamination

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <sys/stat.h>
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
std::vector<std::vector<TH1D *> > tpcProtonPlusMatchedHisto(nRapidityBins,std::vector<TH1D *>(nmTm0Bins,(TH1D *)NULL));
std::vector<std::vector<TF1 *> > tpcProtonPlusMatchedFit(nRapidityBins,std::vector<TF1 *>(nmTm0Bins,(TF1 *)NULL));
std::vector<TF1 *>  tpcProtonWidthFit(nRapidityBins,(TF1 *)NULL);

std::vector<TGraphErrors *>  tpcProtonMeanGraph(nRapidityBins,(TGraphErrors *)NULL);

TCanvas * matchedFittingCanvas;
TCanvas * matchedParCanvas;

//Particle predictions
ParticleInfo *particleInfo;
ParticlePars *proton;
Double_t protonPredict;
TLine * protonLine;

//global vars 
Double_t rapidity, mTm0, mT, pZ, pT, pTotal;
Bool_t DRAW = true;

//Functions
void FitProtons(TString imageDirect);
void SaveFitWidths(TString rootFile, TString imageDirect);
void ComputeBinKinematics(Int_t mTm0Index, Double_t y, Int_t PID);

//-------MAIN-------
void fitZTPCProtonsTofCut(TString ZTPCFILE,TString OUTROOT, TString STARLIB, TString IMAGEDIR){ 

  //Get the Yield Histogram file and Get the 3D Z Histograms 
  TFile *inFile = new TFile(ZTPCFILE,"READ");
  TH3D *protonPlusZTPC3D = (TH3D *)inFile->Get(Form("ProtonPlus/zTPC_ProtonPlus_con_ProtonPlus_AllCent_tofCut"));

  LoadYieldHistograms(protonPlusZTPC3D,&tpcProtonPlusMatchedHisto,"TPC",-1);
     
  //Make Proton prediction... should be zero
  particleInfo = new ParticleInfo(STARLIB,true);
  proton = new ParticlePars();
  DefineParticlePars(proton,"Proton");
  protonLine = new TLine();
  protonLine->SetLineColor(particleInfo->GetParticleColor(PROTON));
  protonLine->SetLineWidth(3);

  matchedFittingCanvas = new TCanvas("matchedFittingCanvas","matchedFittingCanvas",800,600);
  matchedFittingCanvas->SetLogy();

  matchedParCanvas  = new TCanvas("matchedParCanvas","matchedParCanvas",900,300);
  matchedParCanvas->Divide(3,1); 

  //Set Style Options
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(112);
  gStyle->SetTitleFontSize(0.1);

  //Fit the Proton Pluses 
  FitProtons(IMAGEDIR);
  SaveFitWidths(OUTROOT,IMAGEDIR);
	
}

void FitProtons(TString imageDirect){
  
  //Local Variables  
  TH1D *yieldHisto = NULL;
  TF1 *yieldFit = NULL;
  TF1 *widthFit = NULL;
  Int_t chargePID = 1;
  Double_t min = -0.5;
  Double_t max = 0.5;

  for (Int_t yIndex=0; yIndex<19; yIndex++){
    
    for (Int_t mTm0Index=0; mTm0Index<nmTm0Bins; mTm0Index++){
      
      //Make sure the histogram exists
      if (!tpcProtonPlusMatchedHisto.at(yIndex).at(mTm0Index))
	continue;
      
      //Pass the histogram to a local pointer
      yieldHisto = tpcProtonPlusMatchedHisto.at(yIndex).at(mTm0Index);
      yieldHisto->SetAxisRange(-2.,2.5);
      
      //continue if too few entries 
      if (yieldHisto->GetEntries() < 100) continue;

      //Compute bin kinematics
      rapidity = GetRapidityRangeCenter(yIndex);
      ComputeBinKinematics(mTm0Index,rapidity,PROTON);
      protonPredict = particleInfo->PredictZTPC(pTotal,PROTON,PROTON);

      //Only to 1 GeV
      if (mTm0 > 1)
        break;

      matchedFittingCanvas->cd();
      yieldHisto->Draw("E");
      yieldHisto->SetStats(0);
      gPad->Update();
      
      tpcProtonPlusMatchedFit.at(yIndex).at(mTm0Index) = new TF1("","gaus(0)",min,max);
      yieldFit = tpcProtonPlusMatchedFit.at(yIndex).at(mTm0Index);

      yieldFit->SetParName(0,"P_{Amp}");
      yieldFit->SetParName(1,"#mu_{P}");
      yieldFit->SetParName(2,"{#sigma}_{P}");

      yieldFit->SetParameter(0,yieldHisto->GetBinContent(yieldHisto->FindBin(protonPredict)));
      yieldFit->SetParameter(1,protonPredict);
      yieldFit->SetParameter(2,.1);

      //Fit proton
      yieldHisto->Fit(yieldFit,"R");

      //Add fit results to parameter graphs 
      AddPointToParGraph(proton->mean[yIndex],yieldFit,1,mTm0Index);
      AddPointToParGraph(proton->width[yIndex],yieldFit,2,mTm0Index);
      
      //Draw parameter graphs
      matchedParCanvas->cd(1);
      proton->mean[yIndex]->GetXaxis()->SetRangeUser(0,1);
      proton->mean[yIndex]->GetXaxis()->SetTitle("mTm0");
      proton->mean[yIndex]->GetYaxis()->SetRangeUser(-0.06,0.06);
      proton->mean[yIndex]->GetYaxis()->SetTitle("ZTPC mean");
      proton->mean[yIndex]->GetYaxis()->SetTitleOffset(1.5);
      proton->mean[yIndex]->Draw("APZ");
      proton->mean[yIndex]->SetTitle(Form("ProtonMean_%02d",yIndex));
      proton->mean[yIndex]->SetName(Form("ProtonMean_%02d",yIndex));
      
      matchedParCanvas->cd(2);
      proton->width[yIndex]->GetXaxis()->SetRangeUser(0,1);
      proton->width[yIndex]->GetXaxis()->SetTitle("mTm0");
      proton->width[yIndex]->GetYaxis()->SetRangeUser(0.0,0.12);
      proton->width[yIndex]->GetYaxis()->SetTitle("ZTPC width");
      proton->width[yIndex]->GetYaxis()->SetTitleOffset(1.5);
      proton->width[yIndex]->Draw("APZ");
      matchedParCanvas->Update();

      tpcProtonMeanGraph.at(yIndex) = proton->mean[yIndex];
      
    }//End Loop Over mTm0Index
    
    if (DRAW){
      TString fullImageDir = imageDirect;
      TString imageName = TString::Format("ztpc_tofCut_Proton_Widths_Amp_y%02d",yIndex);
      
      matchedParCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
    }
    
    tpcProtonWidthFit.at(yIndex) = new TF1("linearFit","[0]",0,1.0);
    tpcProtonWidthFit.at(yIndex)->SetParameter(0,-999);
    
    if (proton->width[yIndex]->GetN() > 3){
      
      widthFit = tpcProtonWidthFit.at(yIndex);
      widthFit->SetName(Form("matchedWidthFit_yIndex_%02d",yIndex));
      widthFit->SetTitle(Form("matchedWidthFit_yIndex_%02d",yIndex));
      widthFit->SetParameter(0,0.07);
      widthFit->SetParLimits(0,0.02,0.12);
      
      proton->width[yIndex]->Fit(widthFit,"R");
      
      if (DRAW){	  
	matchedParCanvas->Update();
	
	TString fullImageDir = imageDirect;
	TString imageName = TString::Format("zTPC_tofCut_ProtonParameters_y%02d",yIndex);
	
	matchedParCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
      }
      
    }
  }  //End Loop Over yIndex
} //End Fxn


//Saves the Width fit and the Draws Histogram
void SaveFitWidths(TString rootFile, TString imageDirect){

  //setup the output file
  TFile * outroot = new TFile(rootFile,"RECREATE");
  
  Double_t width = -999;
  Double_t widthErr = -999;
  TH1D *yieldHisto = NULL;
  TF1 *yieldFit = NULL;

  for (Int_t yIndex=0; yIndex<19; yIndex++){
    for (Int_t mTm0Index=0; mTm0Index<nmTm0Bins; mTm0Index++){
  
      width =-999;
      widthErr = -999;
      yieldHisto = tpcProtonPlusMatchedHisto.at(yIndex).at(mTm0Index);
      yieldFit = tpcProtonPlusMatchedFit.at(yIndex).at(mTm0Index);

      //Make sure both the histogram and fit exist      
      if (!yieldHisto || !yieldFit) continue;      

      //Draw Fit
      matchedFittingCanvas->Clear();
      matchedFittingCanvas->cd();
      matchedFittingCanvas->SetTopMargin(.05);
      matchedFittingCanvas->SetRightMargin(.05);
      yieldHisto->SetMarkerStyle(kOpenCircle);
      yieldHisto->SetTitle("");
      yieldHisto->SetAxisRange(-2.0,2.0);
      yieldHisto->Draw("E");
      
      yieldFit->SetLineWidth(3);
      yieldFit->SetLineColor(particleInfo->GetParticleColor(PROTON));
      yieldFit->SetLineStyle(3);
      yieldFit->SetNpx(10000);
      yieldFit->Draw("SAME");

      width  = yieldFit->GetParameter(1);
      widthErr = yieldFit->GetParError(1);

      if (DRAW){
	TPaveText label(.12,.8,.35,.92,"BL NDC");
	label.SetTextAlign(11);
	label.AddText(Form("Z(%s) Yield Extraction Fit",particleInfo->GetParticleSymbol(PROTON,1).Data()));
	label.AddText(Form("  %.2f #leq y_{%s} < %.2f | %.2f #leq m_{T}-m_{%s} < %.2f",
			   GetRapidityRangeLow(yIndex),particleInfo->GetParticleSymbol(PROTON,1).Data(),
			   GetRapidityRangeHigh(yIndex),GetmTm0RangeLow(mTm0Index),
			   particleInfo->GetParticleSymbol(PROTON,1).Data(),GetmTm0RangeHigh(mTm0Index)));
	label.GetLine(0)->SetTextSize(.03);
	label.GetLine(1)->SetTextSize(.02);
	label.SetFillColor(kWhite);
	label.SetBorderSize(0);
	label.Draw("SAME");

	TString fullImageDir = imageDirect;
	TString imageName = TString::Format("zTPC_tofCut_ProtonPlus_y%02d_mt%02d",yIndex,mTm0Index);
	matchedFittingCanvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
	
      }
    }
  
    outroot->cd();
    if(tpcProtonWidthFit.at(yIndex)->GetParameter(0) > 0) {
      tpcProtonWidthFit.at(yIndex)->Write();
      tpcProtonMeanGraph.at(yIndex)->Write();
    }
  }
  outroot->Close();
}


void ComputeBinKinematics(Int_t mTm0Index, Double_t y, Int_t PID){

  mTm0 = GetmTm0RangeCenter(mTm0Index);
  mT   = mTm0 + particleInfo->GetParticleMass(PID);
  pZ   = ComputepZ(mT,y);
  pT   = ConvertmTm0ToPt(mTm0,particleInfo->GetParticleMass(PID));
  pTotal = ComputepTotal(pT,pZ);

}


