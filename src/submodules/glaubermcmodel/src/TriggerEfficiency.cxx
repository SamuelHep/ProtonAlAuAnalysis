#include <iostream>

#include <TFile.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveStats.h>

//______________________________________________________________________________________________
TGraphErrors *triggerEfficiency(TH1D *bestFit,TH1D *refMult,
				TString graphName,TString graphTitle){

  //Return pointer to graph of the ratio of the simulated GMC+NB and
  //corrected observed multiplicity distributions. The graph has been fit.

  TGraphErrors *gr = new TGraphErrors();
  gr->SetMarkerStyle(kFullCircle);
  gr->SetName(graphName);
  gr->SetTitle(graphTitle);
  
  for (Int_t i=1; i<refMult->GetNbinsX(); i++){
    
    Double_t refMultVal = refMult->GetBinContent(i);
    Double_t bestFitVal = bestFit->GetBinContent(i);
    Double_t refMultErr = refMult->GetBinError(i);
    Double_t bestFitErr = bestFit->GetBinError(i);
    
    if (refMultVal == 0 || bestFitVal == 0)
      continue;
    
    Double_t ratio = refMultVal/bestFitVal;
    Double_t ratioErr = ratio * sqrt(pow(bestFitErr/bestFitVal,2) + pow(refMultErr/refMultVal,2));
    
    gr->SetPoint(gr->GetN(),refMult->GetBinCenter(i),ratio);
    gr->SetPointError(gr->GetN()-1,0,ratioErr);
    
  }
  
  TF1 *fit = new TF1("","[0]+[1]/x",0,550);
  fit->SetNpx(10000);
  fit->SetParameter(0,1);
  fit->SetParLimits(0,.9,1.1);
  fit->SetParameter(1,.1);
  fit->SetLineColor(kRed);
  gr->Fit(fit,"","",5,100);
  fit->SetRange(0,550);
  
  
  return gr;
  
}

//___________________________________________________________________________________________
void triggerEfficiency_TEST(){

  //Wrapper for the above function useful for block tests

  TFile *file = new TFile("../data/outfile1.root");
  TH1D *refMult = (TH1D *)file->Get("CorrectedRefMult_0");
  TH1D *bestFit = (TH1D *)file->Get("bestFit");
  
  TCanvas *canvas = new TCanvas("canvas","canvas",20,20,800,600);
  canvas->SetLogy();

  refMult->Draw("E");
  bestFit->Draw("SAME");

  TGraphErrors *gr = triggerEfficiency(bestFit,refMult,
				       Form("%s_%s",refMult->GetName(),"_TrigEff"),
				       refMult->GetTitle());
				       
  gStyle->SetOptFit(1);
  TCanvas *canvas1 = new TCanvas("canvas1","canvas1",0,0,800,600);
  canvas1->SetTopMargin(.05);
  canvas1->SetRightMargin(.05);
  TH1F *frame = canvas1->DrawFrame(0,-0.76,1.1*gr->GetX()[gr->GetN()-1],2.1);
  frame->SetTitle(";Corrected RefMult;GMC+NB/Data");
  
  gr->Draw("P");
  canvas1->Update();
  
  TPaveStats *statsBox = (TPaveStats *)gr->FindObject("stats");
  statsBox->SetX1NDC(.15);
  statsBox->SetY1NDC(.15);
  statsBox->SetX2NDC(.5);
  statsBox->SetY2NDC(.31);
  
  canvas1->Modified();
  canvas1->Update();
  
}
