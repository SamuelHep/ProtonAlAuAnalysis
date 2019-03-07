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
#include <TMath.h>

#include "Rtypes.h"
#include "globalDefinitions.h"
#include "ParticleInfo.h"
#include "TFitResult.h"
#include "TrackInfo.h"
#include "PrimaryVertexInfo.h"
#include "EventInfo.h"
#include "utilityFunctions.h"
#include "fitZTPCUtilities.h"



Double_t GetXatMaxY(TGraphErrors* graph);
TGraphErrors * reflectCombineTGraph(TGraphErrors * graph, TString graphName);

void fitdNdy(TString BLASTWAVEFILE, TString OUTFILE, TString IMAGEDIR){

  TFile * inFile = new TFile(BLASTWAVEFILE,"UPDATE");

  TFile * outFile = new TFile(OUTFILE,"RECREATE");

  TGraphErrors * dNdyGraphInput = (TGraphErrors*) inFile->Get("dNdyGraph_ProtonPlus/dNdyGraphStat_Cent00");

  //TGraphErrors * dNdyGraph = reflectCombineTGraph(dNdyGraphInput,"AuAuGraph");
  TGraphErrors * dNdyGraph = dNdyGraphInput;

  gStyle->SetOptFit(1);

  dNdyGraph->SetLineWidth(2);
  dNdyGraph->SetLineColor(kBlack);
  dNdyGraph->SetMarkerColor(kBlack);
  dNdyGraph->SetMarkerStyle(kFullCircle);
  dNdyGraph->SetMarkerSize(1.4);
  


  TCanvas * fittingCanvas = new TCanvas("fittingCanvas","fittingCanvas",800,600);
  fittingCanvas->SetTickx(1);
  fittingCanvas->SetTicky(1);

  TH1F * frame;
  frame = fittingCanvas->DrawFrame(-1.5,0,0,100);
  frame->GetXaxis()->SetLabelSize(.040);
  frame->GetXaxis()->SetTitle("y - y_{CM}");
  frame->GetYaxis()->SetLabelSize(.040);
  frame->GetYaxis()->SetTitle("dN/dy");
  frame->GetYaxis()->SetLabelOffset(0.01);
  frame->GetYaxis()->SetTitleOffset(1.5);

  dNdyGraph->Draw("SAME PZ");
  

  TF1 * doubleGaus = new TF1("doubleGaus","[0]*exp(-0.5*pow( ((x-[1]) /[2]) ,2)) + [3]*exp(-0.5*pow( ((x-[4])/[5]) ,2))",-2,0);

  doubleGaus->SetLineColor(kRed);
  /*
  doubleGaus->SetParameter(0,dNdyGraph->GetMaximum());
  doubleGaus->FixParameter(1,0);
  doubleGaus->SetParameter(2,.7);

  doubleGaus->SetParameter(3,dNdyGraph->GetMaximum());
  doubleGaus->SetParameter(4,GetXatMaxY(dNdyGraph));
  doubleGaus->SetParameter(5,0.5);
  */
  doubleGaus->SetParameter(0,30);
  doubleGaus->SetParameter(1,-1.2);
  //  doubleGaus->SetParLimits(1,0.1,0.6);
  doubleGaus->SetParameter(2,.3);

 doubleGaus->SetParameter(3,40);
 doubleGaus->SetParameter(4,-0.4);
  //  doubleGaus->SetParLimits(1,0.1,0.6);
  doubleGaus->SetParameter(5,.4);

  //  doubleGaus->SetParameter(3,30);
  //  doubleGaus->SetParameter(4,0.4);
  //  doubleGaus->SetParameter(5,0.5);

  doubleGaus->SetParName(0,"Amp. 1");
  doubleGaus->SetParName(1,"#mu 1");
  doubleGaus->SetParName(2,"#sigma 1");
  doubleGaus->SetParName(3,"Amp. 2");
  doubleGaus->SetParName(4,"#mu 2");
  doubleGaus->SetParName(5,"#sigma 2");


  //  doubleGaus->SetParName(3,"Stop Amp.");
  //  doubleGaus->SetParName(4,"Stop #mu");
  //  doubleGaus->SetParName(5,"Stop #sigma");

  doubleGaus->SetLineWidth(3);

  TFitResultPtr fitResult = dNdyGraph->Fit(doubleGaus,"S R");

  Double_t intDoubleGaus = doubleGaus->Integral(-2,0);

  Double_t intDoubleGausErr = doubleGaus->IntegralError(-2,2,fitResult->GetParams());

  cout << "4pi yield=" << intDoubleGaus << " +/- " << intDoubleGausErr << endl;

  doubleGaus->Draw("SAME");

  TF1 * PosFunc = new TF1("positiveGaus","gaus(0)",-2,0);
  TF1 * NegFunc = new TF1("negativeGaus","gaus(0)",-2,0);

  PosFunc->FixParameter(0,doubleGaus->GetParameter(0));
  PosFunc->FixParameter(1,doubleGaus->GetParameter(1));
  PosFunc->FixParameter(2,doubleGaus->GetParameter(2));

  PosFunc->SetParError(0,doubleGaus->GetParError(0));
  PosFunc->SetParError(1,doubleGaus->GetParError(1));
  PosFunc->SetParError(2,doubleGaus->GetParError(2));

  PosFunc->SetLineColor(kBlue);
  PosFunc->SetLineWidth(3);
  PosFunc->SetLineStyle(2);
  PosFunc->SetNpx(1000);
  PosFunc->Draw("SAME");

  NegFunc->FixParameter(0,doubleGaus->GetParameter(3));
  NegFunc->FixParameter(1,doubleGaus->GetParameter(4));
  NegFunc->FixParameter(2,doubleGaus->GetParameter(5));

  NegFunc->SetParError(0,doubleGaus->GetParError(3));
  NegFunc->SetParError(1,doubleGaus->GetParError(4));
  NegFunc->SetParError(2,doubleGaus->GetParError(5));

  NegFunc->SetLineColor(kOrange);
  NegFunc->SetLineWidth(3);
  NegFunc->SetLineStyle(2);
  NegFunc->SetNpx(1000);
  NegFunc->Draw("SAME");

    
  TLegend * leg = new TLegend(.15,.72,.40,.83);
  leg->AddEntry(PosFunc, "Gaus 1","L");
  leg->AddEntry(NegFunc, "Gaus 2","L");
  leg->AddEntry(doubleGaus, "1 + 2","L");
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->Draw("SAME");

  gPad->Update();

   fittingCanvas->Update();
   fittingCanvas->SaveAs(Form("%s/dNdy_2GausFit.pdf",IMAGEDIR.Data()));


   outFile->cd();
   PosFunc->Write();
   NegFunc->Write();
   doubleGaus->Write();

  
}


Double_t GetXatMaxY(TGraphErrors* graph){

  Double_t max = graph->GetMaximum();
  Double_t maxX = -999;
  Double_t dymax = -999;

	Int_t nPoints = graph->GetN();
	if (nPoints <= 0){
		cout<<"WARNING: GetXatMaxY -- TGraph has no points"<<endl;
		return 0.0;
	}
	Double_t *xPointsArr = graph->GetX();
	Double_t *yPointsArr = graph->GetY();
	for (Int_t iPoint = 0; iPoint<nPoints; iPoint++){
	  
	  Float_t xVal = xPointsArr[iPoint];
	  Float_t yVal = yPointsArr[iPoint];
	  
	  Float_t dyVal = yVal-max;
	    
	    if ( dyVal > dymax){
	      dymax = dyVal;
	      maxX = xVal;
	    }
	  
	}

	return maxX;

}

TGraphErrors * reflectCombineTGraph(TGraphErrors * graph, TString graphName){

	Int_t nPoints = graph->GetN();
	if (nPoints <= 0){
		cout<<"WARNING: Combined and Reflect-- TGraph has no points"<<endl;
	}

	TGraphErrors * outGraph = new TGraphErrors();

	Double_t *xPointsArr = graph->GetX();
	Double_t *yPointsArr = graph->GetY();
	
	Double_t *xPointsErArr = graph->GetEX();
	Double_t *yPointsErArr = graph->GetEY();

	for (Int_t iPoint = nPoints-1; iPoint>=0; iPoint--){

	  outGraph->SetPoint(outGraph->GetN(),-xPointsArr[iPoint],yPointsArr[iPoint]);
	  //	  outGraph->SetPointError(outGraph->GetN()-1,-xPointsErArr[iPoint],yPointsErArr[iPoint]);
	  outGraph->SetPointError(outGraph->GetN()-1,-xPointsErArr[iPoint],1);

	}

	for (Int_t iPoint = 0; iPoint<nPoints; iPoint++){

	  outGraph->SetPoint(outGraph->GetN(),xPointsArr[iPoint],yPointsArr[iPoint]);
	  //	  outGraph->SetPointError(outGraph->GetN()-1,xPointsErArr[iPoint],yPointsErArr[iPoint]);
	  outGraph->SetPointError(outGraph->GetN()-1,xPointsErArr[iPoint],0);
	  
	}

	outGraph->SetName(Form("%s",graphName.Data()));

	return outGraph;
}


