#include "../inc/globalDefinitions.h"

void DrawAllSpectra(TString spectraFile, TString imageDir, TString system, Double_t energy, Int_t speciesIndex, 
			Int_t charge, Bool_t corrected=true){

  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gStyle->SetOptFit(1);

  ParticleInfo *particleInfo = new ParticleInfo();
  
  //Load the root file
  TFile *file = new TFile(spectraFile,"READ");

  //Create the Spectra Name
  TString type = "raw";
  TString Type = "Raw";
  TString name = "Uncorrected";

  TString speciesName = "ProtonPlus";

  if (corrected){
    type = "corrected";
    Type = "Corrected";
    name = "Corrected";
  }

  for (int centIndex=0;centIndex<1;centIndex++){
    for (int yIndex=7;yIndex<19;yIndex++){

      TString spectraName(Form("%sSpectra_ProtonPlus_Cent%02d_yIndex%02d_Total",
                               Type.Data(),centIndex,yIndex));
      
      TString spectraFitName(Form("%sSpectra_ProtonPlus_Cent%02d_yIndex%02d_Total_Fit1",
                                  Type.Data(),centIndex,yIndex));

  //Get the spectrum
  TGraphErrors *spectrum = (TGraphErrors *)file->Get(Form("%sSpectra_ProtonPlus/%s",
                                                       Type.Data(),
                                                       spectraName.Data()));

  spectrum->SetMarkerColor(particleInfo->GetParticleColor(speciesIndex));



  //Create the Drawing Canvas
  TCanvas *canvas = new TCanvas("spectrumCanvas","Spectrum Canvas",20,20,800,600);
  canvas->SetTopMargin(.05);
  canvas->SetRightMargin(.05);
  canvas->SetLeftMargin(.12);
  canvas->SetTicks();
  canvas->SetLogy();


  //Determine what the axis ranges should be and then draw the frame
  Double_t maxY = TMath::MaxElement(spectrum->GetN(),spectrum->GetY())*1.6;
  Double_t minY = TMath::MinElement(spectrum->GetN(),spectrum->GetY())*.5;
  Double_t maxX = TMath::MaxElement(spectrum->GetN(),spectrum->GetX())+5*mTm0BinWidth;
  Double_t minX = 0;
  
  TH1F *frame = canvas->DrawFrame(minX,minY,maxX,maxY);
  frame->SetTitle(Form(";m_{T}-m_{%s} (GeV/c^{2});#frac{1}{N_{Evt}}#times#frac{1}{2#pim_{T}}#times#frac{d^{2}N}{dm_{T}dy} (GeV/c^{2})^{2}",
		       particleInfo->GetParticleSymbol(speciesIndex,charge).Data()));
  frame->GetYaxis()->SetTitleOffset(1.3);
  
  spectrum->Draw("PZ");


  //Create the Title
  TPaveText *title = new TPaveText(.52,.81,.9,.9,"brNDC");
  title->SetTextAlign(12);
  title->SetBorderSize(0);
  title->SetFillColor(kWhite);
  title->AddText(Form("%s Spectrum of %s",
		      name.Data(),particleInfo->GetParticleSymbol(speciesIndex,charge).Data()));
  title->GetLine(0)->SetTextSize(.038);
  title->AddText(Form("%s at #sqrt{s_{NN}} = %g GeV",
		      system.Data(),energy));
  title->GetLine(1)->SetTextSize(.03);
  title->Draw("SAME");

  gPad->Update();
  TPaveStats *ps = (TPaveStats*)canvas->GetPrimitive("stats");
  ps->SetX1NDC(.17);
  ps->SetY1NDC(.17);
  ps->SetX2NDC(.45);
  ps->SetY2NDC(.45);

  //Create the Legend
  TLegend *leg = new TLegend(.55,.76,.7,.81);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(.035);
  leg->AddEntry(spectrum,Form("y = %.03g",GetRapidityRangeCenter(yIndex)),"P");
  leg->Draw("SAME");


  TString fullImageDir = TString::Format("%s",imageDir.Data());
  TString imageName = TString::Format("Spetrum_%s_%02d_%02d",
				      particleInfo->GetParticleName(2,charge).Data(),
				      centIndex,yIndex);

  canvas->SaveAs(Form("%s/%s.eps",fullImageDir.Data(),imageName.Data()));
  delete canvas;
    }
  }
}
