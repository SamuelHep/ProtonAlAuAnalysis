#include <iostream>
#include <TH1D.h>
#include <TFile.h>

void CombineCorrectedRefMultHistos(TFile *file, Int_t trigIndex, Int_t zIndexMin, Int_t zIndexMax){

  file->cd("RefMultHistos_Corrected");

  const Int_t nHistos = abs(zIndexMax-zIndexMin)+1;
  
  TH1D *refMultHistos[nHistos];
  for (Int_t i=0; i<nHistos; i++){
    refMultHistos[i] = (TH1D *)gDirectory->Get(Form("refMult_%d_%d_Corrected",trigIndex,i+zIndexMin));
    refMultHistos[i]->Sumw2();
  }
  
  TH1D *combinedRefMultHisto =
    new TH1D(Form("CorrectedRefMult_%d",trigIndex),"",
	     refMultHistos[0]->GetNbinsX(),
	     refMultHistos[0]->GetBinLowEdge(1),
	     refMultHistos[0]->GetBinLowEdge(refMultHistos[0]->GetNbinsX())+refMultHistos[0]->GetBinWidth(1));
  
  combinedRefMultHisto->Sumw2();

  
  for (Int_t i=0; i<nHistos; i++){
    combinedRefMultHisto->Add(refMultHistos[i]);
  }

  combinedRefMultHisto->Draw();

  TFile *outFile = new TFile("AuAu07_CombinedCorrectedRefMult.root","RECREATE");
  combinedRefMultHisto->Write();
  outFile->Close();
  

}
