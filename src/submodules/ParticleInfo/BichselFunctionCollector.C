/*
  Saves the values of the Bichsel curve as a function of beta*gamma to a root file
  Both the curve using the truncated mean and most probable dEdx from the fit are
  stored. 
*/

#include <iostream>
#include <fstream>

#include <math.h>

using namespace std;

//________________________________________________________________________________
void loadLib(){
  if (gClassTable->GetID("Bichsel") < 0) {
    gSystem->Load("St_base");
    gSystem->Load("StBichsel");
  }
}

//________________________________________________________________________________
Double_t bichselFunction(Double_t *x, Double_t *par){

  //The Bichsel curve itself is a parameterization of energy loss as a
  //function of log (base 10) of beta*gamma. (See StBichsel class).
  Bichsel *bichsel = Bichsel::Instance();

  //If par[0] is is zero then return the value obtained from the truncated mean,
  if ((Int_t)par[0] == 0)
    return bichsel->GetI70(log10(x[0]));

  //Otherwise return the value obtained from the most probable fit
  return TMath::Exp(bichsel->GetMostProbableZ(log10(x[0]),1));

}

//____MAIN_______________________________________________________________________
void BichselFunctionCollector(TString fileNameBase){

  //Set the min and max betta*gamma and the the npx values
  Double_t minBetaGamma = 0.1;
  Double_t midBetaGamma = 10.0;
  Double_t maxBetaGamma(10000);
  Int_t npx(100000);

  //Load the Necessary Libraries
  loadLib();

  //Create the Bichsel function for the truncated mean
  TF1 *bichselTruncated_0 = new TF1("bichselTruncated_0",bichselFunction,minBetaGamma,midBetaGamma,1);
  bichselTruncated_0->SetNpx(npx);
  bichselTruncated_0->SetParameter(0,0); //Setting par[0] to 0 means use the truncated method

  TF1 *bichselTruncated_1 = new TF1("bichselTruncated_1",bichselFunction,midBetaGamma,maxBetaGamma,1);
  bichselTruncated_1->SetNpx(npx);
  bichselTruncated_1->SetParameter(0,0); //Setting par[0] to 0 means use the truncated method


  //Create the Bichsel function for the truncated mean
  TF1 *bichselMostProbable_0 = new TF1("bichselMostProbable_0",bichselFunction,minBetaGamma,midBetaGamma,1);
  bichselMostProbable_0->SetNpx(npx);
  bichselMostProbable_0->SetParameter(0,1); //Setting par[0] to 1 means use the most probable method

  TF1 *bichselMostProbable_1 = new TF1("bichselMostProbable_1",bichselFunction,midBetaGamma,maxBetaGamma,1);
  bichselMostProbable_1->SetNpx(npx);
  bichselMostProbable_1->SetParameter(0,1);

  //Create File to Write the Curves to and write them
  TFile *bichselCurveFile = new TFile(Form("PidFunctions_Data/%s_BichselCurves.root",fileNameBase.Data()),"CREATE");
  bichselTruncated_0->Write();
  bichselTruncated_1->Write();
  bichselMostProbable_0->Write();
  bichselMostProbable_1->Write();

  bichselCurveFile->Close();
}
