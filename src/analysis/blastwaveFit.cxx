//Sim Fit the proton pion and lambda spectra
#include <iostream>
#include <vector>
#include <utility>

#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TVirtualFitter.h>
#include <TAxis.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TLine.h>

#include "globalDefinitions.h"
#include "utilityFunctions.h"
#include "StRefMultExtendedCorr.h"
#include "UserCuts.h"
#include "ParticleInfo.h"
#include "SpectraClass.h"
#include "SpectraFitUtilities.h"
#include "SpectraFitFunctions.h"
#include "SimFitter.h"
#include "StyleSettings.h"

void blastwaveFit(){

  ParticleInfo * particleInfo = new ParticleInfo();

  gStyle->SetOptFit(0);  

  //Set charge vector
  std::vector<int> charges;
  charges.push_back(-1);
  charges.push_back(1);
  charges.push_back(1);
  charges.push_back(0);

  std::vector<int> particles;
  particles.push_back(PION);
  particles.push_back(PROTON);
  particles.push_back(KAON);
  particles.push_back(LAMBDA);

  const unsigned int nParticles(particles.size());

  TF1 *heinzBlastWaveFunc = new TF1("heinzBlastWaveFunc",BlastWaveModelFit,0.001,2,6);
  heinzBlastWaveFunc->SetNpx(100);

  TCanvas * canvas = new TCanvas("fittingCanvas","fittingCanvas",1200,800);
  canvas->Divide(2,2);

  //  TF1 *boseBlastWaveFunc = new TF1("boseBlastWaveFunc",BoseBlastWaveFit,0.001,2,8);
  //  boseBlastWaveFunc->SetNpx(100);

  std::vector<TGraphErrors *> spectraToFit;
  std::vector<TF1 *> fitFuncs;

  TFile *protonFile = new TFile("/home/sheppelmann/fxt2015/davisdstanalysis/userfiles/AuAu_4_5GeV_2015/analysis/zTPCProtonSpectra.root","READ");
  TFile *kaonFile = new TFile("../userfiles/comparisonData/kaon/kaonSpectra.root","READ");
  TFile *pionFile = new TFile("../userfiles/comparisonData/pion/pionSpectra.root","READ");
  TFile *lambdaFile = new TFile("../userfiles/comparisonData/lambda/lambdaSpectra.root","READ");

  TGraphErrors * protonSpectra = (TGraphErrors*) protonFile->Get("CorrectedSpectra_ProtonPlus/correctedSpectra_ProtonPlus_Cent00_yIndex07");
  TGraphErrors * pionSpectra = (TGraphErrors*) pionFile->Get("CorrectedSpectra_PionMinus/correctedSpectra_PionMinus_Cent00_yIndex05");
  TGraphErrors * kaonSpectra = (TGraphErrors*) kaonFile->Get("kaonSpectra");
  TGraphErrors * lambdaSpectra = (TGraphErrors*) lambdaFile->Get("lambdaSpectra");
  
  Int_t EXIT = 0;
  if (!protonSpectra){ cout << "No proton spectra graph found" << endl; EXIT = 1; }
  if (!pionSpectra){ cout << "No pion spectra graph found" << endl; EXIT = 1; }
  if (!kaonSpectra){ cout << "No kaon spectra graph found" << endl; EXIT = 1; }
  if (!lambdaSpectra){ cout << "No lambda spectra graph found" << endl; EXIT = 1; }

  if (EXIT) return;
    
  for (unsigned int iParticle=0;iParticle<nParticles;iParticle++){

    TString partSym = particleInfo->GetParticleSymbol(particles.at(iParticle));

     fitFuncs.push_back(new TF1(*heinzBlastWaveFunc));
     fitFuncs.at(iParticle)->SetName(Form("%s_blastwaveFit", particleInfo->GetParticleName(particles.at(iParticle),charges.at(iParticle)).Data()));
     fitFuncs.at(iParticle)->SetLineColor(particleInfo->GetParticleColor(particles.at(iParticle)));
     fitFuncs.at(iParticle)->SetLineStyle(charges.at(iParticle)<0 ? 1 : 9 );
     fitFuncs.at(iParticle)->SetLineWidth(2);
     fitFuncs.at(iParticle)->SetParName(0,Form("A_{%s}",partSym.Data()) );
     fitFuncs.at(iParticle)->SetParName(1,"T_{Kin}");
     fitFuncs.at(iParticle)->SetParName(2,"#beta_{s}");
     fitFuncs.at(iParticle)->SetParName(3,"n");
     fitFuncs.at(iParticle)->SetParName(4,"R");
     fitFuncs.at(iParticle)->SetParName(5,Form("m_{%s}",partSym.Data() ) );

     //Set Particle Mass
     fitFuncs.at(iParticle)->FixParameter(5,particleInfo->GetParticleMass(particles.at(iParticle)));	     
     
    if (particles.at(iParticle) == PROTON) spectraToFit.push_back(new TGraphErrors(*protonSpectra));
    if (particles.at(iParticle) == PION) spectraToFit.push_back(new TGraphErrors(*pionSpectra));
    if (particles.at(iParticle) == KAON) spectraToFit.push_back(new TGraphErrors(*kaonSpectra));
    if (particles.at(iParticle) == LAMBDA) spectraToFit.push_back(new TGraphErrors(*lambdaSpectra));

  }

 enum blastWavePars {BW_TEMP,BW_SURFACEVELOCITY,BW_TRANSVELOCITYPOWER,BW_RADIALINTEGRALLIMIT,
		     BW_PIONMASS,BW_KAONMASS,BW_LAMBDAMASS,BW_PROTONMASS,BW_AMP};//,
 
 int nonAmpPars = 8;
 const int nTotalPars = nonAmpPars + nParticles;

 std::vector<std::vector<int> > parRelations(nParticles,std::vector<int> (0,-1));
 for (unsigned int i=0; i < nParticles; i++){

   parRelations.at(i).resize(fitFuncs.at(i)->GetNpar(),0);
   parRelations.at(i).at(0) = BW_AMP + i;
   parRelations.at(i).at(1) = BW_TEMP;
   parRelations.at(i).at(2) = BW_SURFACEVELOCITY;
   parRelations.at(i).at(3) = BW_TRANSVELOCITYPOWER;
   parRelations.at(i).at(4) = BW_RADIALINTEGRALLIMIT;

 }

 std::vector<double> totalFitPars(nTotalPars);
 
 //Seed values for blastwave
 totalFitPars.at(BW_TEMP) = .13;
 totalFitPars.at(BW_SURFACEVELOCITY) = .64;
 totalFitPars.at(BW_TRANSVELOCITYPOWER) = .73;
 totalFitPars.at(BW_RADIALINTEGRALLIMIT) = 1;

 //Masses
 totalFitPars.at(BW_PIONMASS) = particleInfo->GetParticleMass(PION);
 totalFitPars.at(BW_KAONMASS) = particleInfo->GetParticleMass(KAON);
 totalFitPars.at(BW_LAMBDAMASS) = particleInfo->GetParticleMass(LAMBDA);
 totalFitPars.at(BW_PROTONMASS) = particleInfo->GetParticleMass(PROTON);

 //Amplitudes
 for (unsigned int i=0; i<nParticles; i++){
   totalFitPars.at(BW_AMP+i) = 10000;
 }

 //Set fit ranges
 std::vector<std::pair<float,float> > range;

 range.push_back( make_pair( .1, .5) );
 range.push_back( make_pair( .2, 1) );
 range.push_back( make_pair( .1, 2) );
 range.push_back( make_pair( .4, 2) );


 //Create The Sim Fitter and configure the parameters
 SimFitter blastWaveFit(spectraToFit,fitFuncs);
 blastWaveFit.Init(totalFitPars, &range);
 blastWaveFit.SetPrintLevel(1);
 blastWaveFit.SetParameterRelationships(parRelations);

 //Fix the Masses
 blastWaveFit.FixParameter(BW_PIONMASS);
 blastWaveFit.FixParameter(BW_KAONMASS);
 blastWaveFit.FixParameter(BW_LAMBDAMASS);
 blastWaveFit.FixParameter(BW_PROTONMASS);

 blastWaveFit.SetParLimits(BW_TRANSVELOCITYPOWER,0,1);
 blastWaveFit.FixParameter(BW_RADIALINTEGRALLIMIT);

 //Set Parameter Limits
 blastWaveFit.SetParLimits(BW_TEMP,.08,.3);
 blastWaveFit.SetParLimits(BW_SURFACEVELOCITY,0,1);      

 for (unsigned int i=0; i<nParticles; i++){
   blastWaveFit.SetParLimits(BW_AMP+i,0,100000000);
 }

 ROOT::Fit::FitResult bwFitResults = blastWaveFit.Fit();
 std::vector<TF1 *> bwFitFuncResults = blastWaveFit.GetFitFunctions();




 for ( unsigned int iPart =0 ; iPart<nParticles; iPart++){

   canvas->cd(iPart+1);
   if (PION  == particles.at(iPart) ){ TH1F * h = gPad->DrawFrame(0,0,.6,140);  h->SetTitle("Pion"); }
   if (PROTON  == particles.at(iPart) ){ TH1F * h = gPad->DrawFrame(0,0,1.2,25);  h->SetTitle("Proton"); }
   if (LAMBDA  == particles.at(iPart) ){ TH1F * h = gPad->DrawFrame(0,0,2,.5);  h->SetTitle("Lambda");}
   if (KAON  ==  particles.at(iPart) ){ TH1F * h = gPad->DrawFrame(0,0,2,4); h->SetTitle("Kaon"); }

   if (LAMBDA  == particles.at(iPart) )    spectraToFit.at(iPart)->SetMarkerStyle(kOpenCircle);
   if (KAON  == particles.at(iPart) )    spectraToFit.at(iPart)->SetMarkerStyle(kFullCircle);
   if (LAMBDA  == particles.at(iPart) )    spectraToFit.at(iPart)->SetMarkerSize(1.2);
   if (KAON  == particles.at(iPart) )    spectraToFit.at(iPart)->SetMarkerSize(1.2);

   spectraToFit.at(iPart)->Draw("P");

   //   else spectraToFit.at(iPart)->Draw("AP");
   bwFitFuncResults.at(iPart)->Draw("SAME");
 
 }



}
