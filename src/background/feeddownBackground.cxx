//This code uses the output of the background simulations to determine
//the fractional background due to feed down or secondaries

#include <iostream>
#include <vector>
#include <utility>

#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TF1.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TVirtualFitter.h>
#include <TDirectory.h>
#include <TLatex.h>

#include "globalDefinitions.h"
#include "utilityFunctions.h"
#include "ParticleInfo.h"
#include "UserCuts.h"
#include "StyleSettings.h"

Bool_t save=true;

//___________________________________________________________________________
void RemovePointsWithLargeErrors(TGraphAsymmErrors *spectrum, Double_t maxRelativeError=.1){

  //Loop Over the Points of the spectrum. Remove any point which is found
  //to have a relativeError larger than maxRelativeError

		for (int iPoint=spectrum->GetN()-1; iPoint>=0; iPoint--){
			Float_t  err = spectrum->GetErrorY(iPoint);
			Float_t yval = spectrum->GetY()[iPoint];
			if (yval == 0 && err > 0.1) spectrum->RemovePoint(iPoint);
			else if (yval != 0 && err / spectrum->GetY()[iPoint] > maxRelativeError)
				spectrum->RemovePoint(iPoint);
    }

}


//____MAIN__________________________________________
void feeddownBackground(TString inputFile, TString outputFile){

  //Create an instance of the ParticleInfo Class
  ParticleInfo *particleInfo = new ParticleInfo(GetStarLibraryVersion());

  //Get the Number of Centrality Bins
  const unsigned int nCentralityBins = 1;

  //Get the Collision Energy
  Double_t energy = GetCollisionEnergy();

  //Define the Particle Species to be used
  const int nSpecies(3);
  const int nCharges(2);
  std::vector<int> charge(nCharges);
  charge.at(0) = -1;
  charge.at(1) = 1;
	const int nRapBins = GetMaxRapidityIndexOfInterest()+1;
  
  //Open the inputFile and Load the Histograms
  TFile *inFile = new TFile(inputFile,"READ");
  inFile->cd("FeedDown");
  std::vector< std::vector < std::vector<TH3F *> > > parentIDHisto3D
    (nSpecies, std::vector< std::vector<TH3F *> >
     (nCharges, std::vector<TH3F *>
      (nCentralityBins, (TH3F *)NULL)));
  
  for (int iSpecies=0; iSpecies<nSpecies; iSpecies++){
    for (int iCharge=0; iCharge<nCharges; iCharge++){
      for (unsigned int iCentBin=0; iCentBin<nCentralityBins; iCentBin++){
	
	parentIDHisto3D.at(iSpecies).at(iCharge).at(iCentBin) =
	  (TH3F *)gDirectory->Get(Form("parentGeantIDHisto_%s_cent%02d",
				       particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data(),
				       iCentBin));
	
      }//End Loop Over Centrality Bins
    }//End Loop Over Charges
  }//End Loop Over Species

  //Open the output file
  TFile *outFile = new TFile(outputFile,save ? "UPDATE":"READ");
  outFile->cd();
  
  //Define the Graphs that will hold the background fraction and their fits
	
	//TOTAL Background graphs and fits
	std::vector< std::vector < std::vector< std::vector< TGraphAsymmErrors *> > > > backgroundGraph
     (nSpecies, std::vector< std::vector< std::vector< TGraphAsymmErrors *> > > 
        (nCharges, std::vector< std::vector< TGraphAsymmErrors *> >
           (nCentralityBins, std::vector< TGraphAsymmErrors*>
              (nRapBins, (TGraphAsymmErrors *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TF1 *> > > > backgroundGraphFit
     (nSpecies, std::vector< std::vector< std::vector< TF1 *> > > 
        (nCharges, std::vector< std::vector< TF1 *> >
           (nCentralityBins, std::vector< TF1*>
              (nRapBins, (TF1 *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TF1 *> > > > backgroundGraphFit1
     (nSpecies, std::vector< std::vector< std::vector< TF1 *> > > 
        (nCharges, std::vector< std::vector< TF1 *> >
           (nCentralityBins, std::vector< TF1*>
              (nRapBins, (TF1 *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TGraphErrors *> > > > backgroundGraphFitSys
     (nSpecies, std::vector< std::vector< std::vector< TGraphErrors *> > > 
        (nCharges, std::vector< std::vector< TGraphErrors *> >
           (nCentralityBins, std::vector< TGraphErrors*>
              (nRapBins, (TGraphErrors *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TGraphErrors *> > > > backgroundGraphFitSys1
     (nSpecies, std::vector< std::vector< std::vector< TGraphErrors *> > > 
        (nCharges, std::vector< std::vector< TGraphErrors *> >
           (nCentralityBins, std::vector< TGraphErrors*>
              (nRapBins, (TGraphErrors *)NULL))));
  
  //HADRONIC Background graphs and fits
  std::vector< std::vector < std::vector< std::vector< TGraphAsymmErrors *> > > > hadronicGraph
     (nSpecies, std::vector< std::vector< std::vector< TGraphAsymmErrors *> > > 
        (nCharges, std::vector< std::vector< TGraphAsymmErrors *> >
           (nCentralityBins, std::vector< TGraphAsymmErrors*>
              (nRapBins, (TGraphAsymmErrors *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TF1 *> > > > hadronicGraphFit
     (nSpecies, std::vector< std::vector< std::vector< TF1 *> > > 
        (nCharges, std::vector< std::vector< TF1 *> >
           (nCentralityBins, std::vector< TF1*>
              (nRapBins, (TF1 *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TGraphErrors *> > > > hadronicGraphFitSys
     (nSpecies, std::vector< std::vector< std::vector< TGraphErrors *> > >
        (nCharges, std::vector< std::vector< TGraphErrors *> >
           (nCentralityBins, std::vector< TGraphErrors*>
              (nRapBins, (TGraphErrors *)NULL))));
 
  //FEEDDOWN Background
  std::vector< std::vector < std::vector< std::vector< TGraphAsymmErrors *> > > > feeddownGraph
     (nSpecies, std::vector< std::vector< std::vector< TGraphAsymmErrors *> > >
        (nCharges, std::vector< std::vector< TGraphAsymmErrors *> >
           (nCentralityBins, std::vector< TGraphAsymmErrors*>
              (nRapBins, (TGraphAsymmErrors *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TF1 *> > > > feeddownGraphFit
     (nSpecies, std::vector< std::vector< std::vector< TF1 *> > >
        (nCharges, std::vector< std::vector< TF1 *> >
           (nCentralityBins, std::vector< TF1*>
              (nRapBins, (TF1 *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TF1 *> > > > feeddownGraphFit1
     (nSpecies, std::vector< std::vector< std::vector< TF1 *> > >
        (nCharges, std::vector< std::vector< TF1 *> >
           (nCentralityBins, std::vector< TF1*>
              (nRapBins, (TF1 *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TGraphErrors *> > > > feeddownGraphFitSys
     (nSpecies, std::vector< std::vector< std::vector< TGraphErrors *> > >
        (nCharges, std::vector< std::vector< TGraphErrors *> >
           (nCentralityBins, std::vector< TGraphErrors*>
              (nRapBins, (TGraphErrors *)NULL))));
  std::vector< std::vector < std::vector< std::vector< TGraphErrors *> > > > feeddownGraphFitSys1
     (nSpecies, std::vector< std::vector< std::vector< TGraphErrors *> > >
        (nCharges, std::vector< std::vector< TGraphErrors *> >
           (nCentralityBins, std::vector< TGraphErrors*>
              (nRapBins, (TGraphErrors *)NULL))));

  for (int iSpecies=0; iSpecies<nSpecies; iSpecies++){
    for (int iCharge=0; iCharge<nCharges; iCharge++){

      TString yAxisTitleMod = " ";

      for (unsigned int iCentBin=0; iCentBin<nCentralityBins; iCentBin++){
	for (int yIndex=GetMinRapidityIndexOfInterest(); yIndex<GetMaxRapidityIndexOfInterest()+1; yIndex++){

	  //Initializing Graphs
		//Initializing TOTAL Graphs
	  backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = new TGraphAsymmErrors();
    backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
       SetName(Form("AllBackground_%s_cent%02d_yIndex%02d",
          particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data(),
             iCentBin,yIndex));
    backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
       SetTitle(Form("%s All%sBackground | Cent=[%02d,%02d]%% | y_{%s}=[%.02f,%.02f];(m_{T}-m_{%s})^{Reco} (GeV/c^{2});All Background %sFraction",
           particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
              yAxisTitleMod.Data(), 0, 5,
                 particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
                    GetRapidityRangeLow(yIndex),GetRapidityRangeHigh(yIndex),
                       particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
                          yAxisTitleMod.Data()));

    backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
       SetMarkerStyle(particleInfo->GetParticleMarker(iSpecies,charge.at(iCharge)));
    backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
       SetMarkerColor(kBlack);
	
    //Initializing HADRONIC Graphs
    hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = new TGraphAsymmErrors();
    hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
       SetName(Form("hadronicBackground_%s_cent%02d_yIndex%02d",
          particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data(),iCentBin,yIndex));
    hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
       SetTitle(Form("%s hadronic%sBackground | Cent=[%02d,%02d]%% | y_{%s}=[%.02f,%.02f];(m_{T}-m_{%s})^{Reco} (GeV/c^{2});hadronic %sFraction",
					particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(), yAxisTitleMod.Data(), 0, 5,
						 particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
						    GetRapidityRangeLow(yIndex),GetRapidityRangeHigh(yIndex),
						       particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
						          yAxisTitleMod.Data()));

		hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
		   SetMarkerStyle(particleInfo->GetParticleMarker(iSpecies,charge.at(iCharge)));
	  hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
		   SetMarkerColor(kBlack);

		//Initializing FEEDDOWN Graphs
	  feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = new TGraphAsymmErrors();
	  feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
	     SetName(Form("feeddownBackground_%s_cent%02d_yIndex%02d",
			    particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data(), iCentBin, yIndex));
	  feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
	     SetTitle(Form("%s Feed Down%sBackground | Cent=[%02d,%02d]%% | y_{%s}=[%.02f,%.02f];(m_{T}-m_{%s})^{Reco} (GeV/c^{2});Feed Down%sFraction",
			    particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(), yAxisTitleMod.Data(),
			       iCentBin!=nCentralityBins-1 ? (int)GetCentralityPercents().at(iCentBin+1):0,
			          (int)GetCentralityPercents().at(iCentBin),
			             particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
			                GetRapidityRangeLow(yIndex),GetRapidityRangeHigh(yIndex),
			                   particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(), yAxisTitleMod.Data()));

	  feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
	     SetMarkerStyle(particleInfo->GetParticleMarker(iSpecies,charge.at(iCharge)));
	  feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
	     SetMarkerColor(kBlack);

	  //Defining the Fits
		//TOTAL Background
		backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
				new TF1(Form("%s_FitExp",backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName()),
					"[0]*TMath::Exp([1]*x)+[2]",0.04,2);

    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
       new TF1(Form("%s_FitPower",backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName()),
          "[0]*TMath::Power(x,[1])",0.04,2);
    backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameters(0.08,-2.0);
    backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(0,0,0.7);
    backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(1,-10,0);
		backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameter(2,0.02);
		backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(2,0,.05);

    backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineWidth(3);
    backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineColor(kRed);
    backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineStyle(1);

    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameters(1.0,-1.0);
    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(0,0,5.0);
    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(1,-2,0);

    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineWidth(3);
    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineColor(kBlue);
    backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineStyle(9);

    //HADRONIC Background
    hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
       new TF1(Form("%s_Fit",hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName()),
          "[0]+[1]*x",0.04,2);
		
    hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameters(0.1,-0.1);
    hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParNames("b","m");

    hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineWidth(3);
    hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineColor(kViolet+1);
    hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineStyle(10);

		//FEEDDOWN Background
		feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
			new TF1(Form("%s_FitExp",feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName()),
				"[0]*TMath::Exp([1]*x)+[2]",0.04,2);
	  
	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
	    new TF1(Form("%s_FitPower",feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName()),
		    "[0]*TMath::Power(x,[1])",0.04,2);
	  

	  feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameters(0.08,-2.0);
	  feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(0,0,0.7);
	  feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(1,-10,0);
		feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameter(2,0.02);
		feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(2,0,.05);
	  
	  feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineWidth(3);
	  feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineColor(kRed);
	  feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineStyle(1);

	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParameters(1.0,-1.0);
	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(0,0,5.0);
	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetParLimits(1,-2,0);
	  
	  
	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineWidth(3);
	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineColor(kBlue);
	  feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetLineStyle(9);

	}//End Loop Over yIndex
      }//End Loop Over Centrality Bins
    }//End Loop Over Charges
  }//End Loop Over Species
  
  //Loop Over the species and Compute the Feed Down Background Fraction
  TH1D *parentIDHisto = NULL;
  TH3F *localParentIDHisto3D = NULL;
  for (int iSpecies=0; iSpecies<nSpecies; iSpecies++){
    for (int iCharge=0; iCharge<nCharges; iCharge++){

      outFile->cd();
      outFile->mkdir(Form("%s",particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data()));
      outFile->cd(Form("%s",particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data()));

      gDirectory->mkdir("FeedDownContributionHistos");
      gDirectory->cd("FeedDownContributionHistos");
      
      for (unsigned int iCentBin=0; iCentBin<nCentralityBins; iCentBin++){

				localParentIDHisto3D = parentIDHisto3D.at(iSpecies).at(iCharge).at(iCentBin);
				
				//Loop Over the rapidity bins
				for (int xBin=1; xBin<=localParentIDHisto3D->GetNbinsX(); xBin++){

					//Convert the histogram bin indices to the rapidity indidces
					int yIndex = GetRapidityIndex(localParentIDHisto3D->GetXaxis()->GetBinCenter(xBin));
					if (yIndex > GetMaxRapidityIndexOfInterest() || yIndex < GetMinRapidityIndexOfInterest()) continue;

					TH1D totalHisto("totalHisto","totalHisto",localParentIDHisto3D->GetNbinsY(),mTm0Min,mTm0Max);
					//total pass:
					TH1D passHisto("passHisto","passHisto",localParentIDHisto3D->GetNbinsY(),mTm0Min,mTm0Max);
					TH1D feedPassHisto("feedPassHisto","feedPassHisto",localParentIDHisto3D->GetNbinsY(),mTm0Min,mTm0Max);
					TH1D hadronicPassHisto("hadronicPassHisto","hadronicPassHisto",localParentIDHisto3D->GetNbinsY(),mTm0Min,mTm0Max);
					
					//Loop Over the mTm0 Bins
					for (int yBin=1; yBin<localParentIDHisto3D->GetNbinsY(); yBin++){

						int mTm0Index = GetmTm0Index(localParentIDHisto3D->GetYaxis()->GetBinCenter(yBin));
						parentIDHisto = localParentIDHisto3D->ProjectionZ("pz",xBin,xBin,yBin,yBin);

						parentIDHisto->SetName(Form("feeddownContribution_%s_cent%02d_yIndex%02d_mTm0Index%02d",
							 particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data(), iCentBin, yIndex, mTm0Index));
						parentIDHisto->SetTitle(Form("Parent GeantID Distribution of %s | #sqrt{s_{NN}}=%.01f GeV | Cent=[%02d,%02d]%% | y_{%s}=[%.02f,%.02f] | m_{T}-m_{%s}=[%.03f,%.03f] (GeV/c^{2})",
							 particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(), energy,
							 iCentBin!=nCentralityBins-1 ? (int)GetCentralityPercents().at(iCentBin+1):0,
							 (int)GetCentralityPercents().at(iCentBin),
							 particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
							 localParentIDHisto3D->GetXaxis()->GetBinLowEdge(xBin),
							 localParentIDHisto3D->GetXaxis()->GetBinLowEdge(xBin) + localParentIDHisto3D->GetXaxis()->GetBinWidth(1),
							 particleInfo->GetParticleSymbol(iSpecies,charge.at(iCharge)).Data(),
							 localParentIDHisto3D->GetYaxis()->GetBinLowEdge(yBin),
							 localParentIDHisto3D->GetYaxis()->GetBinLowEdge(yBin) + localParentIDHisto3D->GetYaxis()->GetBinWidth(1)));

						//Make sure there are enough entries
						if (parentIDHisto->GetEntries() <= 5){
							delete parentIDHisto;
							continue;
						}

						//Set Bin Labels
						for (int i=1; i<parentIDHisto->GetNbinsX(); i++){
							parentIDHisto->GetXaxis()->SetBinLabel(i,GetParticleSymbolUsingGeantID(i-1));
						}
						parentIDHisto->GetXaxis()->LabelsOption("d");
						parentIDHisto->GetXaxis()->SetRangeUser(0,32);
						parentIDHisto->SetFillColor(GetCentralityColor(iCentBin));
						parentIDHisto->GetYaxis()->SetTitle("Counts");
						
						Double_t total = parentIDHisto->Integral(parentIDHisto->FindBin(0),parentIDHisto->FindBin(32));
						//total pass:
						Double_t pass = parentIDHisto->Integral(parentIDHisto->FindBin(1),parentIDHisto->FindBin(32));
						Double_t passHadronic = parentIDHisto->Integral(parentIDHisto->FindBin(7),parentIDHisto->FindBin(9)) + parentIDHisto->Integral(parentIDHisto->FindBin(13),parentIDHisto->FindBin(15));
						Double_t passFeed = parentIDHisto->Integral(parentIDHisto->FindBin(1),parentIDHisto->FindBin(6)) + parentIDHisto->Integral(parentIDHisto->FindBin(10),parentIDHisto->FindBin(12)) + parentIDHisto->Integral(parentIDHisto->FindBin(16),parentIDHisto->FindBin(32) );

						totalHisto.SetBinContent(yBin,total);
						passHisto.SetBinContent(yBin,pass);
						hadronicPassHisto.SetBinContent(yBin,passHadronic);
						feedPassHisto.SetBinContent(yBin,passFeed);

						parentIDHisto->Write(parentIDHisto->GetName(),TObject::kOverwrite);	    
						
						delete parentIDHisto;
						parentIDHisto = NULL;
						
					}//End Loop Over mTm0 axis

					passHisto.Rebin(2);
					hadronicPassHisto.Rebin(2);
					feedPassHisto.Rebin(2);
					totalHisto.Rebin(2);

					backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->BayesDivide(&passHisto,&totalHisto);
					hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->BayesDivide(&hadronicPassHisto,&totalHisto);
					feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->BayesDivide(&feedPassHisto,&totalHisto);

					//Skip Graphs with too few points
          if (backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetN() < 2){ 
		            delete backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
		            delete backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
		            delete backgroundGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
		
		            backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
		            backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
		            backgroundGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
		
		            continue;
		      }

		      if (hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetN() < 2){ 
				            delete hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
				            delete hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
				            delete hadronicGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
				
				            hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
				            hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
				            hadronicGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
				
				            continue;
				  }


					if (feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetN() < 2){
						delete feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
						delete feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
						delete feeddownGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex);
						
						feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
						feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
						feeddownGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) = NULL;
						
						continue;
					}

					//Clean up the Graphs
					TGraphChop((TGraphErrors *)backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),.05,true);
					TGraphChop((TGraphErrors *)backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),1.5,false);
					TGraphChop((TGraphErrors *)hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),.05,true);
          TGraphChop((TGraphErrors *)hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),1.5,false);
					TGraphChop((TGraphErrors *)feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),.05,true);
					TGraphChop((TGraphErrors *)feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),1.5,false);

					//Get the min X value of the first point in the graph
					Double_t minX = backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetX()[0] -
					backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetEXlow()[0];

          Double_t hadronicMinX = hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetX()[0] -
          hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetEXlow()[0];

          Double_t feeddownMinX = feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetX()[0] -
          feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetEXlow()[0];

					//Fit the Graph and Get the Confidence Interval
          //Fit for all
					backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetRange(minX,2);
					backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
					   Fit(backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),"EX0R");
					backgroundGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
					   GetConfidenceIntervalOfFit(backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex));

					backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetRange(minX,2);
					backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
					   Fit(backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),"EX0R+");
					backgroundGraphFitSys1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
					   GetConfidenceIntervalOfFit(backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex));

          //fit for hadronic
          hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetRange(hadronicMinX,2);
          hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
             Fit(hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),"EX0R");
          hadronicGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
             GetConfidenceIntervalOfFit(hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex));

          //Remove largest errors to make plot scale more readable
          RemovePointsWithLargeErrors(hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex), 10);

          //fit for feeddown
					feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetRange(feeddownMinX,2);
					feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
					   Fit(feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),"EX0R");
					feeddownGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
						 GetConfidenceIntervalOfFit(feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex));

					feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->SetRange(feeddownMinX,2);
					feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Fit(feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex),"EX0R+");
					feeddownGraphFitSys1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex) =
						 GetConfidenceIntervalOfFit(feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex));
					
				cout<<"yIndex finished: "<<yIndex<<endl;
				}//End Loop Over rapidity axis
			}//End Loop Over centrality bins
    }//End Loop Over charges
  }//End Loop Over species
  //Save Everything if Requested
  if (!save)
    return;

  for (int iSpecies=0; iSpecies<nSpecies; iSpecies++){
    for (int iCharge=0; iCharge<nCharges; iCharge++){

      outFile->cd();
      if (!outFile->cd(Form("%s",particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data()))){
				outFile->mkdir(Form("%s",particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data()));
				outFile->cd(Form("%s",particleInfo->GetParticleName(iSpecies,charge.at(iCharge)).Data()));
      }

			gDirectory->mkdir("BackgroundGraphs");
      gDirectory->mkdir("HadronicBackgroundGraphs");
      gDirectory->mkdir("FeedDownBackgroundGraphs");

      gDirectory->mkdir("BackgroundFits");
      gDirectory->mkdir("HadronicBackgroundFits");
      gDirectory->mkdir("FeedDownBackgroundFits");
      
      for (unsigned int iCentBin=0; iCentBin<nCentralityBins; iCentBin++){
				for (int yIndex=GetMinRapidityIndexOfInterest(); yIndex<GetMaxRapidityIndexOfInterest()+1; yIndex++){
					
          //Saving Total Background Graphs and Fits
          gDirectory->cd("BackgroundGraphs");
					if (!backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)){
						gDirectory->cd("..");
						continue;
					}
					backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
					   Write(backgroundGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(), TObject::kOverwrite);

          gDirectory->cd("..");
					gDirectory->cd("BackgroundFits");
					backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
  					 Write(backgroundGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);
					backgroundGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(backgroundGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);
					backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(backgroundGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);
					backgroundGraphFitSys1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(feeddownGraphFitSys1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);

          gDirectory->cd("..");
					gDirectory->cd("HadronicBackgroundGraphs");
					if (!hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)){
						gDirectory->cd("..");
						continue;
					}
					hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(hadronicGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(), TObject::kOverwrite);

					gDirectory->cd("..");
					gDirectory->cd("HadronicBackgroundFits");
					hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(hadronicGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
								TObject::kOverwrite);
					hadronicGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
					   Write(hadronicGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
								TObject::kOverwrite);

					gDirectory->cd("..");
					gDirectory->cd("FeedDownBackgroundGraphs");
					if (!feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)){
						gDirectory->cd("..");
						continue;
					}
					feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(feeddownGraph.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(), TObject::kOverwrite);
					
					gDirectory->cd("..");
					gDirectory->cd("FeedDownBackgroundFits");
					feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(feeddownGraphFit.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);
					feeddownGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(feeddownGraphFitSys.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);
					feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(feeddownGraphFit1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);
					feeddownGraphFitSys1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->
						 Write(feeddownGraphFitSys1.at(iSpecies).at(iCharge).at(iCentBin).at(yIndex)->GetName(),
						    TObject::kOverwrite);

					gDirectory->cd("..");
					
				}//End Loop Over yIndex
      }//End Loop Over Cent Bin
    }//End Loop Over Charge
  }//End Loop Over Species  
      
}
