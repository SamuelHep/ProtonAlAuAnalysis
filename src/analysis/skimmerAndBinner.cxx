//Applies cuts to skim the data at the event, vertex, and track level. Then bins the
//tracks in centrality, rapidity, and transverse mass for various particle species.
//Can be run with tofCut flag but requires a zTOFGaus.root file 


#include <iostream>
#include <vector>
#include <utility>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>

#include "globalDefinitions.h"
#include "TrackInfo.h"
#include "PrimaryVertexInfo.h"
#include "EventInfo.h"
#include "DavisDstReader.h"
#include "ParticleInfo.h"
#include "StRefMultExtendedCorr.h"
#include "UserCuts.h"
#include "utilityFunctions.h"
#include "fitZTPCUtilities.h"

//____Checks if zTOF variable is in 2 Sigma of fit dEdx_______________________
Bool_t tofMatch(TF1 * gaus,Double_t zTOF);


//___MAIN_____________________________________________
void skimmerAndBinner(TString inputDataFile,TString starLibrary, 
		      Long64_t nEvents=-1, TString outputFile="",TString inputTofFitFile="",
		      Double_t energy=0,Int_t tofCut = 1, TString eventConfig=""){

  //Read the DavisDst
  DavisDstReader *davisDst = NULL;
  if (GetZVertexCuts().first != -999 && GetZVertexCuts().second != 999)
    davisDst = new DavisDstReader(inputDataFile,GetZVertexCuts().first,GetZVertexCuts().second);
  else
    davisDst = new DavisDstReader(inputDataFile);
  //Create Pointers needed for reading the tree
  TrackInfo *track = NULL;
  PrimaryVertexInfo *primaryVertex = NULL;
  EventInfo *event = NULL;
  //If no output file was specified then we won't produce one
  TFile *outFile = NULL;
  if (outputFile.CompareTo(""))
    outFile = new TFile(outputFile,"RECREATE");
    TFile * tofFile = NULL;
  
    if (tofCut == 1){

    tofFile = new TFile(inputTofFitFile,"READ");
    
    if (!tofFile) {
      cout << "No Tof input file found... exiting" << endl;
    return;
  }
  }


  //Create an instance of the Particle Info Class
  ParticleInfo *particleInfo = new ParticleInfo(starLibrary,true,TRUNCATED);

  //************************************************************************
  //Temporary solution since nSigma there is not an nSigma 
  //value for deteruon,triton,helio,alpha,muon
  const int nParticles = 3;
  //************************************************************************

  //Get the centrality vector from the User's function
  const int nCentralityBins = GetNCentralityBins(); 
  
  //Declare Histogram to hold the binned Particles and the Number of Events

  //Contains the number of events in each centrality bin
  TH1D *nEventsHisto = new TH1D("nEvents","Events in Each Centrality Bin",
				nCentralityBins,0,nCentralityBins);
  std::vector<TH1D *> centVarHisto(nCentralityBins,(TH1D *)NULL);

  //The Yield Histograms for the difference centrality bins and particle species
  std::vector<std::vector<TH3D *> > ZTPCHistoPlus(nParticles,vector<TH3D *>(nCentralityBins,(TH3D *)NULL));
  std::vector<std::vector<TH3D *> > ZTOFHistoPlus(nParticles,vector<TH3D *>(nCentralityBins,(TH3D *)NULL));
  std::vector<std::vector<TH3D *> > ZTPCHistoMinus(nParticles,vector<TH3D *>(nCentralityBins,(TH3D *)NULL));
  std::vector<std::vector<TH3D *> > ZTOFHistoMinus(nParticles,vector<TH3D *>(nCentralityBins,(TH3D *)NULL));

  //Number of tracks seen in the TPC and TOF - Will be used for TOF Matching Efficiency
  std::vector<std::vector<TH2D *> > TPCTracksPlus(nParticles,vector<TH2D *>(nCentralityBins,(TH2D *)NULL));
  std::vector<std::vector<TH2D *> > TOFTracksPlus(nParticles,vector<TH2D *>(nCentralityBins,(TH2D *)NULL));
  std::vector<std::vector<TH2D *> > TPCTracksMinus(nParticles,vector<TH2D *>(nCentralityBins,(TH2D *)NULL));
  std::vector<std::vector<TH2D *> > TOFTracksMinus(nParticles,vector<TH2D *>(nCentralityBins,(TH2D *)NULL));

  //Yield Histograms with no centrality binning 
  std::vector<TH3D *>  ZTPCHistoProtonPlusTOFCut(5,(TH3D *)NULL);
  std::vector<TH3D *> ZTOFHistoPlusAllCent(nParticles,(TH3D *)NULL);
  std::vector<TH3D *> ZTPCHistoProtonMinusTOFCut(5,(TH3D *)NULL);
  std::vector<TH3D *> ZTOFHistoMinusAllCent(nParticles,(TH3D *)NULL);

  //proton ToF fits
  std::vector<std::vector<std::vector<TF1 *> > > tofProtonParticleFit
    (5,std::vector<std::vector<TF1*> >
     (nRapidityBins,std::vector<TF1 *>
      (nmTm0Bins,(TF1 *)NULL)));

  
  //Loop Over the Particle Species and centrality bins and define the histograms
  for (Int_t iParticle=0; iParticle<nParticles; iParticle++){
    
    const char *particleSymbolPlus = particleInfo->GetParticleSymbol(iParticle,1).Data();
    const char *particleSymbolMinus = particleInfo->GetParticleSymbol(iParticle,-1).Data();

    //Skip Electrons and Muons
    if (iParticle == ELECTRON || iParticle == MUON)
      continue;
    
    for (Int_t iCentBin=0; iCentBin<nCentralityBins; iCentBin++){
      
      if (!centVarHisto.at(iCentBin)){
	centVarHisto.at(iCentBin) = 
	  new TH1D(Form("centralityVariable_Cent%d",iCentBin),
		   Form("CentralityVariable Distribution in Centrality Bin %d",iCentBin),
		   1000,0,1000);
      }

      //Define Yield Histograms
      ZTPCHistoPlus.at(iParticle).at(iCentBin) = 
	new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-5,10);
      ZTOFHistoPlus.at(iParticle).at(iCentBin) = 
	new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-1,3);

      ZTPCHistoMinus.at(iParticle).at(iCentBin) =
        new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-5,10);
      ZTOFHistoMinus.at(iParticle).at(iCentBin) =
        new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-1,3);

      TPCTracksPlus.at(iParticle).at(iCentBin) =
	new TH2D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max);
      TOFTracksPlus.at(iParticle).at(iCentBin) = 
	new TH2D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max);
      
      TPCTracksMinus.at(iParticle).at(iCentBin) =
	new TH2D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max);
      TOFTracksMinus.at(iParticle).at(iCentBin) =
	new TH2D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max);
 

      //Set the Names and Titles of the Yield Histograms

      ZTPCHistoPlus.at(iParticle).at(iCentBin)->
	SetName(Form("zTPC_%s_Cent%d",particleInfo->GetParticleName(iParticle,1).Data(),iCentBin));
      ZTPCHistoPlus.at(iParticle).at(iCentBin)->
        SetTitle(Form("Z_{TPC} of %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TPC,%s}",
		      particleSymbolPlus,iCentBin,particleSymbolPlus,
		      particleSymbolPlus,particleSymbolPlus));


      ZTOFHistoPlus.at(iParticle).at(iCentBin)->
        SetName(Form("zTOF_%s_Cent%d",particleInfo->GetParticleName(iParticle,1).Data(),iCentBin));
      ZTOFHistoPlus.at(iParticle).at(iCentBin)->
	SetTitle(Form("Z_{TOF} of %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TOF,%s}",
		      particleSymbolPlus,iCentBin,particleSymbolPlus,
		      particleSymbolPlus,particleSymbolPlus));


      ZTPCHistoMinus.at(iParticle).at(iCentBin)->
        SetName(Form("zTPC_%s_Cent%d",particleInfo->GetParticleName(iParticle,-1).Data(),iCentBin));
      ZTPCHistoMinus.at(iParticle).at(iCentBin)->
	SetTitle(Form("Z_{TPC} of %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TPC,%s}",
		      particleSymbolMinus,iCentBin,particleSymbolMinus,
		      particleSymbolMinus,particleSymbolMinus));

      ZTOFHistoMinus.at(iParticle).at(iCentBin)->
	SetName(Form("zTOF_%s_Cent%d",particleInfo->GetParticleName(iParticle,-1).Data(),iCentBin));
      ZTOFHistoMinus.at(iParticle).at(iCentBin)->
	SetTitle(Form("Z_{TOF} of %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TOF,%s}",
                      particleSymbolMinus,iCentBin,particleSymbolMinus,
                      particleSymbolMinus,particleSymbolMinus));


      TPCTracksPlus.at(iParticle).at(iCentBin)->
	SetName(Form("nTPCTracks_%s_Cent%d",particleInfo->GetParticleName(iParticle,1).Data(),iCentBin));
      TPCTracksPlus.at(iParticle).at(iCentBin)->
	SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2})",
		      particleSymbolPlus,iCentBin,particleSymbolPlus,particleSymbolPlus));

      TOFTracksPlus.at(iParticle).at(iCentBin)->
	SetName(Form("nTOFTracks_%s_Cent%d",particleInfo->GetParticleName(iParticle,1).Data(),iCentBin));
      TOFTracksPlus.at(iParticle).at(iCentBin)->
	SetTitle(Form("n TOF Tracks for %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2})",
                      particleSymbolPlus,iCentBin,particleSymbolPlus,particleSymbolPlus));

      TPCTracksMinus.at(iParticle).at(iCentBin)->
	SetName(Form("nTPCTracks_%s_Cent%d",particleInfo->GetParticleName(iParticle,-1).Data(),iCentBin));
      TPCTracksMinus.at(iParticle).at(iCentBin)->
	SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2})",
                      particleSymbolMinus,iCentBin,particleSymbolMinus,particleSymbolMinus));

      TOFTracksMinus.at(iParticle).at(iCentBin)->
        SetName(Form("nTOFTracks_%s_Cent%d",particleInfo->GetParticleName(iParticle,-1).Data(),iCentBin));
      TOFTracksMinus.at(iParticle).at(iCentBin)->
        SetTitle(Form("n TOF Tracks for %s Centrality Bin: %d;y_{%s};m_{T}-m_{%s} (GeV/c^{2})",
                      particleSymbolMinus,iCentBin,particleSymbolMinus,particleSymbolMinus));
    

    }//End Loop Over Centrality Bins

  }


  //Fill the All Centrality ZTOF Histograms 
  for (Int_t iParticle=0; iParticle<nParticles; iParticle++){


    const char *particleSymbolPlus = particleInfo->GetParticleSymbol(iParticle,1).Data();
    const char *particleSymbolMinus = particleInfo->GetParticleSymbol(iParticle,-1).Data();

    TString cutLabel = "_tofCut";

      ZTOFHistoPlusAllCent.at(iParticle) = 
	new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-1,3);

      ZTOFHistoMinusAllCent.at(iParticle) =
        new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-1,3);

      ZTOFHistoPlusAllCent.at(iParticle)->
        SetName(Form("zTOF_%s_AllCent%s",particleInfo->GetParticleName(iParticle,1).Data(),cutLabel.Data()));
      ZTOFHistoPlusAllCent.at(iParticle)->
	SetTitle(Form("Z_{TOF} of %s All Centrality;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TOF,%s}",
		      particleSymbolPlus,particleSymbolPlus,
		      particleSymbolPlus,particleSymbolPlus));

      ZTOFHistoMinusAllCent.at(iParticle)->
	SetName(Form("zTOF_%s_AllCent%s",particleInfo->GetParticleName(iParticle,-1).Data(),cutLabel.Data()));
      ZTOFHistoMinusAllCent.at(iParticle)->
	SetTitle(Form("Z_{TOF} of %s All Centrality;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TOF,%s}",
                      particleSymbolMinus,particleSymbolMinus,
                      particleSymbolMinus,particleSymbolMinus));

  }

  //Fill the Proton Contamination Histograms. One for each confounding particle
  for (Int_t iParticle=0; iParticle<5; iParticle++){


    const char *particleSymbolPlus = particleInfo->GetParticleSymbol(iParticle,1).Data();
    const char *particleSymbolMinus = particleInfo->GetParticleSymbol(iParticle,-1).Data();


    //Fill Histograms for ZTPC fitting width constraint
    ZTPCHistoProtonPlusTOFCut.at(iParticle) = 
	new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-5,10);

      ZTPCHistoProtonMinusTOFCut.at(iParticle) =
        new TH3D("","",nRapidityBins,rapidityMin,rapidityMax,nmTm0Bins,mTm0Min,mTm0Max,1000,-5,10);

      TString cutLabel = "_tofCut";


      ZTPCHistoProtonPlusTOFCut.at(iParticle)->
	SetName(Form("zTPC_ProtonPlus_con_%s_AllCent%s",particleInfo->GetParticleName(iParticle,1).Data(),cutLabel.Data()));
      ZTPCHistoProtonPlusTOFCut.at(iParticle)->
        SetTitle(Form("Z_{TPC} with ToF Cut of ProtonPlus con %s All Centrality;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TPC,%s}",
		      particleSymbolPlus,particleSymbolPlus,
		      particleSymbolPlus,particleSymbolPlus));

      ZTPCHistoProtonMinusTOFCut.at(iParticle)->
	SetName(Form("zTPC_ProtonMinus_con_%s_AllCent%s",particleInfo->GetParticleName(iParticle,-1).Data(),cutLabel.Data()));     
      ZTPCHistoProtonMinusTOFCut.at(iParticle)->
	SetTitle(Form("Z_{TPC} with ToF Cut of %s All Centrality;y_{%s};m_{T}-m_{%s} (GeV/c^{2});Z_{TPC,%s}",
		      particleSymbolMinus,particleSymbolMinus,
		      particleSymbolMinus,particleSymbolMinus));

      //Load tofCut Gaussian Functions 
  if (tofCut == 1){
      for (int rapidityBin=0;rapidityBin<nRapidityBins;rapidityBin++){
	for (int mTm0Bin=0;mTm0Bin<nmTm0Bins;mTm0Bin++){
	  
	  if( !( tofFile->Get(Form("%sTofFit_yIndex%02d_mTm0Index%02d",particleInfo->GetParticleName(iParticle,1).Data(),rapidityBin,mTm0Bin)) ) ) continue;
	  tofProtonParticleFit.at(iParticle).at(rapidityBin).at(mTm0Bin) = (TF1*) tofFile->Get(Form("%sTofFit_yIndex%02d_mTm0Index%02d",particleInfo->GetParticleName(iParticle,1).Data(),rapidityBin,mTm0Bin));
	}
      }
  }
  }//End Loop Over particle species


  //If the user has passed a number of events to process then use it,
  //otherwise use all the entries in the tree
  Long64_t nEntries(0);
  if (nEvents > 0)
    nEntries = nEvents;
  else
    nEntries = davisDst->GetEntries();

  //-------------------------------------------------------------------------------
  //SKIM DATA: Loop Over the tree to apply the event and vertex cuts. Track cuts
  //           are applied during the binning phase below.
  //-------------------------------------------------------------------------------
  std::vector<std::pair<int,std::vector<int> > > goodEntries;
  davisDst->SetBranchStatus("TrackInfo.*",0);
  for (Int_t iEntry=0; iEntry<nEntries; iEntry++){

    //Get the ith entry and check if it passes the cuts
    event = davisDst->GetEntry(iEntry);
    if (!IsGoodEvent(event))
      continue;

    std::vector<int> tempVertexIndex(0);

    //Loop over the primary vertex array of this event
    Int_t nPrimaryVertices = event->GetNPrimaryVertices();
    for (Int_t iPrimaryVertex=0; iPrimaryVertex<nPrimaryVertices; iPrimaryVertex++){

      //Get the ith primary vertex and check if it passes the cuts 
      primaryVertex = event->GetPrimaryVertex(iPrimaryVertex);
      if (!IsGoodVertex(primaryVertex))
        continue;

      //Fill the the tempVertexIndexVector with the index of this vertex
      tempVertexIndex.push_back(iPrimaryVertex);

    }//End Loop Over primary Vertices

    //If there was at least one good primary vertex then add this event to the 
    //vector of goodEntries
    if (tempVertexIndex.size() > 0)
      goodEntries.push_back(std::make_pair(iEntry,tempVertexIndex));

  }//End Loop Over Events
  //END APPLY CUTS

  //--------------------------------------------------------------------------------
  //BIN DATA: Loop Over the entries in the goodEntries vector. Bin the events in 
  //          centrality and the tracks in rapidity and transverse mass if they
  //          pass the track cuts.  
  //--------------------------------------------------------------------------------
  davisDst->SetBranchStatus("*",1);
  for (std::vector<std::pair<int,std::vector<int> > >::iterator goodEntry = goodEntries.begin(); 
       goodEntry != goodEntries.end(); ++goodEntry){
    
    //Get the good Entry
    event = davisDst->GetEntry((*goodEntry).first);
    
    //Loop Over the Primary Vertex Array of this event
    for (std::vector<int>::iterator iPrimaryVertex = ((*goodEntry).second).begin(); 
	 iPrimaryVertex != ((*goodEntry).second).end(); ++iPrimaryVertex){
      
      //Get the ith primary vertex
      primaryVertex = event->GetPrimaryVertex(*iPrimaryVertex);

      //Figure out the centrality bin here
      Int_t centralityVariable = GetCentralityVariable(primaryVertex);
      Int_t centralityBin = GetCentralityBin(centralityVariable);

      //Skip this vertex if it has a bad centrality bin determination
      if (centralityBin < 0)
	continue;

      //Fill the Centrality Histograms
      nEventsHisto->Fill(centralityBin);
      centVarHisto.at(centralityBin)->Fill((double)centralityVariable);

      //Loop Over the Primary Tracks of this primary vertex
      Int_t nPrimaryTracks = primaryVertex->GetNPrimaryTracks();
      for (Int_t iPrimaryTrack=0; iPrimaryTrack<nPrimaryTracks; iPrimaryTrack++){

	//Get the ith primary track and check if it passes the cuts
	track = primaryVertex->GetPrimaryTrack(iPrimaryTrack);
	if (!IsGoodTrack(track))
	  continue;

	//Loop over Each of the Particle Species and Fill the Yield Histograms
	for (Int_t iParticle=0; iParticle<nParticles; iParticle++){

	  //Skip Electrons and Muons
	  if (iParticle == ELECTRON || iParticle == MUON)
	    continue;

	  Double_t rapidity = track->GetRapidity(particleInfo->GetParticleMass(iParticle));
	  Double_t mTm0 = track->GetmTm0(particleInfo->GetParticleMass(iParticle));
	  Double_t nSigmaTPC = track->GetNSigmaOfSpecies(iParticle);
	  
	  //Compute the zTPC variable and fill the yield histogram
	  Double_t zTPC = particleInfo->ComputeZTPC(track->GetdEdx(TRUNCATED),track->GetPt(),
						    track->GetPz(),iParticle);

	  //Compute the zTOF variable and fill the yield Histogram
	  Double_t zTOF = particleInfo->ComputeZTOF(1.0/track->GetTofBeta(),track->GetPt(),
						    track->GetPz(),iParticle);
	  
	  if (track->GetCharge() > 0){
	    ZTPCHistoPlus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0,zTPC);
	    if (fabs(nSigmaTPC) < 2.0)
	      TPCTracksPlus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0);
	  }
	  else if (track->GetCharge() < 0){
	    ZTPCHistoMinus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0,zTPC);
	    if (fabs(nSigmaTPC) < 2.0)
	      TPCTracksMinus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0);
	  }
	  //Fill Histogram for Proton confounding particle Histograms
	  if (tofCut == 1){
	    if (track->GetCharge() > 0){
	      
	      if (iParticle == PROTON ){
		for (Int_t conParticle = 0; conParticle < 5; conParticle++){

		  Int_t rapidityIndex = GetRapidityIndex(rapidity);
		  Int_t mTm0Index = GetmTm0Index(mTm0);
		  Bool_t matchedTrack = false;
		  
		  if( rapidityIndex == -1 || mTm0Index == -1) continue;
		  if(tofProtonParticleFit.at(conParticle).at(rapidityIndex).at(mTm0Index))
		    {
		      matchedTrack = tofMatch( tofProtonParticleFit.at(conParticle).at(rapidityIndex).at(mTm0Index),zTOF ); 
		      if (!matchedTrack) continue;
		      ZTPCHistoProtonPlusTOFCut.at(conParticle)->Fill(rapidity,mTm0,zTPC);
		      
		    }
		}	
	      }
	    }
	  }
	
	  //Check to make sure this track has a TOF hit and a physical 1/Beta
	  if (track->GetTofMatchFlag() < 1 || 1.0/track->GetTofBeta() < 1.0)
	    continue;
	  
	  if (track->GetCharge() > 0){
	    ZTOFHistoPlus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0,zTOF);
	    ZTOFHistoPlusAllCent.at(iParticle)->Fill(rapidity,mTm0,zTOF);
	    if (fabs(nSigmaTPC) < 2.0)
	      TOFTracksPlus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0);
	  }
	  else if (track->GetCharge() < 0){
	    ZTOFHistoMinus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0,zTOF);
	    ZTOFHistoMinusAllCent.at(iParticle)->Fill(rapidity,mTm0,zTOF);
	    if (fabs(nSigmaTPC) < 2.0)
	      TOFTracksMinus.at(iParticle).at(centralityBin)->Fill(rapidity,mTm0);
	  }


	}//End Loop Over Particle Species

      }//End Loop Over Primary Tracks

    }//End Loop Over Primary Vertices

  }//End Loop Over Events/Entry


  //Write the Histograms to File
  outFile->cd();
  nEventsHisto->Write();

  outFile->mkdir("CentralityVariableHistograms");
  outFile->cd("CentralityVariableHistograms");
  for (Int_t iCentBin=0; iCentBin<nCentralityBins; iCentBin++){
    centVarHisto.at(iCentBin)->Write();
  }
  outFile->cd();

  for (Int_t iParticle=0; iParticle<nParticles; iParticle++){

    //Skip Electrons and Muons
    if (iParticle == ELECTRON || iParticle == MUON)
      continue;

    outFile->mkdir(Form("%s",particleInfo->GetParticleName(iParticle,-1).Data()));
    outFile->mkdir(Form("%s",particleInfo->GetParticleName(iParticle,1).Data()));
    
    outFile->cd(Form("%s",particleInfo->GetParticleName(iParticle,-1).Data()));
    ZTOFHistoMinusAllCent.at(iParticle)->Write();
    outFile->cd();
    
    outFile->cd(Form("%s",particleInfo->GetParticleName(iParticle,1).Data()));
    ZTOFHistoPlusAllCent.at(iParticle)->Write();
    outFile->cd();
      
    for (Int_t iCentBin=0; iCentBin<nCentralityBins; iCentBin++){

      outFile->cd(Form("%s",particleInfo->GetParticleName(iParticle,1).Data()));
      ZTPCHistoPlus.at(iParticle).at(iCentBin)->Write();
      ZTOFHistoPlus.at(iParticle).at(iCentBin)->Write();
      TPCTracksPlus.at(iParticle).at(iCentBin)->Write();
      TOFTracksPlus.at(iParticle).at(iCentBin)->Write();

      outFile->cd();

      outFile->cd(Form("%s",particleInfo->GetParticleName(iParticle,-1).Data()));
      ZTPCHistoMinus.at(iParticle).at(iCentBin)->Write();
      ZTOFHistoMinus.at(iParticle).at(iCentBin)->Write();
      TPCTracksMinus.at(iParticle).at(iCentBin)->Write();
      TOFTracksMinus.at(iParticle).at(iCentBin)->Write();
      
      outFile->cd();
    }

  }//End Loop Over particle species
  

  //Filling Confounding particle histograms
  if (tofCut == 1){
    
    for (Int_t iParticle=0; iParticle<5; iParticle++){
      
      //Skip Electrons and Muons
      if (iParticle == ELECTRON || iParticle == MUON)
	continue;
      
      outFile->cd(Form("%s",particleInfo->GetParticleName(2,1).Data()));
      ZTPCHistoProtonPlusTOFCut.at(iParticle)->Write();
      outFile->cd();
    }
  }
  outFile->Close();

}

//Simple function checks if zTOF variable is in Gaussian and returns bool
Bool_t tofMatch(TF1 * gaus,Double_t zTOF)
{

  if (!gaus){ return false;}
 
  Double_t center = gaus->GetParameter(1);
  Double_t range = 2*gaus->GetParameter(2);

  if ( abs(zTOF - center) < range ) { return true;}

  return false;

}
