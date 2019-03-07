#include <iostream>
#include <vector>

using namespace std;

void RunCentralityDetermination(TString DATAFILE, TString DATAHISTO,
				TString GLAUBERFILE,TString OUTFILE,
				Double_t STARTBIN=20, Double_t STOPBIN=-1){

  //Load the Necessary Libraries
  gSystem->Load("../bin/GlauberClass_cxx");
  gSystem->Load("../bin/GlauberModel_cxx");
  gSystem->Load("../bin/FindBestFitNegativeBinomialParameters_cxx");
  gSystem->Load("../bin/FindCentralityBinCuts_cxx");
  gSystem->Load("../bin/DetermineCentralityBin_cxx");
  gSystem->Load("../bin/ReturnStructs_cxx");
  gSystem->Load("../bin/FindNpartNcollDistributions_cxx");
  gSystem->Load("../bin/SystematicErrors_cxx");
  gSystem->Load("../bin/PrintResults_cxx");
  gSystem->Load("../bin/MakeResultsTable_cxx");

  //Load the Data Histogram that we are trying to match
  TFile *dataFile = new TFile(DATAFILE);
  TH1D *dataHisto = (TH1D *)dataFile->Get(DATAHISTO);

  //Load the File from the Glauber Simulation
  GlauberClass *glauberEvent = 0;
  TFile *glauberFile = new TFile(GLAUBERFILE);
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  glauberTree->GetEntry(0);
  
  //Create an output file to save things to
  TFile *outFile = new TFile(OUTFILE,"RECREATE");
  
  //Variables Used in the Loops Below
  Double_t startBinCenter(STARTBIN);    //The value of the bin center where we want to start matching
  Double_t stopBinCenter(STOPBIN);      //The value of the bin center where we want to stop matching
  Int_t nTrials(500);                   //The number of choices of (npp,k) permitted
  Int_t nEvents = glauberTree->GetEntries();
  Bool_t useTrackingEfficiency = false;
  NegBinomialSearchResults bestBinomialParameters;

  //Check to make sure we have everything we need
  if (!dataHisto){
    cout <<"ERROR - RunCentralityDetermination.C -- Data histogram not found! ABORTING!" <<endl;
    exit (EXIT_FAILURE);
  }
  if (!glauberTree){
    cout <<"ERROR - RunCentralityDetermination.C -- Glauber Tree not found! ABORTING!" <<endl;
    exit (EXIT_FAILURE);
  }
  
  //________________________________________________________________________________________________
  //First, we need to find the parameters of the negative binomial distribution
  //which results in a multiplicity distrubution that best matches the data.

  cout <<"*** Step 1 of 5: Searching for the Best Negative Binomial Parameters.... " <<endl;
  bestBinomialParameters =
    FindBestFitNegativeBinomialParameters(dataHisto,glauberFile,
					  startBinCenter,stopBinCenter,
					  useTrackingEfficiency,nTrials);

  cout <<"  Complete!" <<endl;
  cout <<"  Inverse Chi2: " <<bestBinomialParameters.InverseChi2
       <<" npp = " <<bestBinomialParameters.npp
       <<" k = "   <<bestBinomialParameters.k <<endl;
  
  //Save the Default information to the file
  outFile->cd();
  (bestBinomialParameters.bestFitHisto)->SetName("bestFit");
  (bestBinomialParameters.chi2Graph)->SetName("chiSquaredGraph");
  (bestBinomialParameters.bestFitHisto)->Write();
  (bestBinomialParameters.chi2Graph)->Write();
  dataHisto->Write();

  //_______________________________________________________________________________________________
  //Second, Find the centrality bin cuts

  cout <<"*** Step 2 of 5: Determining the centrality bin cuts... " <<endl;
  FindCentralityBinCuts(&bestBinomialParameters);
  cout <<"  Complete!" <<endl;

  for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.centralityBinCuts.size(); iCentBin++){
    cout <<Form("Percent Centrality: %.02F \t Cut: %.05F",
		bestBinomialParameters.centralityBinDefinitions.at(iCentBin),
		bestBinomialParameters.centralityBinCuts.at(iCentBin))
	 <<endl;
  }

  //_______________________________________________________________________________________________
  //Third, Make the Npart, Ncoll, and Impact Parameter Histograms for all the centrality Bins

  cout <<"*** Step 3 of 5: Constructing the Npart, Ncoll, and Impact Parameter Distributions... " <<endl;
  FindNpartNcollDistributions(&bestBinomialParameters,glauberTree);
  cout <<"  Complete!" <<endl;
  
  //_______________________________________________________________________________________________
  //Fourth, Do the systematic Error Investigation

  cout <<"*** Step 4 of 5: Determining the Systematic Errors...." <<endl;
  SystematicErrors(nEvents,glauberEvent->GetNNucleonsNucA(),glauberEvent->GetNNucleonsNucB(),
		   glauberEvent->GetInelasticCrossSection(),&bestBinomialParameters,
		   dataHisto,useTrackingEfficiency);
  cout <<"  Complete!" <<endl;

  //_______________________________________________________________________________________________
  //Fifth, Print Results and Save

  cout <<"*** Step 5 of 5: Printing Results and Saving" <<endl;
  for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.nCentralityBins; iCentBin++){
    PrintResults(&bestBinomialParameters,iCentBin);
  }

  //Make the Results Table
  MakeResultsTable(&bestBinomialParameters);
  
  //Save
  outFile->cd();
  const int nCentBins = bestBinomialParameters.nCentralityBins;
  outFile->WriteObject(&TArrayD(1,&bestBinomialParameters.npp),"npp");
  outFile->WriteObject(&TArrayD(1,&bestBinomialParameters.k),"k");
  outFile->WriteObject(&TArrayD(1,&bestBinomialParameters.InverseChi2),"inverseChi2");
  outFile->WriteObject(&TArrayI(1,&bestBinomialParameters.nCentralityBins),"nCentralityBins");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.centralityBinDefinitions.at(0)),"centralityBinDefinitions");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.centralityBinCuts.at(0)),"centralityBinCuts");

  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nPartMeans.at(0)),"nPartMeans");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nPartStatErrors.at(0)),"nPartStatErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nPartSysErrors.at(0)),"nPartSysErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nCollMeans.at(0)),"nCollMeans");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nCollStatErrors.at(0)),"nCollStatErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nCollSysErrors.at(0)),"nCollSysErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.impactParamMeans.at(0)),"impactParamMeans");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.impactParamStatErrors.at(0)),"impactParamStatErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.impactParamSysErrors.at(0)),"impactParamSysErrors");

  outFile->mkdir("GlauberHistograms");
  outFile->cd("GlauberHistograms");

  bestBinomialParameters.nPartTotalHisto->SetName("nPartTotal");
  bestBinomialParameters.nCollTotalHisto->SetName("nCollTotal");
  bestBinomialParameters.impactParamTotalHisto->SetName("impactParamTotal");
  
  bestBinomialParameters.nPartTotalHisto->Write();
  bestBinomialParameters.nCollTotalHisto->Write();
  bestBinomialParameters.impactParamTotalHisto->Write();
  
  for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.nCentralityBins; iCentBin++){

    bestBinomialParameters.nPartHistos[iCentBin].SetName(Form("nPart_cent%d",iCentBin));
    bestBinomialParameters.nCollHistos[iCentBin].SetName(Form("nColl_cent%d",iCentBin));
    bestBinomialParameters.impactParamHistos[iCentBin].SetName(Form("impactParam_cent%d",iCentBin));
    
    bestBinomialParameters.nPartHistos[iCentBin].Write();
    bestBinomialParameters.nCollHistos[iCentBin].Write();
    bestBinomialParameters.impactParamHistos[iCentBin].Write();
    
  }

  outFile->cd();
  outFile->mkdir("Canvases");
  outFile->cd("Canvases");

  bestBinomialParameters.nPartTotalSysErrCanvas->SetName("nPartTotalSysErrCanvas");
  bestBinomialParameters.nCollTotalSysErrCanvas->SetName("nCollTotalSysErrCanvas");
  bestBinomialParameters.impactParamTotalSysErrCanvas->SetName("impactParamTotalSysErrCanvas");
  bestBinomialParameters.resultsTableCanvas->SetName("resultsTableCanvas");
  
  bestBinomialParameters.nPartTotalSysErrCanvas->Write();
  bestBinomialParameters.nCollTotalSysErrCanvas->Write();
  bestBinomialParameters.impactParamTotalSysErrCanvas->Write();
  bestBinomialParameters.resultsTableCanvas->Write();
  
  //Close the File
  outFile->cd();
  outFile->Close();


  cout <<"FINISHED SUCCESSFULLY." <<endl;
  return;
  
}
