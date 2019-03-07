#include <iostream>
#include <utility>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TThread.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMutex.h>

#include "GlauberClass.h"
#include "ReturnStructs.h"

//Globals
TGraph2D *chi2Graph;
Int_t normStartBin;
Int_t normStopBin;
Int_t nGlauberEntries;
Int_t *nPartArr;
Int_t *nCollArr;

TMutex dataHistoMutex;
TMutex chi2GraphMutex;

Bool_t useTrackingEfficiencyGl;

TH1D *MakeMultiplicityForNB(TH1D *dataHisto,TFile *glauberFile,
			    Double_t startMatchingBinCenter,
			    Double_t stopMatchingBinCenter,
			    Double_t npp, Double_t k);

struct argStruct{
  TH1D *dataHisto;
  Double_t npp;
  Double_t k;
  Double_t chi2;
  Int_t index;
};

//________________________________________________________________
std::pair<double,double> GetNppRange(){

  const Int_t n = chi2Graph->GetN();
  const Int_t top10percent = TMath::CeilNint(n*.10);

  Int_t index[n];
  TMath::Sort(n,chi2Graph->GetZ(),index);

  Double_t nppOfTops[top10percent];
  for (Int_t i=0; i<top10percent; i++){
    nppOfTops[i] = chi2Graph->GetX()[index[i]];
  }

  Double_t nppMin = TMath::MinElement(top10percent,nppOfTops);
  Double_t nppMax = TMath::MaxElement(top10percent,nppOfTops);

  return std::make_pair(nppMin,nppMax);


}

//________________________________________________________________
std::pair<double,double> GetKRange(){

  const Int_t n = chi2Graph->GetN();
  const Int_t top10percent = TMath::CeilNint(n*.10);

  Int_t index[n];
  TMath::Sort(n,chi2Graph->GetZ(),index);

  Double_t kOfTops[top10percent];
  for (Int_t i=0; i<top10percent; i++){
    kOfTops[i] = chi2Graph->GetY()[index[i]];
  }

  Double_t kMin = TMath::MinElement(top10percent,kOfTops);
  Double_t kMax = TMath::MaxElement(top10percent,kOfTops);;

  return std::make_pair(kMin,kMax);
}

//_________________________________________________________________
std::pair<double,double> GetBestPoint(){

  const Int_t n = chi2Graph->GetN();

  Int_t index[n];
  TMath::Sort(n,chi2Graph->GetZ(),index);

  return std::make_pair(chi2Graph->GetX()[index[0]],chi2Graph->GetY()[index[0]]);
}

//_________________________________________________________________
Double_t ComputeHistogramChiSquared(TH1D *hData, TH1D *hSim, Int_t startBin, Int_t stopBin){

  Double_t chi2(0);

  for (Int_t iBin=startBin; iBin<stopBin; iBin++){

    //Don't consider points from data histo that have no error
    if (hData->GetBinError(iBin) == 0.0)
      continue;

    Double_t num = hData->GetBinContent(iBin) - hSim->GetBinContent(iBin);
    Double_t denom = hData->GetBinError(iBin);
    
    chi2 += pow(num/denom,2);
    
  }

  return chi2;
}

//___________________________________________________________________
void *DoParticleProduction(void *args){

  argStruct *argVals = (argStruct *)args;
  Double_t npp = argVals->npp;
  Double_t k   = argVals->k;
  Int_t index  = argVals->index;

  //Create an instance of the GlauberClass
  GlauberClass glauberEvent;
  glauberEvent.SetNegativeBinomialParameters(npp,k);

  dataHistoMutex.Lock();
  TH1D dataHisto(*argVals->dataHisto);
  dataHistoMutex.UnLock();

  Double_t nParticles(0);
  TH1D *simHisto = new TH1D(Form("hTemp_%d",index),Form("hTemp_%d",index),
			 dataHisto.GetNbinsX(),
			 dataHisto.GetBinLowEdge(1),
			 dataHisto.GetBinLowEdge(dataHisto.GetNbinsX())+dataHisto.GetBinWidth(1));

  for (Int_t i=0; i<nGlauberEntries; i++){
    nParticles = glauberEvent.ProduceParticles(nPartArr[i],nCollArr[i],-1,useTrackingEfficiencyGl);
    simHisto->Fill(nParticles);
  }


  //Prepare the simulated histo for comparison to data histo
  simHisto->Sumw2();
  Double_t scaleFactor = simHisto->Integral(normStartBin,simHisto->GetNbinsX());

  //If this simHisto has a bad scale factor then just return
  if (scaleFactor < 1){
    TThread::Lock();
    cout <<"Thread " <<index <<" had a bad scale factor" <<endl;
    TThread::UnLock();
    return (void *)NULL;
  }

  //Scale the Histogram
  simHisto->Scale(1.0/(Double_t)simHisto->Integral(normStartBin,simHisto->GetNbinsX()));


  Double_t chi2 = ComputeHistogramChiSquared(&dataHisto,simHisto,normStartBin,normStopBin);

  chi2GraphMutex.Lock();
  chi2Graph->SetPoint(chi2Graph->GetN(),npp,k,1.0/chi2);
  chi2GraphMutex.UnLock();

  delete simHisto;

  return NULL;
}



//______________________________________________________________________________________________
NegBinomialSearchResults
FindBestFitNegativeBinomialParameters(TH1D *dataHisto,TFile *glauberFile,
				      Double_t startMatchingBinCenter,
				      Double_t stopMatchingBinCenter,
				      Bool_t useTrackingEfficiency,
				      const Int_t nThreads=1000){

  useTrackingEfficiencyGl = useTrackingEfficiency;
  
  Double_t nppMin(0.1); //0.1
  Double_t nppMax(1.0); //2.0
  Double_t kMin(0.0);   //0.0
  Double_t kMax(5.0);     //5

  //The bin for which will will start the matchin procedure
  normStartBin = dataHisto->FindBin(startMatchingBinCenter);

  //This is the bin for which we will end the matching procedure.
  //By default the last bin of the dataHisto is used as the stopping bin.
  //However, if the user has passed in a valid value for the stop bin then use it instead.
  normStopBin = dataHisto->GetNbinsX();
  if (stopMatchingBinCenter > startMatchingBinCenter)
    normStopBin = dataHisto->FindBin(stopMatchingBinCenter);

  //Build the Npart and Ncoll Arrays from the file so that we can do the
  //particle production procedure using multiple threads
  GlauberClass *glauberEvent = 0;
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);

  const Int_t nEntries = glauberTree->GetEntries();
  nGlauberEntries = nEntries;
  nPartArr = new Int_t[nEntries];
  nCollArr = new Int_t[nEntries];

  for (Int_t i=0; i<glauberTree->GetEntries(); i++){
    glauberTree->GetEntry(i);
    nPartArr[i] = glauberEvent->GetNParticipants();
    nCollArr[i] = glauberEvent->GetNBinaryCollisions();
  }//End Loop Over Glauber Tree


  //Make sure we have everything we need
  if (!dataHisto){
    cout <<"ERROR: FindBestFitNegativeBinomialParameters() - dataHisto does not exist!" <<endl;
    exit (EXIT_FAILURE);
  }

  //Prepare the Data Histogram...do this only if sumw2 hasn't been called before
  //this is to protect the histogram in the case of mulitple iterations
  if ( ((TArrayD *)dataHisto->GetSumw2())->GetSize() == 0){
    dataHisto->Sumw2();
  }
  dataHisto->Scale(1.0/(Double_t)dataHisto->Integral(normStartBin,normStopBin));
  dataHisto->SetMarkerColor(kBlack);
  dataHisto->SetMarkerStyle(kFullCircle);
  
  TRandom3 rand(0);
  chi2Graph = new TGraph2D();
  chi2Graph->SetMarkerStyle(kFullCircle);

  const Int_t nLoops(3);
  const Int_t nConcurrentThreads(25);
  for (Int_t iLoop=0; iLoop<nLoops; iLoop++){

    Int_t nSubmittedThreads(0);

    //Decrease the number of threads for each loop
    Int_t nThreadsToSubmit = nThreads/(double)(iLoop+1);

    //While the number of running threads is less than the
    //number of number of threads to be submitted for this loop...
    while (nSubmittedThreads < nThreadsToSubmit){

      //Only submit more threads when there are less than
      //nConcurrentThreads running
      std::vector <argStruct *> argVector;
      std::vector <TThread *>   threadVector;
      while (TThread::Exists() < nConcurrentThreads){

	argVector.push_back(new argStruct);
	argVector.back()->index = nSubmittedThreads;
	argVector.back()->dataHisto = dataHisto;

	//Set the nb parameters based on which loop we are in
	if (iLoop == 0){
	  argVector.back()->npp = rand.Uniform(nppMin,nppMax);
	  argVector.back()->k   = rand.Uniform(kMin,kMax);
	}
	else if (iLoop == 1){
	  std::pair<double,double> nppRange = GetNppRange();
	  argVector.back()->npp = rand.Uniform(nppRange.first,nppRange.second);
	  argVector.back()->k   = rand.Uniform(kMin,kMax);
	}
	else{
	  Double_t tempNpp(0), tempK(0);
	  chi2Graph->GetHistogram()->GetRandom2(tempNpp, tempK);
	  argVector.back()->npp = tempNpp;
	  argVector.back()->k   = tempK;
	}

	//Create the new thread and then run it
	threadVector.push_back(new TThread(Form("thread_%d",nSubmittedThreads),DoParticleProduction,
					   (void *)argVector.back()));
	threadVector.back()->Run();

	//Incremenet the number of threads that have been submitted
	nSubmittedThreads++;


      }//End Loop Over submitting nConcurrent Threads

      //After submitting the nConcurrent Threads wait for them to be done
      //before starting the next batch
      for (unsigned int i=0; i<threadVector.size(); i++){
	threadVector.at(i)->Join();
	threadVector.at(i)->Delete();
	delete argVector.at(i);
      }

    }//End Loop Over submitted Threads

  }//End Loops

  //Best Npp and K
  std::pair<double,double> bestPoint = GetBestPoint();
  
  NegBinomialSearchResults results;
  results.npp = bestPoint.first;
  results.k = bestPoint.second;
  results.InverseChi2 = chi2Graph->Interpolate(results.npp,results.k);
  results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
					       startMatchingBinCenter,stopMatchingBinCenter,
					       results.npp,results.k);
  results.chi2Graph = (TGraph2D *)chi2Graph->Clone(chi2Graph->GetName());

  //Clean Up
  delete nPartArr;
  delete nCollArr;
  delete chi2Graph;
  
  return results;

}


//_______________________________________________________________________________________
TH1D *MakeMultiplicityForNB(TH1D *dataHisto,TFile *glauberFile,
			   Double_t startMatchingBinCenter,
			   Double_t stopMatchingBinCenter,
			   Double_t npp, Double_t k ){

  //The bin for which will will start the matchin procedure
  normStartBin = dataHisto->FindBin(startMatchingBinCenter);

  //This is the bin for which we will end the matching procedure.
  //By default the last bin of the dataHisto is used as the stopping bin.
  //However, if the user has passed in a valid value for the stop bin then use it instead.
  normStopBin = dataHisto->GetNbinsX();
  if (stopMatchingBinCenter > startMatchingBinCenter)
    normStopBin = dataHisto->FindBin(stopMatchingBinCenter);

  //Prepare the Data Histogram...do this only if sumw2 hasn't been called before
  //this is to protect the histogram in the case of mulitple iterations
  if ( ((TArrayD *)dataHisto->GetSumw2())->GetSize() == 0){
    dataHisto->Sumw2();
    dataHisto->Scale(1.0/(Double_t)dataHisto->Integral(normStartBin,dataHisto->GetNbinsX()));
    dataHisto->SetMarkerColor(kBlack);
    dataHisto->SetMarkerStyle(kFullCircle);
  }
  
  //Build the Npart and Ncoll Arrays from the file so that we can do the
  //particle production procedure using multiple threads
  GlauberClass *glauberEvent = 0;
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  glauberEvent->SetNegativeBinomialParameters(npp,k);

  Double_t nParticles(0);
  TH1D *simHisto = new TH1D("htemp","hTemp",
			    dataHisto->GetNbinsX(),
			    dataHisto->GetBinLowEdge(1),
			    dataHisto->GetBinLowEdge(dataHisto->GetNbinsX())+dataHisto->GetBinWidth(1));

  const Int_t nEntries = glauberTree->GetEntries();
  for (Int_t i=0; i<nEntries; i++){
    glauberTree->GetEntry(i);
    nParticles = glauberEvent->ProduceParticles(useTrackingEfficiencyGl);
    simHisto->Fill(nParticles);

  }//End Loop Over Glauber Tree

  //Prepare the simulated histo for comparison to data histo
  simHisto->Sumw2();
  Double_t scaleFactor = simHisto->Integral(normStartBin,simHisto->GetNbinsX());
  simHisto->Scale(1.0/scaleFactor);
  simHisto->SetLineColor(kRed);

  return simHisto;
}
