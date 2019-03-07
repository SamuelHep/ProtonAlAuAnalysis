#include <iostream>
#include <utility>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TThread.h>

#include "GlauberClass.h"
#include "ReturnStructs.h"
#include "DetermineCentralityBin.h"

//Structure for arguments to passed into threads
struct threadArgs{

  Int_t nEvents, index;         //Number of Events to do glauber for, index of thread
  Int_t nNucleonsA, nNucleonsB; //Number of nucleons in each nucleus
  Double_t nnCrossSection;      //Inelastic Cross Section
  Int_t model;                  //Nucleon Distribution Model
  Double_t npp, k;              //Negative Binomial Parameters
  Double_t rHard;               //Nucleon Hard-Core Radius
  std::vector<double> *cuts;    //Array of centrality cuts
  Double_t *woodsSaxonPars;     //Array of Woods-Saxon Parameters (use NULL for default values)
  Double_t x;                   //Hard Soft scaling valule (use -1 for default value);
  Bool_t useTrackingEfficiency; //Whether or not to use tracking efficiency

};

TCanvas *canvas;
TCanvas *canvas1;
TCanvas *canvas2;

TH1D *nPartHistoArray[20][16];
TH1D *nCollHistoArray[20][16];
TH1D *impactParamHistoArray[20][16];

std::vector<double> nPartMeans, nCollMeans, impactParamMeans; 
std::vector<double> nPartSysErrors, nCollSysErrors, impactParamSysErrors;

//___________________________________________________________________________________________________________
void *LoopOverGlauberEvents(void *args){

  threadArgs *argVals = (threadArgs *)args;
  
  GlauberClass glauberEvent(argVals->nNucleonsA,argVals->nNucleonsB,argVals->nnCrossSection);

  if (argVals->woodsSaxonPars)
    glauberEvent.SetWoodsSaxonParameters(argVals->woodsSaxonPars);
  
  glauberEvent.SetNegativeBinomialParameters(argVals->npp,argVals->k);

  const int nCentBins = argVals->cuts->size();
  for (Int_t centBin=0; centBin<nCentBins; centBin++){
    nPartHistoArray[argVals->index][centBin] =
      new TH1D(Form("nPartHisto_%d_%d",argVals->index,centBin),"nPartHisto",400,0,400);
    nCollHistoArray[argVals->index][centBin] =
      new TH1D(Form("nCollHisto_%d_%d",argVals->index,centBin),"nCollHisto",1000,0,1000);
    impactParamHistoArray[argVals->index][centBin] =
      new TH1D(Form("impactParamHisto_%d_%d",argVals->index,centBin),"nCollHisto",100,0,25);
  }
  
  Int_t iEvent(0);
  Int_t nParticles(0);

  Double_t maxImpactParameter =
    1.3 * (glauberEvent.GetNuclearRadius(argVals->nNucleonsA)+glauberEvent.GetNuclearRadius(argVals->nNucleonsB));

  while (iEvent < argVals->nEvents){
    
    glauberEvent.SetImpactParameter(0,maxImpactParameter);

    glauberEvent.GenerateNucleonPositions(argVals->model);

    glauberEvent.CountNpartNcoll();

    if (glauberEvent.GetNParticipants() < 1)
      continue;

    if (argVals->x < 0)
      nParticles = glauberEvent.ProduceParticles(argVals->useTrackingEfficiency);
    else
      nParticles = glauberEvent.ProduceParticles(argVals->x,argVals->useTrackingEfficiency);
    
    Int_t centBin = DetermineCentralityBin(*argVals->cuts, nParticles,nCentBins);

    if (centBin < 0)
      continue;
    
    nPartHistoArray[argVals->index][centBin]->Fill(glauberEvent.GetNParticipants());
    nCollHistoArray[argVals->index][centBin]->Fill(glauberEvent.GetNBinaryCollisions());
    impactParamHistoArray[argVals->index][centBin]->Fill(glauberEvent.GetImpactParameter());
    
    iEvent++;

    glauberEvent.Reset();
  }

  return NULL;
}

//_________________________________________________________________________________________________________
void SystematicErrors(const Int_t nEvents, const Int_t nNucleonsA, const Int_t nNucleonsB,
		      const Double_t defaultInelSigma, const Double_t npp, const Double_t k,
		      std::vector<double> *cuts, Double_t percentCutVary = 0.05,
		      Bool_t useTrackingEfficiency=false){

  //All the histograms in this code will not be associated with a directory
  TH1::AddDirectory(false);
  
  const int nCentBins(cuts->size());
  std::vector<double> cutsLow, cutsHigh;

  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
    cutsLow.push_back(cuts->at(iCentBin) * (1.0 - percentCutVary));
    cutsHigh.push_back(cuts->at(iCentBin) * (1.0 + percentCutVary));
  }

  //Create all threads with default values. The values will be changed below
  const Int_t nThreads(16);
  TThread *threads[nThreads];
  threadArgs *argArray[nThreads];
  for (Int_t iThread=0; iThread<nThreads; iThread++){
    argArray[iThread]                 = new threadArgs;
    argArray[iThread]->nEvents        = nEvents;
    argArray[iThread]->index          = iThread;
    argArray[iThread]->nnCrossSection = defaultInelSigma;
    argArray[iThread]->model           = kWoodsSaxon;
    argArray[iThread]->nNucleonsA     = nNucleonsA;
    argArray[iThread]->nNucleonsB     = nNucleonsB;
    argArray[iThread]->npp            = npp;
    argArray[iThread]->k              = k;
    argArray[iThread]->rHard          = 0;
    argArray[iThread]->cuts           = cuts;
    argArray[iThread]->woodsSaxonPars = NULL;
    argArray[iThread]->x              = -1;
    argArray[iThread]->useTrackingEfficiency = useTrackingEfficiency;
    
    threads[iThread] = new TThread(Form("thread_%d",iThread),LoopOverGlauberEvents,(void *)argArray[iThread]);

  }


  for (Int_t iIndex=0; iIndex<nThreads; iIndex++){
    for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
      nPartHistoArray[iIndex][iCentBin] = NULL;
      nCollHistoArray[iIndex][iCentBin] = NULL;
      impactParamHistoArray[iIndex][iCentBin] = NULL;
    }
  }
  
  //0. Run with all the parameters at their central values
  threads[0]->Run();

  //1. Systematic Errors Due to changing npp to upper val
  argArray[1]->npp = argArray[1]->npp * 1.05;
  threads[1]->Run();
  
  //2. Systematic Error Due to changing npp to lower val
  argArray[2]->npp = argArray[2]->npp * 0.95;
  threads[2]->Run();
 
  //3. Change k to upper val
  argArray[3]->k = argArray[3]->k * 1.05;
  threads[3]->Run();
 
  //4. Change k to lower val
  argArray[4]->k = argArray[4]->k * 0.95;
  threads[4]->Run();
  
  //5. Change sigma to upper val
  argArray[5]->nnCrossSection = argArray[5]->nnCrossSection * 1.05;
  threads[5]->Run();
  
  //6. Change sigma to lower val
  argArray[6]->nnCrossSection = argArray[5]->nnCrossSection * 0.95;
  threads[6]->Run();

  //7. Require nucleon separation of .8 = .4+.4
  argArray[7]->rHard = 0.4;
  threads[7]->Run();

  //Create a GlauberClass object to get the WS parameters
  GlauberClass g(nNucleonsA,nNucleonsB,defaultInelSigma);
  std::vector<double> wsPars8 = g.GetWoodsSaxonParameters(nNucleonsA);
  std::vector<double> wsPars9(wsPars8);
  std::vector<double> wsPars10(wsPars8);
  std::vector<double> wsPars11(wsPars8);
  
  //8. Change half density radius to lower val
  wsPars8[2] = wsPars8[2]*0.98;
  argArray[8]->woodsSaxonPars = &wsPars8.at(0);
  threads[8]->Run();

  //9. Change half density radius to upper error
  wsPars9[2] = wsPars9[2]*1.02;
  argArray[9]->woodsSaxonPars = &wsPars9.at(0);
  threads[9]->Run();

  //10. Change half density radius to lower val
  wsPars10[3] = wsPars10[3]*0.90;
  argArray[10]->woodsSaxonPars = &wsPars10.at(0);
  threads[10]->Run();

  //11. Change half density radius to upper error
  wsPars11[3] = wsPars11[3]*1.10;
  argArray[11]->woodsSaxonPars = &wsPars11.at(0);
  threads[11]->Run();

  Double_t defaultXHardSoft = g.GetHardSoftScale();
  
  //12. Change the hard->soft scaling factor,x to lower val
  argArray[12]->x = defaultXHardSoft * 0.975;
  threads[12]->Run();
    
  //13. Change the hard->soft scaling factor, x to upper val
  argArray[13]->x = defaultXHardSoft * 1.025;
  threads[13]->Run();

  //14. Change the values of the multiplicity cuts by minus 5%
  argArray[14]->cuts = &cutsLow;
  threads[14]->Run();

  //15. Change the values of the multiplicity cuts by plus 5%
  argArray[15]->cuts = &cutsHigh;
  threads[15]->Run();

  //Join all the threads
  for (Int_t iThread=0; iThread < nThreads; iThread++){
    threads[iThread]->Join();
    threads[iThread]->Delete();
    delete argArray[iThread];
  }

  
  canvas = new TCanvas("canvas0","canvas0",20,20,800,600);
  canvas->SetTopMargin(.05);
  canvas->SetRightMargin(.05);
  
  TH1F *frame = canvas->DrawFrame(0,0.6,415,1.4);
  frame->SetTitle(";N_{part};N_{part}^{X}/N_{part}^{Default}");
  frame->SetTitleFont(62,"X");
  frame->SetTitleFont(62,"Y");

  canvas1 = new TCanvas("canvas1","canvas1",20,20,800,600);
  canvas1->SetTopMargin(.05);
  canvas1->SetRightMargin(.05);

  TH1F *frame1 = canvas1->DrawFrame(0,0.6,995,1.4);
  frame1->SetTitle(";N_{coll};N_{coll}^{X}/N_{coll}^{Default}");
  frame1->SetTitleFont(62,"X");
  frame1->SetTitleFont(62,"Y");

  canvas2 = new TCanvas("canvas2","canvas2",20,20,800,600);
  canvas2->SetTopMargin(.05);
  canvas2->SetRightMargin(.05);

  TH1F *frame2 = canvas2->DrawFrame(0,0.6,20,1.4);
  frame2->SetTitle(";b (fm);b^{X}/b^{Default}");
  frame2->SetTitleFont(62,"X");
  frame2->SetTitleFont(62,"Y");
  
  TGraphErrors *nPartGraph[nThreads];
  TGraphErrors *nCollGraph[nThreads];
  TGraphErrors *impactParamGraph[nThreads];
  for (Int_t iIndex=0; iIndex<nThreads; iIndex++){
    nPartGraph[iIndex] = new TGraphErrors();
    nCollGraph[iIndex] = new TGraphErrors();
    impactParamGraph[iIndex] = new TGraphErrors();
  }

  TGraphErrors *nPartTotalSysErr = new TGraphErrors();
  nPartTotalSysErr->SetFillColor(kGray+2);
  nPartTotalSysErr->SetFillStyle(3001);

  nPartGraph[0]->SetMarkerColor(kRed);
  
  nPartGraph[0]->SetMarkerStyle(kFullStar);
  nPartGraph[1]->SetMarkerStyle(kFullSquare);
  nPartGraph[2]->SetMarkerStyle(kFullDiamond);
  nPartGraph[3]->SetMarkerStyle(kFullCross);
  nPartGraph[4]->SetMarkerStyle(kFullTriangleUp);
  nPartGraph[5]->SetMarkerStyle(kFullTriangleDown);
  nPartGraph[6]->SetMarkerStyle(kFullCircle);
  nPartGraph[7]->SetMarkerStyle(kMultiply);
  nPartGraph[8]->SetMarkerStyle(kOpenCircle);
  nPartGraph[9]->SetMarkerStyle(kOpenDiamond);
  nPartGraph[10]->SetMarkerStyle(kOpenSquare);
  nPartGraph[11]->SetMarkerStyle(kOpenCross);
  nPartGraph[12]->SetMarkerStyle(kOpenTriangleUp);
  nPartGraph[13]->SetMarkerStyle(kOpenTriangleDown);
  nPartGraph[14]->SetMarkerStyle(kPlus);
  nPartGraph[15]->SetMarkerStyle(kStar);
  
  nPartGraph[0]->SetMarkerSize(1.5);
  for (Int_t i=1; i<nThreads; i++){
    nPartGraph[i]->SetMarkerSize(1.3);
  }
  
  TGraphErrors *nCollTotalSysErr = new TGraphErrors();
  nCollTotalSysErr->SetFillColor(kGray+2);
  nCollTotalSysErr->SetFillStyle(3001);


  TGraphErrors *impactParamTotalSysErr = new TGraphErrors();
  impactParamTotalSysErr->SetFillColor(kGray+2);
  impactParamTotalSysErr->SetFillStyle(3001);
  
  nCollGraph[0]->SetMarkerColor(kRed);
  impactParamGraph[0]->SetMarkerColor(kRed);
  for (Int_t iGraph=0; iGraph<nThreads; iGraph++){
    nCollGraph[iGraph]->SetMarkerStyle(nPartGraph[iGraph]->GetMarkerStyle());
    nCollGraph[iGraph]->SetMarkerSize(nPartGraph[iGraph]->GetMarkerSize());
    impactParamGraph[iGraph]->SetMarkerStyle(nPartGraph[iGraph]->GetMarkerStyle());
    impactParamGraph[iGraph]->SetMarkerSize(nPartGraph[iGraph]->GetMarkerSize());
  }

  for (Int_t iIndex=0; iIndex<nThreads; iIndex++){

    for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){

      Double_t defaultNpartMean = nPartHistoArray[0][iCentBin]->GetMean();
      Double_t nPartMean = nPartHistoArray[iIndex][iCentBin]->GetMean();      
      nPartGraph[iIndex]->SetPoint(nPartGraph[iIndex]->GetN(),defaultNpartMean,nPartMean/defaultNpartMean);

      Double_t defaultNcollMean = nCollHistoArray[0][iCentBin]->GetMean();
      Double_t nCollMean = nCollHistoArray[iIndex][iCentBin]->GetMean();
      nCollGraph[iIndex]->SetPoint(nCollGraph[iIndex]->GetN(),defaultNcollMean,nCollMean/defaultNcollMean);

      Double_t defaultImpactParamMean = impactParamHistoArray[0][iCentBin]->GetMean();
      Double_t impactParamMean = impactParamHistoArray[iIndex][iCentBin]->GetMean();
      impactParamGraph[iIndex]->SetPoint(impactParamGraph[iIndex]->GetN(),defaultImpactParamMean,
					 impactParamMean/defaultImpactParamMean);
      
    }//End Loop Over iCentBin
    
  }//End Loop Over iIndex
  
  //Compute the Total Systematic Error
  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
    
    Double_t nPartSumSquares(0);
    Double_t nCollSumSquares(0);
    Double_t impactParamSumSquares(0);
    
    for (Int_t iIndex=1; iIndex<nThreads; iIndex++){
      
      nPartSumSquares += pow(nPartGraph[iIndex]->GetY()[iCentBin] - nPartGraph[0]->GetY()[iCentBin],2);
      nCollSumSquares += pow(nCollGraph[iIndex]->GetY()[iCentBin] - nCollGraph[0]->GetY()[iCentBin],2);
      impactParamSumSquares += pow(impactParamGraph[iIndex]->GetY()[iCentBin] -
				   impactParamGraph[0]->GetY()[iCentBin],2);
    }//End Loop Over iIndex
    
    nPartTotalSysErr->SetPoint(nPartTotalSysErr->GetN(),nPartGraph[0]->GetX()[iCentBin],nPartGraph[0]->GetY()[iCentBin]);
    nPartTotalSysErr->SetPointError(nPartTotalSysErr->GetN()-1,2,sqrt(nPartSumSquares));

    nCollTotalSysErr->SetPoint(nCollTotalSysErr->GetN(),nCollGraph[0]->GetX()[iCentBin],nCollGraph[0]->GetY()[iCentBin]);
    nCollTotalSysErr->SetPointError(nCollTotalSysErr->GetN()-1,5,sqrt(nCollSumSquares));

    impactParamTotalSysErr->SetPoint(impactParamTotalSysErr->GetN(),impactParamGraph[0]->GetX()[iCentBin],
				     impactParamGraph[0]->GetY()[iCentBin]);
    impactParamTotalSysErr->SetPointError(impactParamTotalSysErr->GetN()-1,0.25,sqrt(impactParamSumSquares));
    
  }//End Loop Over iCentBins

  //Draw Things
  canvas->cd();
  nPartTotalSysErr->Draw("2");
  nPartGraph[0]->Draw("P");
  nPartGraph[1]->Draw("P");
  nPartGraph[2]->Draw("P");
  nPartGraph[3]->Draw("P");
  nPartGraph[4]->Draw("P");
  nPartGraph[5]->Draw("P");
  nPartGraph[6]->Draw("P");
  nPartGraph[7]->Draw("P");
  nPartGraph[8]->Draw("P");
  nPartGraph[9]->Draw("P");
  nPartGraph[10]->Draw("P");
  nPartGraph[11]->Draw("P");
  nPartGraph[12]->Draw("P");
  nPartGraph[13]->Draw("P");
  nPartGraph[14]->Draw("P");
  nPartGraph[15]->Draw("P");

  nPartGraph[0]->Draw("P");
  
  TPaveText *title = new TPaveText(.60,.83,.94,.93,"brNDC");
  title->SetBorderSize(0);
  title->SetFillColor(kWhite);
  title->SetTextSize(.045);
  title->SetTextFont(62);
  title->SetTextAlign(31);
  title->AddText("Contributions to N_{part} Systematic Error");
  title->Draw("SAME");

  
  TLegend *leg = new TLegend(.13,.13,.92,.41);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetNColumns(3);
  leg->SetTextSize(.035);
  leg->SetTextFont(62);
  leg->AddEntry(nPartGraph[0],"Default ","P");
  leg->AddEntry(nPartGraph[5],"#sigma_{inel}^{NN}+5%","P");
  leg->AddEntry(nPartGraph[6],"#sigma_{inel}^{NN}-5%","P");
  leg->AddEntry(nPartGraph[1],"NB n_{pp}+5%","P");
  leg->AddEntry(nPartGraph[3],"NB k+5%","P");
  leg->AddEntry(nPartGraph[2],"NB n_{pp}-5%","P");
  leg->AddEntry(nPartGraph[4],"NB k-5%","P");
  leg->AddEntry(nPartGraph[7],"NN Hard-Core r=.4","P");
  leg->AddEntry(nPartGraph[8],"Fermi-Dist. c-2%","P");
  leg->AddEntry(nPartGraph[10],"Fermi-Dist. a-10%","P");
  leg->AddEntry(nPartGraph[9],"Fermi-Dist. c+2%","P");
  leg->AddEntry(nPartGraph[11],"Fermi-Dist. a+10%","P");
  leg->AddEntry(nPartGraph[12],"Hard->Soft x-2.5%","P");
  leg->AddEntry(nPartGraph[13],"Hard->Soft x+2.5%","P");
  leg->AddEntry(nPartGraph[14],Form("Mult. Cut -%.02F%%",percentCutVary*100),"P");
  leg->AddEntry(nPartGraph[15],Form("Mult. Cut +%.02F%%",percentCutVary*100),"P");
  
  leg->Draw("SAME");

  canvas1->cd();
  nCollTotalSysErr->Draw("2");
  nCollGraph[0]->Draw("P");
  nCollGraph[1]->Draw("P");
  nCollGraph[2]->Draw("P");
  nCollGraph[3]->Draw("P");
  nCollGraph[4]->Draw("P");
  nCollGraph[5]->Draw("P");
  nCollGraph[6]->Draw("P");
  nCollGraph[7]->Draw("P");
  nCollGraph[8]->Draw("P");
  nCollGraph[9]->Draw("P");
  nCollGraph[10]->Draw("P");
  nCollGraph[11]->Draw("P");
  nCollGraph[12]->Draw("P");
  nCollGraph[13]->Draw("P");
  nCollGraph[14]->Draw("P");
  nCollGraph[15]->Draw("P");

  nCollGraph[0]->Draw("P");

  TPaveText *title1 = new TPaveText(.60,.83,.94,.93,"brNDC");
  title1->SetBorderSize(0);
  title1->SetFillColor(kWhite);
  title1->SetTextSize(.045);
  title1->SetTextFont(62);
  title1->SetTextAlign(31);
  title1->AddText("Contributions to N_{coll} Systematic Error");
  title1->Draw("SAME");

  canvas2->cd();
  impactParamTotalSysErr->Draw("2");
  impactParamGraph[0]->Draw("P");
  impactParamGraph[1]->Draw("P");
  impactParamGraph[2]->Draw("P");
  impactParamGraph[3]->Draw("P");
  impactParamGraph[4]->Draw("P");
  impactParamGraph[5]->Draw("P");
  impactParamGraph[6]->Draw("P");
  impactParamGraph[7]->Draw("P");
  impactParamGraph[8]->Draw("P");
  impactParamGraph[9]->Draw("P");
  impactParamGraph[10]->Draw("P");
  impactParamGraph[11]->Draw("P");
  impactParamGraph[12]->Draw("P");
  impactParamGraph[13]->Draw("P");
  impactParamGraph[14]->Draw("P");
  impactParamGraph[15]->Draw("P");

  impactParamGraph[0]->Draw("P");

  TPaveText *title2 = new TPaveText(.60,.83,.94,.93,"brNDC");
  title2->SetBorderSize(0);
  title2->SetFillColor(kWhite);
  title2->SetTextSize(.045);
  title2->SetTextFont(62);
  title2->SetTextAlign(31);
  title2->AddText("Contributions to Impact Parameter");
  title2->AddText("Systematic Error");
  title2->Draw("SAME");
  
  //Print Loop
  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){

    Double_t meanDefaultNpart = nPartGraph[0]->GetX()[iCentBin];
    Double_t meanDefaultNcoll = nCollGraph[0]->GetX()[iCentBin];
    Double_t meanDefaultImpactParam = impactParamGraph[0]->GetX()[iCentBin];
        
    Double_t percentTotalSys  = nPartTotalSysErr->GetEY()[iCentBin];
    Double_t percentNcollTotalSys = nCollTotalSysErr->GetEY()[iCentBin];
    Double_t percentImpactParamTotalSys = impactParamTotalSysErr->GetEY()[iCentBin];
    
    cout <<Form("%.2F +- %.2F (%.2F%%)\t \t %.2F +- %.2F (%.2F%%)\t \t %.2F +- %.2F (%.2F%%) ",
		meanDefaultNpart,meanDefaultNpart*percentTotalSys,percentTotalSys*100.,
		meanDefaultNcoll,meanDefaultNcoll*percentNcollTotalSys,percentNcollTotalSys*100,
		meanDefaultImpactParam,meanDefaultImpactParam*percentImpactParamTotalSys,
		percentImpactParamTotalSys*100) <<endl;

    nPartMeans.push_back(meanDefaultNpart);
    nCollMeans.push_back(meanDefaultNcoll);
    impactParamMeans.push_back(meanDefaultImpactParam);
    nPartSysErrors.push_back(meanDefaultNpart*percentTotalSys);
    nCollSysErrors.push_back(meanDefaultNcoll*percentNcollTotalSys);
    impactParamSysErrors.push_back(meanDefaultImpactParam*percentImpactParamTotalSys);

  }//End Loop Over Cent Bins

  //Clean Up some globals that won't be used again
  for (Int_t iIndex=0; iIndex<16; iIndex++){
    for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
      if (nPartHistoArray[iIndex][iCentBin]){
	delete nPartHistoArray[iIndex][iCentBin];
      }
      if (nCollHistoArray[iIndex][iCentBin]){
	delete nCollHistoArray[iIndex][iCentBin];
      }
      if (impactParamHistoArray[iIndex][iCentBin]){
	delete impactParamHistoArray[iIndex][iCentBin];
      }
    }
  }
  //Re-enable associated histograms with directories
  TH1::AddDirectory(true);
 
}

//________________________________________________________________________________________
void SystematicErrors(const Int_t nEvents, const Int_t nNucleonsA, const Int_t nNucleonsB,
		      const Double_t defaultInelSigma, NegBinomialSearchResults *fitResults,
		      TH1D *dataHisto, Bool_t useTrackingEfficiency=false){

  //A wrapper for the SystematicErrors Function above that allows passing
  //in the struct used for holding all the values

  Double_t npp = fitResults->npp;
  Double_t k   = fitResults->k;
  std::vector<double> *cuts = &fitResults->centralityBinCuts;
  const int nCentBins = cuts->size();
  
  //To determine by how much to vary the cuts we take the ratio of the integrals of the
  //top 20% of each histogram as defined by the cuts themselves 
  Double_t percentVaryCuts(0);
  Double_t simHistoIntegral =
    fitResults->bestFitHisto->Integral(fitResults->bestFitHisto->FindBin(cuts->at(13)),
				      fitResults->bestFitHisto->GetNbinsX());
  Double_t dataHistoIntegral =
    dataHisto->Integral(dataHisto->FindBin(cuts->at(13)),dataHisto->GetNbinsX());
  percentVaryCuts = fabs((simHistoIntegral/dataHistoIntegral) -1.0 );
			
  SystematicErrors(nEvents, nNucleonsA, nNucleonsB,defaultInelSigma, npp, k,cuts,
		   percentVaryCuts, useTrackingEfficiency);

  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
    //The means may have already been filled in FindNpartNcollDistributions
    if (fitResults->nPartMeans.size() == 0){
      fitResults->nPartMeans.push_back(nPartMeans.at(iCentBin));
      fitResults->nCollMeans.push_back(nCollMeans.at(iCentBin));
      fitResults->impactParamMeans.push_back(impactParamMeans.at(iCentBin));
    }
    fitResults->nPartSysErrors.push_back(nPartSysErrors.at(iCentBin));
    fitResults->nCollSysErrors.push_back(nCollSysErrors.at(iCentBin));
    fitResults->impactParamSysErrors.push_back(impactParamSysErrors.at(iCentBin));
    
  }

  fitResults->nPartTotalSysErrCanvas = (TCanvas *)canvas->Clone(canvas->GetName());
  fitResults->nCollTotalSysErrCanvas = (TCanvas *)canvas1->Clone(canvas->GetName());
  fitResults->impactParamTotalSysErrCanvas = (TCanvas *)canvas2->Clone(canvas->GetName());

  //Clean Up Globals
  delete canvas;
  delete canvas1;
  delete canvas2;

  

  
}


//________________________________________________________________________________________
void SystematicErrors_Test(){

  //A wrapper for the SystematicErrors Function above that is useful for testing

  const Int_t nCentBins(16);
  Int_t nEvents(500);
  Int_t nNucleonsA(197), nNucleonsB(197);
  Double_t defaultInelSigma=30.5;
  Double_t npp(0.31), k(1.58);
  Double_t percentVaryCuts(0.07);
  Double_t cutsArr[nCentBins] = {5,8,11,14,19,26,33,42,53,65,79,96,115,137,162,192};
  std::vector<double> cuts;
  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
    cuts.push_back(cutsArr[iCentBin]);
  }
  Bool_t useTrackingEfficiency = false;
  
  
  SystematicErrors(nEvents, nNucleonsA, nNucleonsB, defaultInelSigma, npp, k,
		   &cuts, percentVaryCuts,useTrackingEfficiency);

}
