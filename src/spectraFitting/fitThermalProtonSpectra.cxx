#include <iostream>
#include <vector>
#include <utility>
#include <numeric>

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
#include <TThread.h>

#include "globalDefinitions.h"
#include "utilityFunctions.h"
#include "StRefMultExtendedCorr.h"
#include "TFitResult.h"
#include "UserCuts.h"
#include "ParticleInfo.h"
#include "SpectraClass.h"
#include "SpectraFitUtilities.h"
#include "SpectraFitFunctions.h"
#include "SimFitter.h"
#include "StyleSettings.h"


std::pair<double,double> computeTotaldNdy(TF1 *thisFunc, TMatrixDSym *thisCovMatrix, double min, double max);
void removeTpcPoints(TGraphErrors * spectraGraph, Int_t yIndex);

void fitThermalProtonSpectra(TString spectraFileName,TString outfile, TString imgDir, int centBin= -999){

  //Draw and Save options
  bool draw = true;
  bool save = true;
  bool reflect = false;
  //constants 
  ParticleInfo *particleInfo = new ParticleInfo();
  const int nCentBins(GetNCentralityBins());
  const int pid = PROTON;
  const int charge = 1;

  //Directories 
  TFile *spectraFile = new TFile(spectraFileName,save ? "UPDATE" : "READ");
  TString saveDir = TString::Format("FitSpectraClass_%s",particleInfo->GetParticleName(pid).Data());
  TString saveDirFit = TString::Format("%s_FitFuncs",saveDir.Data());

  TFile *outFile = new TFile(outfile,"RECREATE");
  TString outDir = TString::Format("dNdyGraph_%s",particleInfo->GetParticleName(pid,charge).Data());
  outFile->mkdir(Form("%s",outDir.Data()));
  outFile->mkdir(Form("CorrectedSpectra_%s",particleInfo->GetParticleName(pid,charge).Data()));
  outFile->mkdir(Form("BlastWaveFit_%s",particleInfo->GetParticleName(pid,charge).Data()));
  outFile->mkdir(Form("ThermalFit_%s",particleInfo->GetParticleName(pid,charge).Data()));

  //Initialize Draw Canvas
  TCanvas * spectraCanvas = NULL;
  TCanvas * spectraSysCanvas = NULL;

  TCanvas * parameterCanvas = NULL;
  TCanvas * parameterSysCanvas = NULL;

  spectraCanvas = new TCanvas("spectraCanvas","spectraCanvas",1.3*700,700);
  spectraSysCanvas = new TCanvas("spectraSysCanvas","spectraSysCanvas",1.3*700,700);
  parameterCanvas = new TCanvas("parameterCanvas","parameterCanvas",1000,500);
  parameterSysCanvas = new TCanvas("parameterSysCanvas","parameterSysCanvas",900,900);

  //Divide specta canvas by number of fits and parameter canvas by number of varying parameters
  spectraCanvas->Divide(4,3);
  spectraSysCanvas->Divide(4,3);
  parameterCanvas->Divide(2,1);
  parameterSysCanvas->Divide(2,2);

  //set nominal and systematic fit
  //******************Normal funcs***********************
  TF1 *nominalFit = new TF1("nominalFit",ThermalFitFunc,0.001,2.0,3);
  TF1 *sysFit = new TF1("sysFit",BlastWaveModelFit,0.0,2.0,6);

  //******************testing***********************
  nominalFit->SetNpx(500);

  //enum the thermal fit parameters
  enum simFitPars {AMP_TH,TSLOPE_TH,PMASS_TH};

  //define parameter graph vectors 
   std::vector<std::vector< TGraphErrors *> >  parameterGraphs
    (nCentBins, std::vector<TGraphErrors *> 
     (2, (TGraphErrors *) NULL));

  std::vector<std::vector< TGraphErrors *> >  parameterGraphsBW
    (nCentBins, std::vector<TGraphErrors *> 
     (4, (TGraphErrors *) NULL));

  std::vector<std::vector< TGraphErrors *> >  dNdyGraph
    (nCentBins, std::vector<TGraphErrors *> 
     (2, (TGraphErrors *) NULL));

  std::vector<std::vector< TGraphErrors *> >  dNdyGraphReflect
    (nCentBins, std::vector<TGraphErrors *> 
     (2, (TGraphErrors *) NULL));

  //make parameter graphs

  for (int iCentBin=0; iCentBin<nCentBins; iCentBin++){
    
    //if user is only running one centrality bin
    if ( !(centBin == -999 || iCentBin == centBin) ) continue;  
    
    for (int i = 0;i<2;i++){
      parameterGraphs.at(iCentBin).at(i) = new TGraphErrors();
      parameterGraphs.at(iCentBin).at(i)->SetName(Form("parameterGraph_Cent%02d_par%d",iCentBin,i));
      parameterGraphs.at(iCentBin).at(i)->SetTitle(Form("parameterGraph_Cent%02d_par%d",iCentBin,i));
    }
    
    for (int i = 0;i<4;i++){
      parameterGraphsBW.at(iCentBin).at(i) = new TGraphErrors();
      parameterGraphsBW.at(iCentBin).at(i)->SetName(Form("parameterGraph_Cent%02d_par%d",iCentBin,i));
      parameterGraphsBW.at(iCentBin).at(i)->SetTitle(Form("parameterGraph_Cent%02d_par%d",iCentBin,i));
    }
  } 

  //define spectra vectors and spectra graph vectors                                  
  std::vector<std::vector<SpectraClass *> >  spectra
    (nCentBins,std::vector<SpectraClass *> 
     (0,(SpectraClass *)NULL));
  std::vector<std::vector<TGraphErrors *> > spectraToFit
    (nCentBins,std::vector<TGraphErrors *> 
     (0,(TGraphErrors *)NULL));
  std::vector<std::vector<TF1 *> > fitFuncs
    (nCentBins,std::vector<TF1 *> 
     (0, (TF1 *)NULL));
  std::vector<std::vector<TF1 *> >  sysFuncs
    (nCentBins,std::vector<TF1 *> 
     (0, (TF1 *)NULL));
  std::vector<std::vector<double> >  rapidityVec
    (nCentBins,std::vector<double> 
     (0, (double)0));


  //Load spectra and make fitFuncs
  for (int iCentBin=0; iCentBin<nCentBins; iCentBin++){
    if ( !(centBin == -999 || iCentBin == centBin) ) continue;  
    for (unsigned int yIndex=0; yIndex<=18; yIndex++){  

      //if user is only running one centrality bin

      
	    TString spectraName = TString::Format("CorrectedSpectra_%s_Cent%02d_yIndex%02d_Total",
						  particleInfo->GetParticleName(pid,charge).Data(),
						  iCentBin,yIndex);
	    TString spectraClassDir = TString::Format("SpectraClass_%s",
						      particleInfo->GetParticleName(pid,charge).Data());
	    TString spectrumDir = TString::Format("CorrectedSpectra_%s",
						  particleInfo->GetParticleName(pid,charge).Data());
	    TString spectrumFitDir = TString::Format("%s_FitFuncs",
						     spectrumDir.Data());
	    TString dNdyDir = TString::Format("dNdyGraph_%s",
					      particleInfo->GetParticleName(pid,charge).Data());
	    spectraFile->cd();
	    spectraFile->cd(spectraClassDir.Data());
	    SpectraClass *temp = (SpectraClass *)gDirectory->Get(spectraName.Data());
	    
	    if (temp){

	      cout << "Loaded iCent " << iCentBin << " and yIndex " << yIndex << endl;

	      rapidityVec.at(iCentBin).push_back(GetRapidityRangeCenter(yIndex));
	      spectra.at(iCentBin).push_back(temp);
	      spectra.at(iCentBin).back()->SetSpectraFile(spectraFile);
	      spectra.at(iCentBin).back()->SetSpectraClassDir(spectraClassDir);
	      spectra.at(iCentBin).back()->SetSpectraDir(spectrumDir);
	      spectra.at(iCentBin).back()->SetSpectraNamePrefix("correctedSpectra");
	      spectra.at(iCentBin).back()->LoadSpectra();
	      
	      spectraToFit.at(iCentBin).push_back(spectra.at(iCentBin).back()->GetTGraph());
	      //	      removeTpcPoints(spectraToFit.at(iCentBin).back(),yIndex);

	      fitFuncs.at(iCentBin).push_back(new TF1(*nominalFit));
	      sysFuncs.at(iCentBin).push_back(new TF1(*sysFit));
	      sysFuncs.at(iCentBin).back()->SetLineStyle(9);
	      
	      fitFuncs.at(iCentBin).back()->SetName(Form("%s_Fit0",spectraName.Data()));
	      sysFuncs.at(iCentBin).back()->SetName(Form("%s_Fit1",spectraName.Data()));

	    } //Enod Loop Over temp
    }//End Loop Over yIndex
  //*************************
  //***** NOMINAL FIT *******
  //*************************

    //Seeds for thermal fit      
    Float_t ampThermal = 40;
    Float_t tSlope = .19;
    Float_t pMass = particleInfo->GetParticleMass(pid);

      //Loop through the yIndex
    for (std::vector<int>::size_type index = 0; index != fitFuncs.at(iCentBin).size(); index++){
      
      Double_t rapidity = rapidityVec.at(iCentBin).at(index);
      //      rapidity = (-1.0)*rapidity;
      TF1 * spectraFuncTemp = fitFuncs.at(iCentBin).at(index);
      TGraphErrors * tempGraphToFit = spectraToFit.at(iCentBin).at(index);

      //Set function
      spectraFuncTemp->SetParameter(AMP_TH,ampThermal);
      spectraFuncTemp->SetParameter(TSLOPE_TH,tSlope);
      spectraFuncTemp->FixParameter(PMASS_TH,pMass);

      if ( abs(rapidity) > 1.8 ) {
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ high rap fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	tempGraphToFit->Fit(spectraFuncTemp,"","",0.4,1);
      }
      else {
	tempGraphToFit->Fit(spectraFuncTemp,"","",0.2,1);
      }

      //Add fit results to parameter graphs
      parameterGraphs.at(iCentBin).at(0)->SetPoint(parameterGraphs.at(iCentBin).at(0)->GetN(),rapidity,spectraFuncTemp->GetParameter(0));
      parameterGraphs.at(iCentBin).at(0)->SetPointError(parameterGraphs.at(iCentBin).at(0)->GetN()-1,0,spectraFuncTemp->GetParError(0));
      parameterGraphs.at(iCentBin).at(1)->SetPoint(parameterGraphs.at(iCentBin).at(1)->GetN(),rapidity,spectraFuncTemp->GetParameter(1));	
      parameterGraphs.at(iCentBin).at(1)->SetPointError(parameterGraphs.at(iCentBin).at(1)->GetN()-1,0,spectraFuncTemp->GetParError(1));
    }
  
    //If option set, draw parameters
    if (draw){
      
      parameterCanvas->cd(1);
      parameterGraphs.at(iCentBin).at(0)->SetMarkerStyle(20);
      parameterGraphs.at(iCentBin).at(0)->Draw("AP");
      parameterCanvas->cd(2);
      parameterGraphs.at(iCentBin).at(1)->SetMarkerStyle(20);
      parameterGraphs.at(iCentBin).at(1)->Draw("AP");
      parameterCanvas->SaveAs(Form("%s/thermalParameters_cent_%02d.png",imgDir.Data(),iCentBin));

    }

      
  //*************************
  //***** SYSTEMATIC FIT *******
  //*************************

    //Define initial set values
    Float_t ampBlastwave = 50000;
    Float_t tempBlastwave = .100;
    Float_t surfaceVelocity = .1; 
    Float_t transVelocityPower = .2;
    Float_t radialIntegral = 1;

    //Make temperature vector for later averaging 
    std::vector <Float_t> FitTemp;

    //enum the three fit stages
    enum FitStage{FREEFIT,FIXTEMP,FINALFIT};

    //enum the six blastwave fit parameters
    enum sysFitPars {AMP,TEMP,SURFACEVELOCITY,TRANSVELOCITYPOWER,
                   RADIALINTEGRALLIMIT,PMASS};
      
    //make function for fitting the surface velocity
    TF1 * surfaceVelocityFit = new TF1("","[0]*cos((x - [2])*[1] )",-3,3);
    
    TCanvas * canTest = new TCanvas("c1","c1",800,600);
    
     for (int iFit=0; iFit <=FINALFIT; iFit++){
       

       //before the final fit, fit the surface velocity graph with a gaussian 
       if (iFit == FINALFIT)	       
	 {

	   //Set intial values for surface velocity fit; Fix centroid at mid-rapidity
	   surfaceVelocityFit->SetParameter(0,0.5);
	   surfaceVelocityFit->SetParameter(1,.5);
	   surfaceVelocityFit->FixParameter(2,-1.5);
	 
	   //Fit surface velocity Gaussian


	   parameterGraphsBW.at(iCentBin).at(3)->Fit(surfaceVelocityFit,"R");

	   canTest->cd();
	   parameterGraphsBW.at(iCentBin).at(3)->Draw();
	   
	   //Save Surface velocityFit
	   canTest->SaveAs(Form("%s/surfaceVelocityFit.png",imgDir.Data(),iCentBin));

	   //Reinitialize the parameterGraphsBw for surface velocity
	   parameterGraphsBW.at(iCentBin).at(3) = new TGraphErrors();
	   parameterGraphsBW.at(iCentBin).at(3)->SetName(Form("parameterGraph_Cent%02d_par%d",iCentBin,2));
	   parameterGraphsBW.at(iCentBin).at(3)->SetTitle(Form("parameterGraph_Cent%02d_par%d",iCentBin,2));

	 }
       
       for (std::vector<int>::size_type index = 0; index != sysFuncs.at(iCentBin).size(); index++){
	 
	 //Get the relative rapidity for the sysFunction 
	 Double_t rapidity = rapidityVec.at(iCentBin).at(index);

	 //	 rapidity = rapidity*(-1.0);

	 //Pass sysFunc to local pointer
	 TF1 * spectraFuncTemp = sysFuncs.at(iCentBin).at(index);

	 //Pass Graph to local pointer	 
	 TGraphErrors * tempGraphToFit = spectraToFit.at(iCentBin).at(index);
	 Float_t averageTemp = ( accumulate(FitTemp.begin(), FitTemp.end(), 0.0 ) / FitTemp.size() );
	 spectraFuncTemp->SetNpx(100);
	 
	 //Set function
	 spectraFuncTemp->FixParameter(TRANSVELOCITYPOWER,transVelocityPower);
	 spectraFuncTemp->FixParameter(RADIALINTEGRALLIMIT,radialIntegral);
	 spectraFuncTemp->FixParameter(PMASS,pMass);
	 
	 //Set fit initial values for different fit stage 
	 switch(iFit)
	   {
	   case FREEFIT: //Unconstrained     
	     spectraFuncTemp->SetParameter(AMP,ampBlastwave);
	     spectraFuncTemp->SetParameter(SURFACEVELOCITY,surfaceVelocity);
	     spectraFuncTemp->SetParLimits(SURFACEVELOCITY,.01,1.0);
	     spectraFuncTemp->SetParameter(TEMP,tempBlastwave);
	     spectraFuncTemp->SetParLimits(TEMP,0.020,.160);
	     break;
	     
	   case FIXTEMP: //Fix Temp to average value of FREEFIT
	     spectraFuncTemp->SetParameter(AMP,ampBlastwave);
	     spectraFuncTemp->SetParameter(SURFACEVELOCITY,spectraFuncTemp->GetParameter(SURFACEVELOCITY));
	     spectraFuncTemp->FixParameter(TEMP,averageTemp);
	     break;
	     
	   case FINALFIT: //Fix Temp to average value of FREEFIT and Fix surfaceVelocity to value from surfaceVelocity fit
	     spectraFuncTemp->SetParameter(AMP,spectraFuncTemp->GetParameter(AMP));
	     spectraFuncTemp->FixParameter(SURFACEVELOCITY,surfaceVelocityFit->Eval(rapidity));
	     spectraFuncTemp->FixParameter(TEMP,spectraFuncTemp->GetParameter(TEMP));
	     break;
	   }	  
	 
	 TFitResultPtr blastWaveFit = tempGraphToFit->Fit(spectraFuncTemp,"R S");
	 
	 //Add temp to temp vector 
	 FitTemp.push_back(spectraFuncTemp->GetParameter(TEMP));
	 
	 //sum blastwave fit to get total dndny
	 TMatrixDSym blastWaveCov = blastWaveFit->GetCovarianceMatrix();
	 std::pair<double,double> blastWavedNdy = computeTotaldNdy(spectraFuncTemp,&blastWaveCov,0.001,3);

	 //Set points on parameter graph for final round fitting 
	 if (iFit == FIXTEMP){
	   parameterGraphsBW.at(iCentBin).at(3)->SetPoint(parameterGraphsBW.at(iCentBin).at(3)->GetN(),rapidity,spectraFuncTemp->GetParameter(2));
	   parameterGraphsBW.at(iCentBin).at(3)->SetPointError(parameterGraphsBW.at(iCentBin).at(3)->GetN()-1,0,spectraFuncTemp->GetParError(2));
	 }
	 
	 if (iFit == FINALFIT){

	   //Add fit results to parameter graphs
	   parameterGraphsBW.at(iCentBin).at(0)->SetPoint(parameterGraphsBW.at(iCentBin).at(0)->GetN(),rapidity,spectraFuncTemp->GetParameter(0));
	   parameterGraphsBW.at(iCentBin).at(0)->SetPointError(parameterGraphsBW.at(iCentBin).at(0)->GetN()-1,0,spectraFuncTemp->GetParError(0));
	   
	   parameterGraphsBW.at(iCentBin).at(1)->SetPoint(parameterGraphsBW.at(iCentBin).at(1)->GetN(),rapidity,blastWavedNdy.first);
	   parameterGraphsBW.at(iCentBin).at(1)->SetPointError(parameterGraphsBW.at(iCentBin).at(1)->GetN()-1,0,blastWavedNdy.second);
	   
	   parameterGraphsBW.at(iCentBin).at(2)->SetPoint(parameterGraphsBW.at(iCentBin).at(2)->GetN(),rapidity,spectraFuncTemp->GetParameter(1));	
	   parameterGraphsBW.at(iCentBin).at(2)->SetPointError(parameterGraphsBW.at(iCentBin).at(2)->GetN()-1,0,spectraFuncTemp->GetParError(1));
	   
	   parameterGraphsBW.at(iCentBin).at(3)->SetPoint(parameterGraphsBW.at(iCentBin).at(3)->GetN(),rapidity,spectraFuncTemp->GetParameter(2));
	   parameterGraphsBW.at(iCentBin).at(3)->SetPointError(parameterGraphsBW.at(iCentBin).at(3)->GetN()-1,0,spectraFuncTemp->GetParError(2));
	 }
       }

       //Only Draw if Final fit
       if (iFit == FINALFIT){

	 //Check draw option
	 if (draw){

	   //Draw each parameter
	   parameterSysCanvas->cd(1);
	   parameterGraphsBW.at(iCentBin).at(0)->SetMarkerStyle(20);
	   parameterGraphsBW.at(iCentBin).at(0)->Draw("AP");
	   cout << "Normalization" << endl;
	   parameterGraphsBW.at(iCentBin).at(0)->Print();
	   parameterSysCanvas->cd(2);
	   parameterGraphsBW.at(iCentBin).at(1)->SetMarkerStyle(20);
	   parameterGraphsBW.at(iCentBin).at(1)->Draw("AP");
	   cout << "dNdy" << endl;
	   parameterGraphsBW.at(iCentBin).at(1)->Print();
	   parameterSysCanvas->cd(3);
	   parameterGraphsBW.at(iCentBin).at(2)->SetMarkerStyle(20);
	   parameterGraphsBW.at(iCentBin).at(2)->Draw("AP");
	   cout << "Temp" << endl;
	   parameterGraphsBW.at(iCentBin).at(2)->Print();
	   parameterSysCanvas->cd(4);
	   parameterGraphsBW.at(iCentBin).at(3)->SetMarkerStyle(20);
	   parameterGraphsBW.at(iCentBin).at(3)->Draw("AP");
	   cout << "Velocity" << endl;
	   parameterGraphsBW.at(iCentBin).at(3)->Print();
	   
	   //Save parameter canvas
	   parameterSysCanvas->SaveAs(Form("%s/blastwaveParameters_cent_%02d.png",imgDir.Data(),iCentBin));
	 }
       }
     }
   
     //Pass to local pointer for Systematic Error estimation
     TGraphErrors * dNdyGraphFit1 = parameterGraphs.at(iCentBin).at(0);
     TGraphErrors * dNdyGraphFit2 = parameterGraphsBW.at(iCentBin).at(1);
     TGraphErrors * dNdyGraphStat = dNdyGraphFit1;
     TGraphErrors * dNdyGraphSys = new TGraphErrors();

     //Difference Squared of dNdy will be passed to vector
     std::vector <Float_t> ErrorDiffSq;

     //Check if number of points is different
     if (dNdyGraphFit2->GetN() != dNdyGraphFit1->GetN()) 
       {
	 cout << "Sys graph and Nominal Graph have different number of points!!!! Exiting... " << endl;
	 break;
       }

     //Loop through sys and nominal graph and take the difference if yield
     for ( int i=0; i < dNdyGraphFit1->GetN(); i++)
       {
	 Float_t diff = dNdyGraphFit2->GetY()[i] - dNdyGraphFit1->GetY()[i];
	 ErrorDiffSq.push_back( diff*diff ); 

       }

     //average the squared difference
     Float_t meanSquare = ( accumulate(ErrorDiffSq.begin(), ErrorDiffSq.end(), 0.0 ) / ErrorDiffSq.size() );

     //compute RMS
     Float_t ErrorRMS = sqrt(meanSquare);

     //Set ErrorRMS to Error on dNdy
     for ( int i=0; i < dNdyGraphFit1->GetN(); i++)
       {
	 dNdyGraphSys->SetPoint(i,dNdyGraphFit1->GetX()[i],dNdyGraphFit1->GetY()[i]);
	 dNdyGraphSys->SetPointError(i,0,ErrorRMS);

       }


     for (std::vector<int>::size_type index = 0; index != sysFuncs.at(iCentBin).size(); index++){
       
       outFile->cd();
       outFile->cd(Form("BlastWaveFit_%s",particleInfo->GetParticleName(pid,charge).Data()));
       sysFuncs.at(iCentBin).at(index)->Write();
       outFile->cd();
       outFile->cd(Form("ThermalFit_%s",particleInfo->GetParticleName(pid,charge).Data()));
       fitFuncs.at(iCentBin).at(index)->Write();

       outFile->cd(Form("CorrectedSpectra_%s",particleInfo->GetParticleName(pid,charge).Data()));
	
	if (spectraToFit.at(iCentBin).at(index)) spectraToFit.at(iCentBin).at(index)->Write();
	
       
     }


     //Save dNdySys and dNdyStat to file

     if ( reflect == true ){ 

	 dNdyGraphReflect.at(iCentBin).at(0) = new TGraphErrors();
	 dNdyGraph.at(iCentBin).at(0) = new TGraphErrors();
	 dNdyGraphReflect.at(iCentBin).at(1) = new TGraphErrors();
	 dNdyGraph.at(iCentBin).at(1) = new TGraphErrors();

       for (int i=0; i < dNdyGraphStat->GetN();i++){
	 dNdyGraphReflect.at(iCentBin).at(0)->SetPoint(i,-(dNdyGraphStat->GetX()[i]+1.5),dNdyGraphStat->GetY()[i]);
	 dNdyGraphReflect.at(iCentBin).at(0)->SetPointError(i,0,(dNdyGraphStat->GetErrorY(i)));

	 dNdyGraphReflect.at(iCentBin).at(1)->SetPoint(i,-(dNdyGraphSys->GetX()[i]+1.5),dNdyGraphSys->GetY()[i]);
	 dNdyGraphReflect.at(iCentBin).at(1)->SetPointError(i,0,(dNdyGraphSys->GetErrorY(i)));

	 dNdyGraph.at(iCentBin).at(0)->SetPoint(i,(dNdyGraphStat->GetX()[i]+1.5),dNdyGraphStat->GetY()[i]);
	 dNdyGraph.at(iCentBin).at(0)->SetPointError(i,0,(dNdyGraphStat->GetErrorY(i)));

	 dNdyGraph.at(iCentBin).at(1)->SetPoint(i,(dNdyGraphSys->GetX()[i]+1.5),dNdyGraphSys->GetY()[i]);
	 dNdyGraph.at(iCentBin).at(1)->SetPointError(i,0,(dNdyGraphSys->GetErrorY(i)));

       }

     }
     else {

       dNdyGraph.at(iCentBin).at(0) = dNdyGraphStat; 
       dNdyGraph.at(iCentBin).at(1) = dNdyGraphSys; 

     }

     dNdyGraph.at(iCentBin).at(0)->SetName(Form("dNdyGraphStat_Cent%02d",iCentBin)); 
     dNdyGraph.at(iCentBin).at(0)->SetTitle(Form("dNdyGraphStat_Cent%02d",iCentBin)); 

     dNdyGraph.at(iCentBin).at(1)->SetName(Form("dNdyGraphSys_Cent%02d",iCentBin)); 
     dNdyGraph.at(iCentBin).at(1)->SetTitle(Form("dNdyGraphSys_Cent%02d",iCentBin)); 
     
     if ( reflect == true ){ 

     dNdyGraphReflect.at(iCentBin).at(0)->SetName(Form("dNdyGraphStatRef_Cent%02d",iCentBin)); 
     dNdyGraphReflect.at(iCentBin).at(0)->SetTitle(Form("dNdyGraphStatRef_Cent%02d",iCentBin)); 

     dNdyGraphReflect.at(iCentBin).at(1)->SetName(Form("dNdyGraphSysRef_Cent%02d",iCentBin)); 
     dNdyGraphReflect.at(iCentBin).at(1)->SetTitle(Form("dNdyGraphSysRef_Cent%02d",iCentBin)); 
     
     }

     outFile->cd();
     outFile->cd( Form("%s",outDir.Data()) );
     dNdyGraph.at(iCentBin).at(0)->Write();
     dNdyGraph.at(iCentBin).at(1)->Write();

     if ( reflect == true ){ 
       dNdyGraphReflect.at(iCentBin).at(0)->Write();
       dNdyGraphReflect.at(iCentBin).at(1)->Write();
     }

  } //End of Centrality Loop
}



/******************************************************************
   INTEGRATE THE FIT CURVE TO FIND THE DNDY
 ******************************************************************/
std::pair<double,double> computeTotaldNdy(TF1 *thisFunc, TMatrixDSym *thisCovMatrix,
                                                        double min, double max){

  if (!thisFunc){
    cout <<"WARNING - SpectraClass::ComputeTotaldNdy() - The Function does not exist. Returning zeros!\n";
    return std::make_pair(0,0);
  }

  ParticleInfo * particleInfo = new ParticleInfo("SL16a",true);

  Double_t functionRangeMin(0), functionRangeMax(0);
  thisFunc->GetRange(functionRangeMin,functionRangeMax);
  Double_t stepWidth = fabs(functionRangeMax-functionRangeMin)/
    (double) thisFunc->GetNpx();

  Double_t mTm0 = min;
  Double_t sum(0);
  Double_t sumErr(0);
  Double_t pMass = particleInfo->GetParticleMass(2);
  Double_t percent;


  while (mTm0 <= max){

    Double_t mT = mTm0 + pMass;
    sum += mT * thisFunc->Integral(mTm0,mTm0+stepWidth);

    if (thisCovMatrix)
      sumErr += TMath::Power(mT*thisFunc->IntegralError(mTm0,mTm0+stepWidth,
    							thisFunc->GetParameters(),
							thisCovMatrix->GetMatrixArray()),2);
    
    mTm0 += stepWidth;
    percent = 100*mTm0/functionRangeMax;

  }

  Double_t preFactors = 2.0*TMath::Pi();

  Double_t dNdy = preFactors * sum;
  Double_t dNdyErr = preFactors * sqrt(sumErr);

  return std::make_pair(dNdy,dNdyErr);


}

/******************************************************************
   REMOVE LOW MTM0 POINTS WITH POOR EFFICIENCY CALCULATION
 ******************************************************************/
void removeTpcPoints(TGraphErrors * spectraGraph, Int_t yIndex){
  
    if (yIndex  == 8) {
     for (int n =0 ; n < 1; n++){
       spectraGraph->RemovePoint(0);
      }
   }
  
  if (yIndex == 17) {
    for (int n =0 ; n < 3; n++){
      spectraGraph->RemovePoint(spectraGraph->GetN()-1);
    }
  }


  if (yIndex == 18) {
    for (int n =0 ; n < 5; n++){
      spectraGraph->RemovePoint(spectraGraph->GetN()-1);
    }
  }

}

