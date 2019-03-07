#include <iostream>

#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include "globalDefinitions.h"
#include "utilityFunctions.h"
#include "ParticleInfo.h"
#include "SpectraCorrectionCurves.h"

//___________________________________________________________________
CorrectionCurves::CorrectionCurves(TFile *CORRECTIONFILE, int PID, int CHARGE, int CENTINDEX, int YINDEX){
  corrFile = CORRECTIONFILE;
  pid = PID;
  charge = CHARGE;
  iCentBin = CENTINDEX;
  yIndex = YINDEX;
  particleInfo = new ParticleInfo();

  //Set All the Pointers to NULL
  interBinMeanCurve = NULL;
  interBinMeanCurveConf = NULL;
  energyLossCurve = NULL;
  energyLossCurveConf = NULL;
  efficiencyCurve = NULL;
  efficiencyCurveConf = NULL;
  efficiencyCurveSys = NULL;
  efficiencyCurveSysConf = NULL;
  knockoutProtonCurve = NULL;
  knockoutProtonCurveConf = NULL;
  muonBackgroundCurvePowerLaw = NULL;
  muonBackgroundCurveExpo = NULL;
  feedDownBackgroundCurveExpo = NULL;
  feedDownBackgroundCurvePowerLaw = NULL;
  feedDownBackgroundCurveExpoConf = NULL;
  feedDownBackgroundCurvePowerLawConf = NULL;
  
  LoadCurves();
}

//___________________________________________________________________
CorrectionCurves::~CorrectionCurves(){

  if (particleInfo)
    delete particleInfo;
  
  if (interBinMeanCurve)
    delete interBinMeanCurve;
  if (interBinMeanCurveConf)
    delete interBinMeanCurveConf;
  if (energyLossCurve)
    delete energyLossCurve;
  if (energyLossCurveConf)
    delete energyLossCurveConf;
  if (efficiencyCurve)
    delete efficiencyCurve;
  if (efficiencyCurveConf)
    delete efficiencyCurveConf;
  if (efficiencyCurveSys)
    delete efficiencyCurveSys;
  if (efficiencyCurveSysConf)
    delete efficiencyCurveSysConf;
  if (knockoutProtonCurve)
    delete knockoutProtonCurve;
  if (knockoutProtonCurveConf)
    delete knockoutProtonCurveConf;
  if (muonBackgroundCurvePowerLaw)
    delete muonBackgroundCurvePowerLaw;
  if (muonBackgroundCurveExpo)
    delete muonBackgroundCurveExpo;
  if (muonBackgroundCurvePowerLawConf)
    delete muonBackgroundCurvePowerLawConf;
  if (muonBackgroundCurveExpoConf)
    delete muonBackgroundCurveExpoConf;
  if (feedDownBackgroundCurveExpo)
    delete feedDownBackgroundCurveExpo;
  if (feedDownBackgroundCurvePowerLaw)
    delete feedDownBackgroundCurvePowerLaw;
  if (feedDownBackgroundCurveExpoConf)
    delete feedDownBackgroundCurveExpoConf;
  if (feedDownBackgroundCurvePowerLawConf)
    delete feedDownBackgroundCurvePowerLawConf;
}

//____________________________________________________
void CorrectionCurves::LoadCurves(){
    corrFile->cd();

    //Inter-Bin Transverse Mass Average
    interBinMeanCurve =
      (TF1 *)corrFile->Get(Form("%s/InterBinmTm0Fits/avgmTm0_%s_Cent%02d_yIndex%02d_Fit",
				particleInfo->GetParticleName(pid,charge).Data(),
				particleInfo->GetParticleName(pid,charge).Data(),
				iCentBin,yIndex));
    interBinMeanCurveConf =
      (TGraphErrors *)corrFile->Get(Form("%s/InterBinmTm0Fits/avgmTm0_%s_Cent%02d_yIndex%02d_Fit_Conf",
					 particleInfo->GetParticleName(pid,charge).Data(),
					 particleInfo->GetParticleName(pid,charge).Data(),
					 iCentBin,yIndex));

    if (!interBinMeanCurve){
      cout <<"No Inter-Bin Transverse Mass Curve Available. Using bin center.\n";
      interBinMeanCurve = new TF1("tempInterBinMeanCurve","pol0",0,2);
      interBinMeanCurve->SetParameter(0,mTm0BinWidth/2.0);
    }

    //Energy Loss
    
    energyLossCurve =
      (TF1 *)corrFile->Get(Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d",
				particleInfo->GetParticleName(pid,charge).Data(),
				particleInfo->GetParticleName(pid,charge).Data(),
				yIndex));
    energyLossCurveConf =
      (TGraphErrors*)corrFile->Get(Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Conf",
					particleInfo->GetParticleName(pid,charge).Data(),
					particleInfo->GetParticleName(pid,charge).Data(),
					yIndex));
    
    
    //Tracking Efficiency
    efficiencyCurve =
      (TF1 *)corrFile->Get(Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d",
				particleInfo->GetParticleName(pid,charge).Data(),
				particleInfo->GetParticleName(pid,charge).Data(),
				iCentBin,yIndex));
    efficiencyCurveConf =
      (TGraphErrors *)corrFile->Get(Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys_Conf",
					 particleInfo->GetParticleName(pid,charge).Data(),
					 particleInfo->GetParticleName(pid,charge).Data(),
					 iCentBin,yIndex));

    efficiencyCurveSys =
      (TF1 *)corrFile->Get(Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d",
				particleInfo->GetParticleName(pid,charge).Data(),
				particleInfo->GetParticleName(pid,charge).Data(),
				iCentBin,yIndex));

    efficiencyCurveSysConf =
      (TGraphErrors *)corrFile->Get(Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys_Conf",
					 particleInfo->GetParticleName(pid,charge).Data(),
					 particleInfo->GetParticleName(pid,charge).Data(),
					 iCentBin,yIndex));

    //Proton Knockout -- For Protons Only
    if (pid == PROTON && charge > 0){
      knockoutProtonCurve = (TF1 *)corrFile->Get(Form("%s/KnockoutBackgroundFits/knockoutProton_%d_%d_Fit",
						      particleInfo->GetParticleName(pid,charge).Data(),
						      iCentBin,yIndex));
    }//End Knockout Protons
    
    //Muon Background -- For Pions Only
    if (pid == PION){
      muonBackgroundCurvePowerLaw =
	(TF1*)corrFile->Get(Form("%s/MuonBackgroundFits/muonBackground_%s_cent%02d_yIndex%02d_FitPower",
				 particleInfo->GetParticleName(pid,charge).Data(),
				 particleInfo->GetParticleName(pid,charge).Data(),
				 iCentBin,yIndex));
      muonBackgroundCurveExpo =
	(TF1*)corrFile->Get(Form("%s/MuonBackgroundFits/muonBackground_%s_cent%02d_yIndex%02d_FitExp",
				 particleInfo->GetParticleName(pid,charge).Data(),
				 particleInfo->GetParticleName(pid,charge).Data(),
				 iCentBin,yIndex));
      
      muonBackgroundCurvePowerLawConf =
	(TGraphErrors*)corrFile->Get(Form("%s/MuonBackgroundFits/muonBackground_%s_cent%02d_yIndex%02d_FitPower_Conf",
					  particleInfo->GetParticleName(pid,charge).Data(),
					  particleInfo->GetParticleName(pid,charge).Data(),
					  iCentBin,yIndex));
      muonBackgroundCurveExpoConf =
	(TGraphErrors*)corrFile->Get(Form("%s/MuonBackgroundFits/muonBackground_%s_cent%02d_yIndex%02d_FitExp_Conf",
					  particleInfo->GetParticleName(pid,charge).Data(),
					  particleInfo->GetParticleName(pid,charge).Data(),
					  iCentBin,yIndex));
      
    }//End Muon Background

      
    //Feed Down - For Pions and Protons
    if (pid == PION || pid == PROTON){

     
	feedDownBackgroundCurveExpo =
	 	  (TF1 *)corrFile->Get(Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp",
				    particleInfo->GetParticleName(pid,charge).Data(),
				    particleInfo->GetParticleName(pid,charge).Data(),
				    0,yIndex));
	feedDownBackgroundCurvePowerLaw =
	  (TF1 *)corrFile->Get(Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitPower",
				    particleInfo->GetParticleName(pid,charge).Data(),
				    particleInfo->GetParticleName(pid,charge).Data(),
				    0,yIndex));
	
	feedDownBackgroundCurveExpoConf =
	  (TGraphErrors *)corrFile->Get(Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp_Conf",
					     particleInfo->GetParticleName(pid,charge).Data(),
					     particleInfo->GetParticleName(pid,charge).Data(),
					     0,yIndex));
	feedDownBackgroundCurvePowerLawConf =
	  (TGraphErrors *)corrFile->Get(Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitPower_Conf",
					     particleInfo->GetParticleName(pid,charge).Data(),
					     particleInfo->GetParticleName(pid,charge).Data(),
					     0,yIndex));
      

					     /*
	feedDownBackgroundCurveExpo =
	 	  (TF1 *)corrFile->Get(Form("%s/HadronicBackgroundFits/hadronicBackground_%s_cent%02d_yIndex%02d_Fit",
				    particleInfo->GetParticleName(pid,charge).Data(),
				    particleInfo->GetParticleName(pid,charge).Data(),
				    0,yIndex));
	feedDownBackgroundCurvePowerLaw =
	  (TF1 *)corrFile->Get(Form("%s/HadronicBackgroundFits/hadronicBackground_%s_cent%02d_yIndex%02d_Fit",
				    particleInfo->GetParticleName(pid,charge).Data(),
				    particleInfo->GetParticleName(pid,charge).Data(),
				    0,yIndex));
	
	feedDownBackgroundCurveExpoConf =
	  (TGraphErrors *)corrFile->Get(Form("%s/HadronicBackgroundFits/hadronicBackground_%s_cent%02d_yIndex%02d_Fit_Conf",
					     particleInfo->GetParticleName(pid,charge).Data(),
					     particleInfo->GetParticleName(pid,charge).Data(),
					     0,yIndex));
	feedDownBackgroundCurvePowerLawConf =
	  (TGraphErrors *)corrFile->Get(Form("%s/HadronicBackgroundFits/hadronicBackground_%s_cent%02d_yIndex%02d_Fit_Conf",
					     particleInfo->GetParticleName(pid,charge).Data(),
					     particleInfo->GetParticleName(pid,charge).Data(),
					     0,yIndex));

      */

    }//End Feed Down
      
}
    
