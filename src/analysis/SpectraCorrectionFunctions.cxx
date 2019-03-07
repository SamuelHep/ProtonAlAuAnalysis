#include <iostream>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include "globalDefinitions.h"
#include "utilityFunctions.h"
#include "ParticleInfo.h"
#include "SpectraCorrectionCurves.h"
#include "SpectraCorrectionFunctions.h"



//________________________________________________________________
TGraphErrors *ApplyCorrectionChain(TGraphErrors *rawSpectrum, CorrectionCurves *corrCurve,
				   Double_t minTPCEfficiency, Bool_t doSystematics){

  //Apply all of the corrections below to the raw spectrum using the
  //the corrections curves in the CorrectionCurves object

  //If a correction does not apply to the present particle species then
  //the correction is skipped and the spectra is simply copied and propagated
  //to the next step in the chain.

  //If a correction curve does not exist then the correction is skipped

  //Make the raw spectrum exists
  if (!rawSpectrum){
    cout <<"WARNING - SpectraCorrectionFunctions::ApplyCorrectionChain() reports that the raw spectrum does not exist! Returning NULL! \n";
    return NULL;
  }

  if (!corrCurve){
    cout <<"WARNING - SpectraCorrectionFunctions::ApplyCorrectionChain() reports that the CorrectionCurve pointer is NULL. Returning NULL \n";
    return NULL;
  }
  
  //Get the Particle Mass
  ParticleInfo particleInfo;
  Double_t pidMass = particleInfo.GetParticleMass(corrCurve->GetPid());
  
  //Create the pointers for the Spectra at various stages
  TGraphErrors *interBinMeanSpectrum = NULL;
  TGraphErrors *muonCorrSpectrum = NULL;
  TGraphErrors *feedDownCorrSpectrum = NULL;
  TGraphErrors *finiteBinWidthCorrSpectrum = NULL;
  TGraphErrors *energyLossCorrSpectrum = NULL;
  TGraphErrors *efficiencyCorrSpectrum = NULL;

  //----NOMINAL CORRECTION MODE----
  if (!doSystematics){

    //Inter-Bin Mean Location
    interBinMeanSpectrum =
      ApplyInterBinTransverseMassLocation(rawSpectrum,corrCurve->GetInterBinMeanCurve());

    //Muon Background Correction
    if (corrCurve->GetPid() == PION)
      muonCorrSpectrum =
	ApplyMuonBackgroundCorrection(interBinMeanSpectrum,corrCurve->GetMuonBackgroundCurvePowerLaw());
    else
      muonCorrSpectrum = new TGraphErrors(*interBinMeanSpectrum);

    //Feed Down Correction
    if (corrCurve->GetPid() == PION || corrCurve->GetPid() == PROTON)
     feedDownCorrSpectrum =
    	ApplyFeedDownBackgroundCorrection(muonCorrSpectrum,corrCurve->GetFeedDownBackgroundCurveExpo());
    else
    feedDownCorrSpectrum = new TGraphErrors(*muonCorrSpectrum);

    //Finite Bin Width
    finiteBinWidthCorrSpectrum =
           ApplyFiniteBinWidthCorrection(feedDownCorrSpectrum,corrCurve->GetEnergyLossCurve(),pidMass);

    //Energy Loss
    energyLossCorrSpectrum =
      ApplyEnergyLossCorrection(finiteBinWidthCorrSpectrum,corrCurve->GetEnergyLossCurve(),pidMass);

    //Efficiency
    efficiencyCorrSpectrum =
      ApplyEfficiencyCorrection(energyLossCorrSpectrum,corrCurve->GetEfficiencyCurve(),minTPCEfficiency);

     
  }//End Nominal Correction Mode

  //----SYSTEMATIC ERROR MODE-----
  else{

    //Inter-Bin Mean Location
    interBinMeanSpectrum =
      ApplyInterBinTransverseMassLocation(rawSpectrum,corrCurve->GetInterBinMeanCurve(),
					  corrCurve->GetInterBinMeanCurveConf());

    //Muon Background Correction
    if (corrCurve->GetPid() == PION)
    muonCorrSpectrum =
      ApplyMuonBackgroundCorrection(interBinMeanSpectrum,corrCurve->GetMuonBackgroundCurvePowerLaw(),
				    corrCurve->GetMuonBackgroundCurveExpo(),
				    corrCurve->GetMuonBackgroundCurvePowerLawConf(),
				    corrCurve->GetMuonBackgroundCurveExpoConf());
    else
      muonCorrSpectrum = new TGraphErrors(*interBinMeanSpectrum);

    //Feed Down Correction
    if (corrCurve->GetPid() == PION || corrCurve->GetPid() == PROTON)
    feedDownCorrSpectrum =
    ApplyFeedDownBackgroundCorrection(muonCorrSpectrum,corrCurve->GetFeedDownBackgroundCurveExpo(),
    	corrCurve->GetFeedDownBackgroundCurvePowerLaw(),
    					corrCurve->GetFeedDownBackgroundCurveExpoConf(),
    					corrCurve->GetFeedDownBackgroundCurvePowerLawConf());
  else
    feedDownCorrSpectrum = new TGraphErrors(*muonCorrSpectrum);


    //Finite Bin Width Correction
    finiteBinWidthCorrSpectrum =
        ApplyFiniteBinWidthCorrection(feedDownCorrSpectrum,corrCurve->GetEnergyLossCurve(),pidMass,
				    			    corrCurve->GetEnergyLossCurveConf());


    //Energy Loss Correction
    energyLossCorrSpectrum =
      ApplyEnergyLossCorrection(finiteBinWidthCorrSpectrum,corrCurve->GetEnergyLossCurve(),pidMass,
				corrCurve->GetEnergyLossCurveConf());


    //Efficiency Correction
    efficiencyCorrSpectrum =
      ApplyEfficiencyCorrection(energyLossCorrSpectrum,corrCurve->GetEfficiencyCurve(),minTPCEfficiency,
				corrCurve->GetEfficiencyCurveConf(),corrCurve->GetEfficiencyCurveSys(),
				corrCurve->GetEfficiencyCurveSysConf());

  }//End Systematic Error Mode
    
  //Delete the Intermediate spectra before returning
  delete interBinMeanSpectrum;
  delete muonCorrSpectrum;
  delete feedDownCorrSpectrum;
  delete finiteBinWidthCorrSpectrum;
  delete energyLossCorrSpectrum;


  //Return the fully corrected spectrum which is the last spectrum
  //  return efficiencyCorrSpectrum;
  return efficiencyCorrSpectrum;
}


//_________________________________________________________________
TGraphErrors *ApplyInterBinTransverseMassLocation(TGraphErrors *spectrum, TF1 *interBinAverageCurve,
						  TGraphErrors *interBinAverageCurveSys){

  //Before this function the points in the spectrum are locted at the bin centers of the
  //transverse mass bins. This function uses the parameterizations to the inter-bin averages
  //to set the spectra points at the location of the population mean.

  //If the curve does not exist then return a copy and print a warning
  if (!interBinAverageCurve){

    cout <<"INFO - SpectraCorrectionFunctions::ApplyInterBinTransverseMassLocation() - No function\n"
	 <<"was passed to set the bin locations. This may be because you have not yet done the\n"
	 <<"inter-bin transverse mass study.\n"
	 <<"Correction curves will be evaluated at bin centers.\n" <<endl;
    return (new TGraphErrors(*spectrum));
    
  }

  //Create a new Spectrum
  TGraphErrors *corrSpectrum = new TGraphErrors();

  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){
    
    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint];

    //The interBinAverage Curve is parameterized the mean with respect
    //to the start of the bin. So we need to add the low bin edge
    Double_t lowBinEdge = GetmTm0RangeLow(GetmTm0Index(mTm0));
    Double_t populationMean =  lowBinEdge + interBinAverageCurve->Eval(mTm0);
    Double_t populationMeanStdErr = interBinAverageCurve->GetParError(0);

    //If we are doing systematics, then use the one sigma band
    //to generate a random deviate to use as the population Mean
    if (interBinAverageCurveSys){

      TRandom3 rand(0);
      populationMean = rand.Gaus(populationMean,populationMeanStdErr);
      
    }

    //Set the Point
    corrSpectrum->SetPoint(corrSpectrum->GetN(),populationMean,yield);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,populationMeanStdErr,yErr);
    
  }//End Loop Over Points

  return corrSpectrum;

}

//_________________________________________________________________
std::vector<double> GetCorrectedBinWidthVector(TGraphErrors *spectrum, TF1 *interBinMeanCurve, TF1 *energyLossCurve,
					       Double_t pidMass, TF1 *efficiencyCurve, Double_t minTPCEfficiency,
					       TF1 *feedDownCurve){

  //This function returns a vector containing the corrected bin widths.
  //It makes NO CHANGES to the incomming spectrum itself
  ParticleInfo particleInfo;
  
  //Apply the Inter-Bin Transverse Mass Location
  TGraphErrors *spectrumCopy = ApplyInterBinTransverseMassLocation(spectrum,interBinMeanCurve);

  std::vector<double> correctedBinWidths;

  for (int iPoint=0; iPoint<spectrumCopy->GetN(); iPoint++){

    Double_t mTm0 = spectrumCopy->GetX()[iPoint];
    
    Double_t originalBinLowEdge = ConvertmTm0ToPt(GetmTm0RangeLow(GetmTm0Index(mTm0)),pidMass);   //In pT
    Double_t originalBinHighEdge = ConvertmTm0ToPt(GetmTm0RangeHigh(GetmTm0Index(mTm0)),pidMass); //In pT
    Double_t originalBinWidth = originalBinHighEdge - originalBinLowEdge;                         //In pT

    //Make sure the bin Center is in the range of the energy loss fit
    //if not then don't add a point to the corrected spectrum
    Double_t minFitRange(0), maxFitRange(0);
    energyLossCurve->GetRange(minFitRange,maxFitRange);
    Double_t originalBinCenter = originalBinLowEdge + originalBinWidth/2.0;
    if (originalBinCenter < minFitRange || originalBinCenter > maxFitRange)
      continue;

    //Make sure the TPC tracking efficiency is above the minimum, if not
    //do not add a point to the corrected spectrum
    Double_t pTOfPoint = ConvertmTm0ToPt(mTm0,pidMass);
    Double_t correctionShift = energyLossCurve->Eval(pTOfPoint);
    Double_t correctedpT = pTOfPoint - correctionShift;
    Double_t correctedmTm0 = ComputemTm0(correctedpT,pidMass);
    Double_t effCorrFactor = efficiencyCurve->Eval(correctedmTm0);
    if (effCorrFactor < minTPCEfficiency)
      continue;

    //Make sure the point is in the range of the feed down correction curve
    //Only Do this if the particle is a not a kaon
    if (pidMass != particleInfo.GetParticleMass(KAON)){
      Double_t minFeedDownRange(0), maxFeedDownRange(0);
      if (!feedDownCurve)
	continue;
      feedDownCurve->GetRange(minFeedDownRange,maxFeedDownRange);
      if (mTm0 < minFeedDownRange)
	continue;
    }

    //Compute the Shift
    Double_t lowEdgeShift = energyLossCurve->Eval(originalBinLowEdge);
    Double_t highEdgeShift = energyLossCurve->Eval(originalBinHighEdge);

    Double_t correctedBinLowEdge = originalBinLowEdge - lowEdgeShift;
    Double_t correctedBinHighEdge = originalBinHighEdge - highEdgeShift;
    Double_t correctedBinWidth = correctedBinHighEdge - correctedBinLowEdge; //Corrected Bin Width in pT
    Double_t correctedBinWidthmTm0 = ComputemTm0(correctedBinHighEdge,pidMass) -
      ComputemTm0(correctedBinLowEdge,pidMass); //Corrected Bin Width in mT-m0

    correctedBinWidths.push_back(correctedBinWidthmTm0);
    
  }//End Loop Over Points

  delete spectrumCopy;

  return correctedBinWidths;

}

//_________________________________________________________________
TGraphErrors *ApplyFiniteBinWidthCorrection(TGraphErrors *spectrum, TF1 *energyLossCurve,
					    Double_t pidMass,TGraphErrors *energyLossCurveSys){

  //Loop Over the points of the spectrum and apply the Finite Bin Width (dmT) correction.
  //Return a new spectrum

  if (!energyLossCurve){
    cout <<"INFO - SpectraCorrectionFunctions::ApplyFiniteBinWidthCorrection() - No energy loss curve \n was passed for this spectrum. Returning a copy of the input spectrum.\n";
    return (new TGraphErrors(*spectrum));
  }
  
  TGraphErrors *corrSpectrum = new TGraphErrors();

  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint];

    //Compute the shift for the Lower and Upper Edges of the mTm0 Bin
    //keeping in mind that the energy loss curve is in pT and the
    //spectra is in mT-m0...so we need to convert everything first to pT.
    
    //Original Bin Values
    Double_t originalBinLowEdge = ConvertmTm0ToPt(GetmTm0RangeLow(GetmTm0Index(mTm0)),pidMass);   //In pT
    Double_t originalBinHighEdge = ConvertmTm0ToPt(GetmTm0RangeHigh(GetmTm0Index(mTm0)),pidMass); //In pT
    Double_t originalBinWidth = originalBinHighEdge - originalBinLowEdge;                         //In pT

    //Make sure the bin Center is in the range of the energy loss fit
    //if not then don't add a point to the corrected spectrum
    Double_t minFitRange(0), maxFitRange(0);
    energyLossCurve->GetRange(minFitRange,maxFitRange);
    Double_t originalBinCenter = originalBinLowEdge + originalBinWidth/2.0;
    if (originalBinCenter < minFitRange || originalBinCenter > maxFitRange)
      continue;
    
    //Compute the Shift
    Double_t lowEdgeShift = energyLossCurve->Eval(originalBinLowEdge);
    Double_t highEdgeShift = energyLossCurve->Eval(originalBinHighEdge);

    //The systematics associated with the energy loss correction are computed
    //by varying the value of the curve within its one sigma limits using
    //the errors on the parameters. (See Above)
    if (energyLossCurveSys){

      TRandom3 rand(0);
      lowEdgeShift = rand.Gaus(energyLossCurve->Eval(originalBinLowEdge),
			       GetErrorAtValue(energyLossCurveSys,originalBinLowEdge));
      highEdgeShift = rand.Gaus(energyLossCurve->Eval(originalBinHighEdge),
				GetErrorAtValue(energyLossCurveSys,originalBinHighEdge));
      
    }
    
    //Compute the New Bin Edges by applying the shift
    Double_t correctedBinLowEdge = originalBinLowEdge - lowEdgeShift;
    Double_t correctedBinHighEdge = originalBinHighEdge - highEdgeShift;
    Double_t correctedBinWidth = correctedBinHighEdge - correctedBinLowEdge; //Corrected Bin Width in pT
    Double_t correctedBinWidthmTm0 = ComputemTm0(correctedBinHighEdge,pidMass) -
      ComputemTm0(correctedBinLowEdge,pidMass); //Corrected Bin Width in mT-m0
    
    //Compute the Correction Factor
    Double_t correctionFactor = originalBinWidth/correctedBinWidth; //pT Cancels in ratio
						   
    //Apply the dmT Correction
    corrSpectrum->SetPoint(corrSpectrum->GetN(),mTm0,yield*correctionFactor);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,xErr,yErr*correctionFactor);

  }//End Loop Over Points

  corrSpectrum->SetMarkerStyle(spectrum->GetMarkerStyle());
  corrSpectrum->SetMarkerColor(spectrum->GetMarkerColor());

  return corrSpectrum;
 
}

//_________________________________________________________________
TGraphErrors *ApplyEnergyLossCorrection(TGraphErrors *spectrum, TF1 *energyLossCurve,
					Double_t pidMass, TGraphErrors *energyLossCurveSys){

  //Apply the correction to the measured momentum and return a new spectrum

  if (!energyLossCurve){
    cout <<"INFO - SpectraCorrectionFunctions::ApplyEnergyLossCorrection() - No energy loss curve \n was passed for this spectrum. Returning a copy of the input spectrum.\n";
    return (new TGraphErrors(*spectrum));
  }
    
  TGraphErrors *corrSpectrum = new TGraphErrors();

  //Loop Over the points of the spectrum an apply the correction
  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint];

    //Find the Shift by evaluating the Energy Loss Curve
    //keeping in mind that the energy loss curve is in pT and the
    //spectra is in mT-m0...so we need to convert
    Double_t pTOfPoint = ConvertmTm0ToPt(mTm0,pidMass);

    //Make sure the pT of the point is in the energy Loss fit range
    //if not then do not add a point to the corrected spectrum
    Double_t minFitRange(0), maxFitRange(0);
    energyLossCurve->GetRange(minFitRange,maxFitRange);
    if (pTOfPoint < minFitRange || pTOfPoint > maxFitRange)
      continue;

    //Compute the correction shift. By default  the energy
    //loss curve is evaluated at the pT of the point. Otherwise a
    //random pT in the bin is chosen for systematic error determination.
    Double_t correctionShift = energyLossCurve->Eval(pTOfPoint);
    if (energyLossCurveSys){

      TRandom3 rand(0);
      correctionShift = rand.Gaus(energyLossCurve->Eval(pTOfPoint),
				  GetErrorAtValue(energyLossCurveSys,pTOfPoint));

    }
    
    //Apply the Correction and convert back to mTm0
    Double_t correctedpT = pTOfPoint - correctionShift;
    Double_t correctedmTm0 = ComputemTm0(correctedpT,pidMass);

    //Apply the correction, this correction does not effect the errors
    //so they are propagated through
    corrSpectrum->SetPoint(corrSpectrum->GetN(),correctedmTm0,yield);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,xErr,yErr);

  }

  corrSpectrum->SetMarkerStyle(spectrum->GetMarkerStyle());
  corrSpectrum->SetMarkerColor(spectrum->GetMarkerColor());

  return corrSpectrum;
}

//_________________________________________________________________
TGraphErrors *ApplyEfficiencyCorrection(TGraphErrors *spectrum, TF1 *efficiencyCurve, Double_t minTPCEfficiency,
					TGraphErrors *efficiencyCurveConf, TF1 *efficiencyCurveSys,
					TGraphErrors *efficiencyCurveSysConf){

  //Apply the TPC Tracking Efficiency Correction
  //and return a new spectrum

  if (!efficiencyCurve){
    cout <<"INFO - SpectraCorrectionFunctions::ApplyEfficiencyCorrection() - No efficiency curve \n was passed for this spectrum. Returning a copy of the input spectrum.\n";
    return (new TGraphErrors(*spectrum));
  }
  
  TGraphErrors *corrSpectrum = new TGraphErrors();

  //If we are going to do systematic errors then generate
  //a one sigma deviate function from the systematic curve
  //  TF1 *deviateCurve = NULL;
  //if (efficiencyCurveSys){
  //  deviateCurve = GetOneSigmaEfficiencyDeviateFunction(efficiencyCurveSys);
  //}

  //Loop Over the Points
  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint]; 

    //Get the Efficiency Correction Factor. By default the correction factor
    //is determined by evaluating the efficiency curve. Otherwise, the parameters
    //of the curve are varied (see above) to determine the systematic error.
    Double_t effCorrFactor = efficiencyCurve->Eval(mTm0);

    //Make sure the TPC tracking efficiency is above the minimum, if not
    //do not add a point to the corrected spectrum
    if (effCorrFactor < minTPCEfficiency)
      continue;    

    //If Systematic Errors are requested
    if (efficiencyCurveConf){
      //TRandom3 rand(0);
      //effCorrFactor = rand.Gaus(efficiencyCurve->Eval(mTm0),
      //			GetErrorAtValue(efficiencyCurveConf,mTm0));
      
      TRandom3 rand(0);
      
      //First, Choose which function to use
      int whichFunc = TMath::Nint(rand.Uniform(0,1));

      //Get a new efficiency correction factor
      if (whichFunc == 0)
	effCorrFactor = rand.Gaus(efficiencyCurve->Eval(mTm0),
				  GetErrorAtValue(efficiencyCurveConf,mTm0));
      else
	effCorrFactor = rand.Gaus(efficiencyCurveSys->Eval(mTm0),
				  GetErrorAtValue(efficiencyCurveSysConf,mTm0));
      
    }

    //If the efficiency factor is greater than 1 or less than zero
    //set the correction factor to either 1 or 0 respectively. This imposes
    //the physical boundaries conditions
    if (effCorrFactor > 1.0)
      effCorrFactor = 1.0;
    else if (effCorrFactor <= 0.0)
      effCorrFactor = 0.0001;
    
    //Apply the correction
    corrSpectrum->SetPoint(corrSpectrum->GetN(),mTm0,yield/effCorrFactor);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,xErr,yErr/effCorrFactor);
    
  }//End Loop Over Points

  corrSpectrum->SetMarkerStyle(spectrum->GetMarkerStyle());
  corrSpectrum->SetMarkerColor(spectrum->GetMarkerColor());

  //  if (deviateCurve)
  //delete deviateCurve;
  
  return corrSpectrum;
  
}

//__________________________________________________________________________________
TGraphErrors *ApplyKnockoutProtonCorrection(TGraphErrors *spectrum, TF1 *backgroundCurve, Bool_t doSystematics){

  //Apply the Knockout Proton Correction and return a new spectrum

  if (!backgroundCurve){
    cout <<"INFO - SpectraCorrectionFunctions::ApplyKnockoutCorrection() - No knockout background curve \n was passed for this spectrum. Returning a copy of the input spectrum.\n";
    return (new TGraphErrors(*spectrum));
  }
  
  TGraphErrors *corrSpectrum = new TGraphErrors();

  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint];

    //Find the Background fraction by evaluating the backgroundCurve
    Double_t backgroundFraction = backgroundCurve->Eval(mTm0);

    //If we are doing the systematics then evaluate the curve
    //at a random location within the bin
    if (doSystematics){

      TRandom3 rand(0);
      Double_t randommTm0 = rand.Uniform(mTm0-xErr,mTm0+xErr);
      backgroundFraction = backgroundCurve->Eval(randommTm0);

    }//End Do Systematics

    //Compute the Correction Factor
    Double_t correctionFactor = 1.0 - backgroundFraction;

    //Make sure the correction factor is good
    //if bad then set it to 1 so no change will be made
    if (TMath::IsNaN(correctionFactor))
      correctionFactor = 1.0;

    //Apply the Correction
    Double_t correctedYield = yield * correctionFactor;
    Double_t correctedError = yErr * correctionFactor;

    //Set the Point in the corrected Spectrum
    corrSpectrum->SetPoint(corrSpectrum->GetN(),mTm0,correctedYield);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,xErr,correctedError);

  }//End Loop Over Points

  corrSpectrum->SetMarkerStyle(spectrum->GetMarkerStyle());
  corrSpectrum->SetMarkerColor(spectrum->GetMarkerColor());

  return corrSpectrum;

}

//_________________________________________________________________________________________
TGraphErrors *ApplyMuonBackgroundCorrection(TGraphErrors *spectrum, TF1 *muonBackgroundCurveDefault,
					    TF1 *muonBackgroundCurveAlt,
					    TGraphErrors *muonBackgroundCurveDefaultSys,
					    TGraphErrors *muonBackgroundCurveAltSys){

  //Check that the background curves exists
  if (!muonBackgroundCurveDefault){
    cout <<"INFO -- CorrectSpectra::ApplyMuonBackgroundCorrection() -- NO Muon Background Curve EXISTS! Returning copy of input spectrum. \n";
    return (new TGraphErrors(*spectrum));    
  }
  
  //Apply the Muon Background Correction and return a new spectrum
  TGraphErrors *corrSpectrum = new TGraphErrors();

  //Get the Min and Max Range of the fit function and the spectrum
  Double_t minFitDefault(0), maxFitDefault(0);
  muonBackgroundCurveDefault->GetRange(minFitDefault,maxFitDefault);

  //Loop Over the points of the spectrum
  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint];

    //If the mTm0 Value of the point is less than the minFitRange of correction curve skip
    if (mTm0 < minFitDefault)
      continue;

    //Find the Background Fraction
    Double_t backgroundFraction = muonBackgroundCurveDefault->Eval(mTm0);

    //If we are doing systematics
    if (muonBackgroundCurveAlt){

      //Organize the Curves and One Sigma Bands
      TF1 *curves[2] = {muonBackgroundCurveDefault,muonBackgroundCurveAlt};
      TGraphErrors *sigmaBands[2] = {muonBackgroundCurveDefaultSys,muonBackgroundCurveAltSys};

      TRandom3 rand(0);

      //Choose which curve to use;
      int whichCurve = TMath::Nint(rand.Uniform(0,1));

      //Get the One sigma value from the chosen function at this mT-m0
      double oneSigma = GetErrorAtValue(sigmaBands[whichCurve],mTm0);

      //Set a new background fraction which is a random deviate
      //of a gaussian defined by the value of the curve and
      //its one sigma band.
      backgroundFraction = rand.Gaus(curves[whichCurve]->Eval(mTm0),oneSigma);

    }//End Do Systematics

    //Compute the Correction Factor
    Double_t correctionFactor = 1.0 - backgroundFraction;

    //Apply the Correction Factor
    Double_t correctedYield = yield * correctionFactor;
    Double_t correctedError = yErr * correctionFactor;

    //Set the Point in the Corrected Spectrum
    corrSpectrum->SetPoint(corrSpectrum->GetN(),mTm0,correctedYield);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,xErr,correctedError);
    
  }//End Loop Over Points

  corrSpectrum->SetMarkerStyle(spectrum->GetMarkerStyle());
  corrSpectrum->SetMarkerColor(spectrum->GetMarkerColor());

  return corrSpectrum;

}

//_______________________________________________________________________________
TGraphErrors *ApplyFeedDownBackgroundCorrection(TGraphErrors *spectrum, TF1 *feedDownBackgroundCurveDefault,
						TF1 *feedDownBackgroundCurveAlt,
						TGraphErrors *feedDownBackgroundCurveDefaultSys,
						TGraphErrors *feedDownBackgroundCurveAltSys){
  
  //Check that the background curve exists
  if (!feedDownBackgroundCurveDefault){
    cout <<"INFO -- CorrectSpectra::ApplyFeedDownBackgroundCorrection -- No Feed Down Background Curve EXISTS! Returning a copy of the input spectrum! \n";
    return (new TGraphErrors(*spectrum));
  }


  //Apply the Feed Down Correction and return a new spectrum
  TGraphErrors *corrSpectrum = new TGraphErrors();

  //Get the Min and Max Range of the fit function and the spectrum
  Double_t minFitDefault(0), maxFitDefault(0);
  feedDownBackgroundCurveDefault->GetRange(minFitDefault,maxFitDefault);

  //Loop Over the points of the spectrum
  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){
    
    //Get the point info
    Double_t mTm0  = spectrum->GetX()[iPoint];
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t xErr  = spectrum->GetEX()[iPoint];
    Double_t yErr  = spectrum->GetEY()[iPoint];

    //If the mTm0 Value of the point is less than the minFitRange of the correction curve skip
    if (mTm0 < minFitDefault)
      continue;
    
    //Find the Background Fraction
    Double_t backgroundFraction = feedDownBackgroundCurveDefault->Eval(mTm0);
    
    //If we are doing systematics
    if (feedDownBackgroundCurveAlt){

      //Organize the Curves and one sigma bands
      TF1 *curves[2] = {feedDownBackgroundCurveDefault,feedDownBackgroundCurveAlt};
      TGraphErrors *sigmaBands[2] = {feedDownBackgroundCurveDefaultSys,feedDownBackgroundCurveAltSys};

      TRandom3 rand(0);

      //Choose which curve to use
      int whichCurve = TMath::Nint(rand.Uniform(0,1));

      //Get the One Sigma value from the chosen function at this mT-m0
      double oneSigma = GetErrorAtValue(sigmaBands[whichCurve],mTm0);

      //Set a new feed down fraction which is a random deviate
      //of a gaussian defined by the value of the curve and its
      //one sigma band
      backgroundFraction = rand.Gaus(curves[whichCurve]->Eval(mTm0),oneSigma);
      
      
    }//End Do Systematics
    
    //Compute the Correction Factor
    Double_t correctionFactor = 1.0 - backgroundFraction;

    //Apply the Correction Factor
    Double_t correctedYield = yield * correctionFactor;
    Double_t correctedError = yErr * correctionFactor;

    //Set the Point in the Corrected Spectrum
    corrSpectrum->SetPoint(corrSpectrum->GetN(),mTm0,correctedYield);
    corrSpectrum->SetPointError(corrSpectrum->GetN()-1,xErr,correctedError);
    

  }//End Loop Over Points

  corrSpectrum->SetMarkerStyle(spectrum->GetMarkerStyle());
  corrSpectrum->SetMarkerColor(spectrum->GetMarkerColor());

  return corrSpectrum;

}

//____________________________________________________________________
//TGraphErrors *MakeSystematicErrorSpectrum(TGraphErrors *nominalSpectrum, TGraphErrors *widthSpectrum,
//					  const int nSpectra, std::vector<TGraphErrors *> *sysSpectra){

TGraphErrors *MakeSystematicErrorSpectrum(TGraphErrors *nominalSpectrum,
					  const int nSpectra, std::vector<TGraphErrors *> *sysSpectra){
  
  //Return a graph whose errors are the systematic errors
  TGraphErrors *sysSpectrum = new TGraphErrors();
  sysSpectrum->SetName(Form("%s_SysErr",nominalSpectrum->GetName()));
  sysSpectrum->SetMarkerSize(0);
  sysSpectrum->SetFillStyle(3002);
  sysSpectrum->SetFillColor(nominalSpectrum->GetMarkerColor()-3);

  //For each point in the nominal spectrum compute the systematic errors from the
  //ensemeble of spectra. Use the central values of the points in the nominal spectrum
  //as the central values in the systematic spectrum. Compute the systematic error
  //and the mean of the difference between the nominal central value and the distribution
  //of the ensemble points.
  for (int iPoint=0; iPoint<nominalSpectrum->GetN(); iPoint++){
    //Get the Nominal Point Value
    Double_t mTm0 = nominalSpectrum->GetX()[iPoint];
    Double_t yield = nominalSpectrum->GetY()[iPoint];
  //Double_t mTm0_width = widthSpectrum->GetX()[iPoint];
  // Double_t yield_width = widthSpectrum->GetY()[iPoint];

    std::vector<double> ySysVec(0);
    std::vector<double> xSysVec(0);
    for (int iSpectra=0; iSpectra<nSpectra; iSpectra++){
      if (sysSpectra->at(iSpectra)->GetN() == 0)
	continue;
      int jPoint = TGraphFindPoint(sysSpectra->at(iSpectra),mTm0);
      if (fabs(mTm0 - sysSpectra->at(iSpectra)->GetX()[jPoint]) > mTm0BinWidth/2.0)
	continue;
      ySysVec.push_back(fabs(yield - sysSpectra->at(iSpectra)->GetY()[iPoint]));
      xSysVec.push_back(fabs(mTm0  - sysSpectra->at(iSpectra)->GetX()[iPoint]));
    }//End Loop Over Ensemble
    
    //    ySysVec.push_back(fabs(yield - yield_width));
    //    xSysVec.push_back(fabs(mTm0  - mTm0_width));



    if (ySysVec.size() == 0)
      continue;

    double xErrSys = TMath::Mean(xSysVec.size(),&xSysVec.at(0));
    double yErrSys = TMath::Mean(ySysVec.size(),&ySysVec.at(0));

    sysSpectrum->SetPoint(sysSpectrum->GetN(),mTm0,yield);
    sysSpectrum->SetPointError(sysSpectrum->GetN()-1,xErrSys,yErrSys);
    
  }//End Loop Over Points

  return sysSpectrum;

}

//__________________________________________________________________
TGraphErrors *MakeTotalErrorSpectrum(TGraphErrors *nominalSpectrum, TGraphErrors *systematicSpectrum){

  //Return a spectrum whose errors are the quadrature sum of the
  //statistical error in the nominal spectrum and the systematic error
  //in the systematic spectrum
  TGraphErrors *totalErrorSpectrum = new TGraphErrors();
  totalErrorSpectrum->SetName(Form("%s_TotalErr",nominalSpectrum->GetName()));
  totalErrorSpectrum->SetMarkerStyle(nominalSpectrum->GetMarkerStyle());
  totalErrorSpectrum->SetMarkerColor(nominalSpectrum->GetMarkerColor());

  if (systematicSpectrum->GetN() == 0){
    cout << "returing total error spectrum" << endl;
    return totalErrorSpectrum;
  }

  for (int iPoint=0; iPoint<nominalSpectrum->GetN(); iPoint++){

    double mTm0 = nominalSpectrum->GetX()[iPoint];
    int jPoint = TGraphFindPoint(systematicSpectrum,mTm0);
    if (fabs(mTm0 - systematicSpectrum->GetX()[jPoint]) > mTm0BinWidth/2.0)
      continue;

    double totalErrX = TMath::Sqrt(TMath::Power(nominalSpectrum->GetEX()[iPoint],2) +
				   TMath::Power(systematicSpectrum->GetEX()[jPoint],2));
    double totalErrY = TMath::Sqrt(TMath::Power(nominalSpectrum->GetEY()[iPoint],2) +
				   TMath::Power(systematicSpectrum->GetEY()[jPoint],2));

    totalErrorSpectrum->SetPoint(totalErrorSpectrum->GetN(),
				 nominalSpectrum->GetX()[iPoint],
				 nominalSpectrum->GetY()[iPoint]);
    totalErrorSpectrum->SetPointError(totalErrorSpectrum->GetN()-1,
				      totalErrX,
				      totalErrY);

    
  }//End Loop Over Points
  
  return totalErrorSpectrum;

}


//__________________________________________________________________
TF1 *GetOneSigmaEfficiencyDeviateFunction(TF1 *func){

  //Return a function of the same functional form as func, but whose
  //parameters have been set by sampling gaussians for each variable.

  if (!func){
    cout <<"WARNING - SpectraCorrectionFunctions::GetOneSigmaDevaiteFunction - The input function does not exist! \n Returning NULL\n";
    return NULL;
  }
  
  Double_t minRange(0), maxRange(0);
  func->GetRange(minRange,maxRange);
  
  TF1 *deviateFunc = new TF1("temp",Form("%s",func->GetExpFormula().Data()),
			     minRange,maxRange);

  Bool_t goodFunc=false;
  TRandom3 rand(0);

  //Generate a random function and make sure it is good. There is some
  //posibility that some values of the paramters for the efficiency function
  //will make it decrease as a function of transverse mass rather than increase.
  //This this is non-physical we reject those functions
  while (!goodFunc){
    for (int iPar=0; iPar<func->GetNpar(); iPar++){
      
      deviateFunc->FixParameter(iPar,rand.Gaus(func->GetParameter(iPar),
					       func->GetParError(iPar)));    
    }//End Loop Over Parameters

    if (deviateFunc->Eval(minRange) < deviateFunc->Eval(maxRange) &&
	deviateFunc->Eval(minRange) >= 0 &&
	deviateFunc->Eval(maxRange) <= 1)
      goodFunc = true;

  }//End Make Sure the Function is good
  
  return deviateFunc;
  
}

//___________________________________________________________________
std::vector<int> CleanSpectrum(TGraphErrors *gr, int pid){

  //Remove points in the spectrum that have a relative error
  //greater than 3sigma of the mean relative error.

  //Remove points in the spectrum that are different from their
  //neightbor by more than 3sigma of the mean difference

  //Return vector of point indices that have been removed
  std::vector<int> removedPoints(0);

  if (gr->GetN() < 1)
    return removedPoints;
  
  std::vector<double> relativeError(gr->GetN(),0);
  std::vector<double> percentDiff(gr->GetN(),0);

  for (int iPoint=0; iPoint<gr->GetN(); iPoint++){

    relativeError.at(iPoint) = gr->GetEY()[iPoint] / gr->GetY()[iPoint];

    int neighbor(iPoint+1);
    if (iPoint == gr->GetN()-1)
      neighbor = iPoint-1;
    percentDiff.at(iPoint) = fabs(gr->GetY()[iPoint] - gr->GetY()[neighbor]);
    percentDiff.at(iPoint) /= (gr->GetY()[iPoint] + gr->GetY()[neighbor])/2.0;
		   
  }//End Loop Over Points

  //Compute the Means and RMS Values
  double meanRelativeError = TMath::Mean(relativeError.size(),&relativeError.at(0));
  double rmsRelativeError = TMath::RMS(relativeError.size(),&relativeError.at(0));

  double meanPercentDiff = TMath::Mean(percentDiff.size(),&percentDiff.at(0));
  double rmsPercentDiff  = TMath::RMS(percentDiff.size(),&percentDiff.at(0));

  //Remove Points that have values greater than nSigma Sigma compared to the mean
  //or that have errors which extend to below zero
  double nSigma(5);
  for (int iPoint=0; iPoint<gr->GetN(); iPoint++){

    //For Kaons, only go up to mTm0=1GeV
    if (pid == KAON){
      if (gr->GetX()[iPoint] > 1.0){
	removedPoints.push_back(iPoint);
	cout <<iPoint <<" This is a Kaon Spectrum - This point is removed because it is over 1GeV." <<endl;
	continue;
      }
    }

    if (TMath::IsNaN(gr->GetEY()[iPoint])){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed due to y error being NAN" <<endl;
      continue;
    }

    if (TMath::IsNaN(gr->GetEX()[iPoint])){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed due to x error being NAN" <<endl;
      continue;
    }
    
    if (gr->GetY()[iPoint] < pow(10,-12)){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed due to being zero." <<endl;
      continue;
    }

    if (gr->GetY()[iPoint] - gr->GetEY()[iPoint] < 0){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed since Y error extended below zero" <<endl;
      continue;
    }

    if (gr->GetX()[iPoint] - gr->GetEX()[iPoint] < 0){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed since X error extended below zero" <<endl;
      continue;
    }
    
    if (fabs(relativeError.at(iPoint) - meanRelativeError) > nSigma*rmsRelativeError){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed due to error." <<endl;
      continue;
    }

    if (fabs(percentDiff.at(iPoint) - meanPercentDiff) > nSigma*rmsPercentDiff){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed due to percentdiff being to big" <<endl;
      continue;
    }

    if (gr->GetEY()[iPoint]/gr->GetY()[iPoint] < .00001){
      removedPoints.push_back(iPoint);
      cout <<iPoint <<" Removed due to relative error being too small." <<endl;
      continue;
    }

  }//End Loop Over Points

  if (removedPoints.size() == 0)
    return removedPoints;
  
  for (unsigned int i=removedPoints.size()-1; (int)i>=0; i--){
    gr->RemovePoint(removedPoints.at(i));
    cout <<"RemovedPoint: " <<removedPoints.at(i) <<endl;
  }
  
  return removedPoints;  
 
}
