#include <iostream>
#include <utility>
#include <vector>

#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom3.h>

#include "utilityFunctions.h"
#include "ParticleInfo.h"
#include "SpectraClass.h"
#include "SpectraFitUtilities.h"

/*******************************************************************
Compute the Total dNdy - the total dNdy is computed as the sum of
three dNdy values: 
dNdy_Total = dNdy_LowExtrap + dNdy_Fiducial + dNdy_HighExtrap
 *******************************************************************/
std::pair<double,double> ComputeTotaldNdy(SpectraClass *spectrum, TF1 *lowmTFit, TF1 *highmTFit){

  //Notes:
  //The dNdy Parameter of the fits is assumed to be parameter 0 - which is my convention.
  //The pairs below contain <dNdy,dNdy_Error>
  
  std::pair<double,double> dNdy_LowExtrap = std::make_pair(0,0);
  std::pair<double,double> dNdy_Fiducial;
  std::pair<double,double> dNdy_HighExtrap = std::make_pair(0,0);

  //------Low mT Extrapolation
  if (lowmTFit){
    dNdy_LowExtrap.first = lowmTFit->GetParameter(0);
    dNdy_LowExtrap.second = lowmTFit->GetParError(0);
  }
  else {
    cout <<"WARNING - SpectraFitUtilities::ComputeTotaldNdy() - the lowmTFit does not exist. dNdy will not contain yield extrapolated to mT=0.\n";
  }

  //------Fiducial
  dNdy_Fiducial = spectrum->ComputeFiducialdNdy();

  //------High mT Extrapolation
  if (highmTFit){
    dNdy_HighExtrap.first = highmTFit->GetParameter(0);
    dNdy_HighExtrap.second = highmTFit->GetParError(0);
  }
  else {
    cout <<"WARNING - SpectraFitUtilties::ComputeTotaldNdy() = the highmTFit does not exist. dNdy will not contain yield extrapolated to mT values higher than the last point of the spectrum.\n";
  }

  double totalYield = dNdy_LowExtrap.first + dNdy_Fiducial.first + dNdy_HighExtrap.first;
  double totalYieldErr = sqrt( pow(dNdy_LowExtrap.second,2) +
			       pow(dNdy_Fiducial.second,2) +
			       pow(dNdy_HighExtrap.second,2) );

  return std::make_pair(totalYield,totalYieldErr);
  
}



/******************************************************************
Vary the fit length to until a good Chi^2/NDF is found
 ******************************************************************/
//____Vary the Upper range_____________________________________________________________
double VaryFitLengthUpper(SpectraClass *spectrum, TF1 *fitFunc, TGraphErrors *chi2Graph){

  //Vary the Upper range of the fit over the length of the spectrum.
  //Return the length which gives the chi2/NDF closest to 1

  std::vector<double> chi2Vals;
  std::vector<double> maxRangeVals;
  Double_t minRange = spectrum->GetX()[0] - spectrum->GetEX()[0];
  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    Double_t maxRange = spectrum->GetX()[iPoint] + spectrum->GetBinWidthVector()->at(iPoint)/2.0;

    spectrum->Fit(fitFunc,"NQ","",minRange,maxRange);
    
    //Require NDF to be at least 3
    if (fitFunc->GetNDF() < 8)
      continue;
    
    chi2Vals.push_back(fabs( (fitFunc->GetChisquare()/(double)fitFunc->GetNDF()) - 1.0) );
    maxRangeVals.push_back(maxRange);
    
    if (chi2Graph){
      chi2Graph->SetPoint(chi2Graph->GetN(),maxRange,chi2Vals.back());
    }
    
  }//End Loop Over Points of spectrum

  int minIndex = TMath::LocMin(chi2Vals.size(),&chi2Vals.at(0));

  return maxRangeVals.at(minIndex);
}

//____Vary The Lower Range________________________________________________________________
double VaryFitLengthLower(SpectraClass *spectrum, TF1 *fitFunc, TGraphErrors *chi2Graph){


  //Vary the Lower range of the fit over the length of the spectrum.
  //Return the length that gives the chi2/NDF closest 1.

  std::vector<double> chi2Vals;
  std::vector<double> minRangeVals;
  Double_t maxRange = spectrum->GetX()[spectrum->GetN()-1] + spectrum->GetEX()[spectrum->GetN()-1];
  for (int iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    Double_t minRange = spectrum->GetX()[iPoint] - spectrum->GetBinWidthVector()->at(iPoint)/2.0;

    spectrum->Fit(fitFunc,"NQ","",minRange,maxRange);

    if (fitFunc->GetNDF() < 8)
      continue;

    chi2Vals.push_back(fabs( (fitFunc->GetChisquare()/(double)fitFunc->GetNDF()) -1.0 ) );
    minRangeVals.push_back(minRange);

    if (chi2Graph)
      chi2Graph->SetPoint(chi2Graph->GetN(),minRange,chi2Vals.back());
    
  }//End Loop Over Points of the Spectrum

  int minIndex = TMath::LocMin(chi2Vals.size(),&chi2Vals.at(0));

  return minRangeVals.at(minIndex);

}


/*********************************************************************************
Get Errors associated with fixing slopt parameters
**********************************************************************************/
Double_t GetPercentErrorFromFixedSlope(SpectraClass *spectrum, TF1 *fitFunc,
				       TGraphErrors *slopeConfInterval, Double_t rapidity,
				       Double_t minFitRange, Double_t maxFitRange){

  //Vary the fixed slope using the confidence interval to determine
  //the error associated with fixing the slope

  //Create a Copy of the Fit Function
  TF1 *tempFitFunc = new TF1(*fitFunc);

  Double_t defaultYield = fitFunc->GetParameter(0);

  tempFitFunc->FixParameter(1,slopeConfInterval->Eval(rapidity));
  spectrum->Fit(tempFitFunc,"NQ","",minFitRange,maxFitRange);

  Double_t percentDiff = fabs( tempFitFunc->GetParameter(0) - defaultYield ) / defaultYield;

  return percentDiff;

  /*
  const int nTrials(25);
  std::vector<double> yieldFractions(nTrials);
  TRandom3 rand(0);
  for (unsigned int iTrial=0; iTrial<nTrials; iTrial++){

    tempFitFunc->FixParameter(1,rand.Gaus(slopeConfInterval->Eval(rapidity),
					  GetErrorAtValue(slopeConfInterval,rapidity)));

    spectrum->Fit(tempFitFunc,"NQ","",minFitRange,maxFitRange);

    yieldFractions.at(iTrial) = fabs( tempFitFunc->GetParameter(0) - defaultYield ) / defaultYield;
    
  }//End Loop Over Trials
  
  return TMath::MaxElement(nTrials,&yieldFractions.at(0));
  */
}
