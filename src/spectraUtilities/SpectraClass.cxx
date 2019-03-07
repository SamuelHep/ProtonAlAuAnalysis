/********************************************
             !!!!WARNING!!!!
This class is constructed assumming that the
spectra it holds are ploted as a function
of mT-m0 and with normalization factors
(1/2pi)x(1/mT)x(d^2N/dmTdy). All other types
of spectra will result in erroneous values!
*********************************************/

#include <iostream>
#include <utility>

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TMatrixDSym.h>

#include "utilityFunctions.h"
#include "ParticleInfo.h"
#include "SpectraClass.h"

ClassImp(SpectraClass);

//_________________________________________________________________________________
SpectraClass::SpectraClass(){

  pid = -999;
  charge = -999;
  centIndex = -999;
  yIndex = -999;
  nPoints = -999;
  name = "";
  nSpectraFits = 0;
  
  dNdyResult=-999;
  dNdyResultErr = -999;
  dNdyResultStatErr = -999;
  dNdyResultSysErrHigh = -999;
  dNdyResultSysErrLow = -999;
  meanpTResult = -999;
  meanpTResultErr = -999;
  meanmTResult = -999;
  meanmTResultErr = -999;
 
}

//_________________________________________________________________________________
SpectraClass::SpectraClass(int pidIndex, int pidCharge, int centBin, int yBin,
			   TGraphErrors *tpc, TGraphErrors *tof){

  tpcSpectrum = NULL;
  tofSpectrum = NULL;
  
  if (tpc){
    if (tpc->GetN() > 0)
      tpcSpectrum = new TGraphErrors(*tpc);
  }


  if (tof){
    if (tof->GetN() > 0)
      tofSpectrum = new TGraphErrors(*tof);
  }

  
  pid = pidIndex;
  charge = pidCharge;
  centIndex = centBin;
  yIndex = yBin;
  nSpectraFits = 0;
  
  name = TString::Format("combinedSpectra_%s_Cent%02d_yIndex%02d",
			 particleInfo.GetParticleName(pid,charge).Data(),
			 centIndex,yIndex);

  dNdyResult=-999;
  dNdyResultErr = -999;
  dNdyResultStatErr = -999;
  dNdyResultSysErrHigh = -999;
  dNdyResultSysErrLow  = -999;
  meanpTResult = -999;
  meanpTResultErr = -999;
  meanmTResult = -999;
  meanmTResultErr = -999;
    
  MakeMultiGraph();
  if (binWidths.size() > 0)
    ComputeFiducialdNdy();
  
}

//_________________________________________________________________________________
SpectraClass::~SpectraClass(){

  //The tpc and tof spectra are deleted internally
  //by the call to the multigraph destructor
  
  if (totalSpectrum)
    delete totalSpectrum;
  if (xVals)
    delete xVals;
  if (yVals)
    delete yVals;

  for (unsigned int i=0; i<spectraFits.size(); i++){
    if (spectraFits.at(i))
      delete spectraFits.at(i);
    if (spectraFitCovarianceMatrix.at(i))
      delete spectraFitCovarianceMatrix.at(i);
  }
  spectraFits.clear();
  spectraFitCovarianceMatrix.clear();
  
}

//_________________________________________________________________________________
int SpectraClass::LoadSpectra(){

  if (!spectraFile){
    cout <<"ERROR - SpectraClass::LoadSpectra() - SpectraFile pointer not set!. Please use SpectraClass::SetSpectraFile() and SetSpectraDir() to set the spectra file and spectra directory. EXITING!\n";
    exit (EXIT_FAILURE);
  }

  if (!spectraDir.CompareTo("")){
    cout <<"ERROR - SpectraClass::LoadSpectra() - SpectraDirectory string is not set. Please use SpectraClass::SetSpectraDir() to set the location of the spectrum directory. EXITING\n";
    exit (EXIT_FAILURE);
  }

  TString baseName = "correctedSpectra";
  if (namePrefix.CompareTo(""))
    baseName = namePrefix;
  
  spectraFile->cd();
  tpcSpectrum = (TGraphErrors *)spectraFile->Get(Form("%s/%s_%s_Cent%02d_yIndex%02d",
						      spectraDir.Data(),
						      baseName.Data(),
						      particleInfo.GetParticleName(pid,charge).Data(),
						      centIndex,yIndex));
  
  tofSpectrum = (TGraphErrors *)spectraFile->Get(Form("%s/%sTOF_%s_Cent%02d_yIndex%02d",
						      spectraDir.Data(),
						      baseName.Data(),
						      particleInfo.GetParticleName(pid,charge).Data(),
						      centIndex,yIndex));

  if (!tpcSpectrum && ! tofSpectrum){
    cout <<"WARNING: SpectraClass::LoadSpectra() - No Spectra Loaded!" <<endl;
    return 0 ;
  }
  
  if (tpcSpectrum){
    tpcSpectrum->SetMarkerStyle(kFullCircle);
  }
  if (tofSpectrum){
    tofSpectrum->SetMarkerStyle(kFullSquare); 
  }

  MakeMultiGraph();
  if (dNdyFiducial.first < .0001 && binWidths.size() > 0)
    ComputeFiducialdNdy();

  //Load the Spectra Fit Functions here
  for (unsigned int i=0; i<nSpectraFits; i++){
    //    cout <<Form("%s/%s_Fit%d",
    //		spectraFitDir.Data(),
    //		name.Data(),i+1) <<endl;
    spectraFits.push_back((TF1 *)spectraFile->Get(Form("%s/%s_Fit%d",
						       spectraFitDir.Data(),
						       name.Data(),i)));
  }

  return 1;
}

//_________________________________________________________________________________
void SpectraClass::MakeMultiGraph(){

  int nPointsTPC(0);
  int nPointsTOF(0);

  //If neither Spectra Exist, and the Dir variable does exist
  //then load them from the directory
  if (!tpcSpectrum && !tofSpectrum){
    LoadSpectra();
  }
  
  totalSpectrum = new TMultiGraph();
  if (tpcSpectrum){
    totalSpectrum->Add(tpcSpectrum);
    nPointsTPC = tpcSpectrum->GetN();
  }
  if (tofSpectrum){
    totalSpectrum->Add(tofSpectrum);
    nPointsTOF = tofSpectrum->GetN();
  }

  nPoints = nPointsTPC + nPointsTOF;
  xVals = new double[nPoints];
  yVals = new double[nPoints];
  exVals = new double[nPoints];
  eyVals = new double[nPoints];
  for (int i=0; i<nPointsTPC; i++){
    xVals[i]  = tpcSpectrum->GetX()[i];
    exVals[i] = tpcSpectrum->GetEX()[i];
    yVals[i]  = tpcSpectrum->GetY()[i];
    eyVals[i] = tpcSpectrum->GetEY()[i];
    
  }
  for (int i=0; i<nPointsTOF; i++){
    xVals[nPointsTPC+i] = tofSpectrum->GetX()[i];
    exVals[nPointsTPC+i] = tofSpectrum->GetEX()[i];
    yVals[nPointsTPC+i] = tofSpectrum->GetY()[i];
    eyVals[nPointsTPC+i] = tofSpectrum->GetEY()[i];
  }

  totalSpectrum->SetTitle(Form(";m_{T}-m_{%s} (GeV/c^{2});%s",
			       particleInfo.GetParticleSymbol(pid,charge).Data(),
			       GetSpectraAxisTitle().Data()));
  
}

//________________________________________________________________________________
void SpectraClass::Fit(TF1 *fitFunc, TString opts, TString gopts, Axis_t minRange, Axis_t maxRange){

  totalSpectrum->Fit(fitFunc,Form("%s",opts.Data()),gopts.Data(),minRange,maxRange);
  
}

//_______________________________________________________________________________
void SpectraClass::WriteSpectrum(bool overwriteSpectra, bool overwriteSpectraClass){

  //NOTE: Objects are overwritten if overwrite values are set to true.
  //      If overwrites are set to false, then the object is only written
  //      if an object of the same name does not already exist.

  if (!spectraFile){
    cout <<"WARNING - SpectraClass::WriteSpectrum() - The SpectraFile is not defined. Please call SpectraClass::SetSpectraFile() before trying to write the spectrum. NOTHING WRITTEN!\n";
    return;
  }
  
  //If the necessary Directories do not exist, then make them
  spectraFile->cd();
  if (!gDirectory->GetKey(spectraClassDir.Data()))
    gDirectory->mkdir(spectraClassDir.Data());
  if (!gDirectory->GetKey(spectraDir.Data()) && (tpcSpectrum || tofSpectrum) )
    gDirectory->mkdir(spectraDir.Data());
  if (!gDirectory->GetKey(spectraDir.Data()) && (spectraFits.size() > 0) )
    gDirectory->mkdir(spectraFitDir.Data());
  
  //Write the Spectra Class Object
  spectraFile->cd();
  spectraFile->cd(spectraClassDir.Data());
  if (overwriteSpectraClass){
    ((TObject *)this)->Write(Form("%s",name.Data()),
  			     TObject::kOverwrite);
  }
  else {
    if (!gDirectory->GetKey(name.Data()))
      ((TObject *)this)->Write(Form("%s",name.Data()));
  }//End Write spectra class

  //Write the Individual Spectras
  if (tpcSpectrum || tofSpectrum){
    spectraFile->cd();
    spectraFile->cd(spectraDir.Data());
    if (overwriteSpectra){
      if (tpcSpectrum)
	tpcSpectrum->Write(Form("%s",tpcSpectrum->GetName()),TObject::kOverwrite);
      if (tofSpectrum)
	tofSpectrum->Write(Form("%s",tofSpectrum->GetName()),TObject::kOverwrite);
    }
    else {
      if (tpcSpectrum)
	if (!gDirectory->GetKey(tpcSpectrum->GetName()))
	  tpcSpectrum->Write();
      if (tofSpectrum)
	if (!gDirectory->GetKey(tofSpectrum->GetName()))
	  tofSpectrum->Write();
    }
  }//End Write individual spectra
  
  //Write the Spectra Fits
  if (spectraFits.size() > 0){
    spectraFile->cd();
    spectraFile->cd(spectraFitDir.Data());
    for (unsigned int i=0; i<spectraFits.size(); i++){
      if (overwriteSpectraClass){
	if (spectraFits.at(i))
	  spectraFits.at(i)->Write(Form("%s",spectraFits.at(i)->GetName()),TObject::kOverwrite);
      }
      else{
	if (spectraFits.at(i))
	  if (!gDirectory->GetKey(spectraFits.at(i)->GetName()))
	    spectraFits.at(i)->Write();
      }
    }
  }//End Write Spectra Fits
 
}

//______________________________________________________________________________
void SpectraClass::DrawSpectrum(TString opts){

  totalSpectrum->Draw(Form("%s",opts.Data()));

}


//____________________________________________________________
Double_t SpectraClass::GetdNdyFiducial(){

  //If the Fiducial dNdy has already been set then return it
  if (dNdyFiducial.first > 0.00001)
    return dNdyFiducial.first;

  //Otherwise calculate it
  ComputeFiducialdNdy();

  return GetdNdyFiducial();
    
}

//____________________________________________________________
Double_t SpectraClass::GetdNdyFiducialErr(){

  //If the Fiducial dNdy Err has already been set, then return it
  if (dNdyFiducial.second > 0.00000001)
    return dNdyFiducial.second;

  //Otherwise compute it
  ComputeFiducialdNdy();

  return GetdNdyFiducialErr();
  
}

//____________________________________________________________
double SpectraClass::GetdNdy(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetdNdy() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return dNdyVec.at(which).first;
  
}

//____________________________________________________________
double SpectraClass::GetdNdyErr(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetdNdy() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return dNdyVec.at(which).second;
  
}

//____________________________________________________________
double SpectraClass::GetMeanmT(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetMeanmT() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  if (meanmTVec.at(which).first > 0.000001)
    return meanmTVec.at(which).first;

  else
    ComputeMeanmTm0(which);

  return GetMeanmT(which);
  
}

//____________________________________________________________
double SpectraClass::GetMeanmTErr(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetMeanmTErr() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }
  
  if (meanmTVec.at(which).first > 0.00000001)
    return meanmTVec.at(which).second;

  else
    ComputeMeanmTm0(which);

  return GetMeanmTErr(which);
  
}

//____________________________________________________________
double SpectraClass::GetMeanpT(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetMeanpT() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return ConvertmTm0ToPt(GetMeanmT(which),particleInfo.GetParticleMass(pid));
  
}

//____________________________________________________________
double SpectraClass::GetMeanpTErr(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetMeanpTErr() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }
  
  return GetMeanmTErr(which);
  
}

//___________________________________________________________
TF1 *SpectraClass::GetSpectrumFit(unsigned int which){

  if (which > spectraFits.size()){
    cout <<"WARNING - SpectraClass::GetSpectrumFit - The index of the requested spectrum fit does not exist!\n Returning NULL.\n";
    return NULL;
  }

  return spectraFits.at(which);  
}

//___________________________________________________________
double SpectraClass::GetdNdyLowExtrapolation(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetdNdyLowExtrapolation() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return dNdyLowExtrap.at(which).first;
  
}

//___________________________________________________________
double SpectraClass::GetdNdyLowExtrapolationErr(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetdNdyLowExtrapolationErr() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return dNdyLowExtrap.at(which).second;
}

//___________________________________________________________
double SpectraClass::GetdNdyHighExtrapolation(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetdNdyHighExtrapolation() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return dNdyHighExtrap.at(which).first;
  
}

//___________________________________________________________
double SpectraClass::GetdNdyHighExtrapolationErr(unsigned int which){

  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::GetdNdyHighExtrapolationErr() - Fit index does not exist. Returning 0.\n";
    return 0.0;
  }

  return dNdyHighExtrap.at(which).second;
}

//___________________________________________________________
void SpectraClass::SetdNdy(double val, unsigned int which){
  
  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::SetdNdy() - Fit index does not exist. Nothing done.\n";
  }

  dNdyVec.at(which).first = val;
}

//___________________________________________________________
void SpectraClass::SetdNdyStatErr(double val, unsigned int which){
  
  if (which >= nSpectraFits){
    cout <<"WARNING - SpectraClass::SetdNdyStatErr() - Fit index does not exist. Nothing done.\n";
  }

  dNdyVec.at(which).second = val;
}



//___________________________________________________________
void SpectraClass::AddSpectraFit(TF1 *fit, TMatrixDSym *matrix, int which){
  
  //  cout << "**DEBUG A 1" << endl;

  //Add the function fit to the list of spectra fits
  if (which < 0){
    spectraFits.push_back(new TF1(*fit));
    which = spectraFits.size() -1;
  }
  else{
    if (spectraFits.size() > (unsigned int)which){
      if (spectraFits.at(which))
	delete spectraFits.at(which);
      spectraFits.at(which) = new TF1(*fit);
    }      
    else
      spectraFits.insert(spectraFits.begin()+which,new TF1(*fit));
  }

  //  cout << "**DEBUG A 1.1" << endl;
  //Set the Chi2 and NDF
  spectraFits.at(which)->SetChisquare(fit->GetChisquare());
  spectraFits.at(which)->SetNDF(fit->GetNDF());

  //  cout << "**DEBUG A 1.2" << endl;

  //Resize all of the internal vectors
  spectraFitCovarianceMatrix.resize(spectraFits.size(),(TMatrixDSym *)NULL);
  //    cout << "**DEBUG A 1.3" << endl;
  dNdyLowExtrap.resize(spectraFits.size(),std::make_pair(0,0));
  dNdyHighExtrap.resize(spectraFits.size(),std::make_pair(0,0));
  //  cout << "**DEBUG A 1.4" << endl;
  dNdyVec.resize(spectraFits.size(),std::make_pair(0,0));
  meanpTVec.resize(spectraFits.size(),std::make_pair(0,0));
  meanmTVec.resize(spectraFits.size(),std::make_pair(0,0));

  //  cout << "**DEBUG A 2" << endl;

  //Reset the number of spectra fits
  nSpectraFits = (int)spectraFits.size();

  //  cout << "**DEBUG A 3" << endl;
  //Add the Covariance Matrix of the fit
  AddSpectraFitCovMatrix(matrix,which);

  //  cout << "**DEBUG A 4" << endl;
  //Compute the Quantities Associated with the fit
  ComputeInternalVectorQuantities(which);

  //  cout << "**DEBUG A 5" << endl;



}

 
//______________________________________________________________
void SpectraClass::AddSpectraFitCovMatrix(TMatrixDSym *matrix, int which){

  //Protected Again NULL Covariance Matrix
  if (!matrix){
    return;
  }

  //  cout << "**DEBUG B 1" << endl;

  if (which < 0)
    spectraFitCovarianceMatrix.push_back(new TMatrixDSym(*matrix));
  

  else 
  
    spectraFitCovarianceMatrix.at(which) = new TMatrixDSym(*matrix);
  
  //  cout << "**DEBUG B 2" << endl;

  //Recompute the internal quantities becuase some of them depend on
  //the covariance Matrix
  ComputeInternalVectorQuantities(which);

  //  cout << "**DEBUG B 3" << endl;
}

//_______________________________________________________________
void SpectraClass::PrintSummary(int which, bool printPoints){

  if (printPoints){

    for (int iPoint=0; iPoint<nPoints; iPoint++){
      cout <<Form("%.04g\t%.04g\t%.04g\t%.04g\n",
		  xVals[iPoint], exVals[iPoint],
		  yVals[iPoint], eyVals[iPoint]);
    }
    
  }
  
  else {

    Double_t rapidity = GetRapidityRangeCenter(yIndex);
    if (fabs(rapidity) < .01)
      rapidity = 0.0;
    
    cout <<Form("%.04g\t%.06g\t%.06g\t%.04g\t%.04g\t%.04g\t%.04g\n",
		rapidity,
		dNdyResult,
		dNdyResultErr,
		meanmTResult,
		meanmTResultErr,
		meanpTResult,
		meanpTResultErr);
  }
}



/**************************************************************************
FUNCTIONS THAT COMPUTE QUANTITIES
 **************************************************************************/
//____________________________________________________________
void SpectraClass::ComputeInternalVectorQuantities(int which){

  //Compute and set the quantities using the public functions
  //  cout << "***DEBUG C 1 " << endl;
  meanmTVec.at(which) = ComputeMeanmTm0(which);
  //  cout << "***DEBUG C 2 " << endl;
  meanpTVec.at(which) = std::make_pair(GetMeanpT(which),GetMeanpTErr(which));
  //  cout << "***DEBUG C 3 " << endl;
  dNdyVec.at(which) = ComputeTotaldNdy(which);
  //  cout << "***DEBUG C 4 " << endl;
  //Set the Internal quantities with the private functions
  ComputedNdyHighExtrapolation(which);    
  //  cout << "***DEBUG C 5 " << endl;
  ComputedNdyLowExtrapolation(which);
  //  cout << "***DEBUG C 6 " << endl;
}

//_______________________________________________________________
void SpectraClass::ComputedNdyHighExtrapolation(int which){

  //Compute the amount of dNdy that results from the extrapolation
  //of the fit function to transverse mass greater than the fiducial region.

  //  cout << "##DEBUG E 1 " << endl; 

  //Get the Fit Function and its limits
  TF1 *fitFunc = spectraFits.at(which);
  Double_t minFitRange(0), maxFitRange(0);
  fitFunc->GetRange(minFitRange,maxFitRange);
  
  //  cout << "##DEBUG E 2 " << endl; 

  //Get the maximum mT-m0 of the fiducial region
  //  Double_t maxFiducialmTm0 = GetX()[nPoints-1] + binWidths.at(nPoints-1)/2.0;
  Double_t maxFiducialmTm0 = 1.0;

  //  cout << "##DEBUG E 3 " << endl; 

  //Compute the dNdy in the range [maxFiducialmTm0,maxFitRange]
  std::pair<double,double> dNdyTemp = ComputedNdyInRange(which,maxFiducialmTm0,maxFitRange);
  
  //  cout << "##DEBUG E 4 " << endl; 

  dNdyHighExtrap.at(which) = dNdyTemp; 
  
}

//_______________________________________________________________
void SpectraClass::ComputedNdyLowExtrapolation(int which){

  //Compute the amount of dNdy that results from the extrapolation
  //of the fit function to transverse mass less than the fiducial region
  //  cout << "##DEBUG F 1 " << endl; 
  //Get the Fit Function and its limits
  TF1 *fitFunc = spectraFits.at(which);
  Double_t minFitRange(0), maxFitRange(0);
  fitFunc->GetRange(minFitRange,maxFitRange);

  //  cout << "##DEBUG F 2 " << endl; 

  //Get the minimum mT-m0 of the fiducial region
  //  Double_t minFiducialmTm0 = GetX()[0] - binWidths.at(0)/2.0;
  Double_t minFiducialmTm0 = 0.2;

  //  cout << "##DEBUG F 3 " << endl; 

  //Compute the dNdy in the range [minFitRange,minFiducialRange]
  std::pair<double,double> dNdyTemp = ComputedNdyInRange(which,minFitRange,minFiducialmTm0);

  //  cout << "##DEBUG F 4 " << endl; 

  dNdyLowExtrap.at(which) = dNdyTemp; 
}

//____________________________________________________________
std::pair<double,double> SpectraClass::ComputeTotaldNdy(TF1 *thisFunc, TMatrixDSym *thisCovMatrix,
							double min, double max){

  if (!thisFunc){
    cout <<"WARNING - SpectraClass::ComputeTotaldNdy() - The Function does not exist. Returning zeros!\n";
    return std::make_pair(0,0);
  }

  //  cout << "****DEBUG D 1" << endl;

  Double_t functionRangeMin(0), functionRangeMax(0);
  thisFunc->GetRange(functionRangeMin,functionRangeMax);
  Double_t stepWidth = fabs(functionRangeMax-functionRangeMin)/
    (double) thisFunc->GetNpx();

  //  cout << "****DEBUG D 2" << endl;

  Double_t mTm0 = min;
  Double_t sum(0);
  Double_t sumErr(0);
  Double_t pMass = particleInfo.GetParticleMass(pid);

  //  cout << "****DEBUG D 3" << endl;

  while (mTm0 <= max){

    Double_t mT = mTm0 + pMass;
    sum += mT * thisFunc->Integral(mTm0,mTm0+stepWidth);

    if (thisCovMatrix)
      sumErr =0;
      //      sumErr += TMath::Power(mT*thisFunc->IntegralError(mTm0,mTm0+stepWidth,
      //						thisFunc->GetParameters(),
      //						thisCovMatrix->GetMatrixArray()),2);
     mTm0 += stepWidth;
  }

  //  cout << "****DEBUG D 4" << endl;

  Double_t preFactors = 2.0*TMath::Pi();
  
  Double_t dNdy = preFactors * sum;
  Double_t dNdyErr = preFactors * sqrt(sumErr);
   
  //  cout << "****DEBUG D 5" << endl;

  return std::make_pair(dNdy,dNdyErr);  

  
}

//____________________________________________________________
std::pair<double,double> SpectraClass::ComputeTotaldNdy(int which){

  return ComputedNdyInRange(which);
  
}

//____________________________________________________________
std::pair<double,double> SpectraClass::ComputedNdyInRange(int which, double min, double max){

  //Check which spectrum is requested, assume first fit if default
  if (which < 0)
    which = 0;

  //  cout << "compute total" << endl;

  return ComputeTotaldNdy(spectraFits.at(which),spectraFitCovarianceMatrix.at(which),min,max);
  
}

//____________________________________________________________
std::pair<double,double> SpectraClass::ComputeMeanmTm0(int whichFunc){

  //Return the mean transverse mass of the spectrum. It is
  //computed as a weighted mean.

  //Always use the default unless specified
  if (whichFunc < 0)
    whichFunc = 0;

  TF1 *fitFunc = spectraFits.at(whichFunc);

  Double_t functionRangeMin(0), functionRangeMax(0);
  fitFunc->GetRange(functionRangeMin,functionRangeMax);
  Double_t stepWidth = fabs(functionRangeMax-functionRangeMin)/
    (double) fitFunc->GetNpx();
  
  Double_t sum(0);
  Double_t sumWeights(0);
  Double_t sumErr(0); //No statistical error is ascribed (all systematic)
  Double_t mTm0(functionRangeMin);
  Double_t pMass = particleInfo.GetParticleMass(pid);
  Double_t factor = 2.0*TMath::Pi() * stepWidth;

  while (mTm0 < functionRangeMax){
    Double_t mT = mTm0 + pMass;
    Double_t weight = mT * fitFunc->Eval(mTm0) * factor;

    sum += weight * mTm0;
    sumWeights += weight;
    
    mTm0 += stepWidth;
  }

  std::pair<double,double> vals = std::make_pair(sum/sumWeights, sqrt(sumErr));
  
  return vals;

}

//______________________________________________________________________________
std::pair<double,double> SpectraClass::ComputeFiducialdNdy(double min, double max){

  //Return the dNdy by summing up the yield in the spectra from min to max.
  //If min and max are both set then the sum is done only in the range [min, max].
  //If only min is set then the sum is done from min to the last point in the spectrum.
  //If only max is set then the sum is done from the first point to max.

  //Check to make sure the binWidths vector is set
  if (binWidths.size() != (unsigned int)nPoints){
    cout <<"ERROR - SpectraClass::ComputeFiducialdNdy() - the BinWidths vector is not currently set properly.\n You may need to call SpectraClass::SetBinWidthVector() before calling this function or if you have\n already called it you may have filled the vector with the incorrect number of values. EXITING! \n";
    cout <<binWidths.size() <<" " <<nPoints <<endl;
    exit (EXIT_FAILURE);
  }

  double dNdy(0);
  double dNdyErr(0);
  double errSum(0);

  double pidMass = particleInfo.GetParticleMass(pid);

  for (int iPoint=0; iPoint < nPoints; iPoint++){

    double mTm0 = GetX()[iPoint];
    if (min > 0 && mTm0 < min)
      continue;
    if (max > 0 && mTm0 > max)
      continue;

    double mT = mTm0 + pidMass;
    double yield = GetY()[iPoint];
    double yieldErr = GetEY()[iPoint];
    double binWidth = binWidths.at(iPoint);
    
    dNdy += yield *binWidth * mT;
    errSum += pow(yieldErr * binWidth * mT, 2);
      
  }//End Loop Over points of the spectrum
  
  dNdyErr = sqrt(errSum);

  //Apply the Normalization Factors
  dNdy = dNdy * 2.0 * TMath::Pi();
  dNdyErr = dNdyErr * 2.0 * TMath::Pi();

  //Set the Interval Value
  if (dNdyFiducial.first < 0.001){
    dNdyFiducial.first = dNdy;
    dNdyFiducial.second = dNdyErr;
  }

  return std::make_pair(dNdy,dNdyErr);
  
  
}


/**************************************************************************
FUNCTIONS USEFUL FOR GRAPHING
 **************************************************************************/
TGraphErrors *SpectraClass::GetTGraph(){

  TGraphErrors *g = new TGraphErrors(nPoints,xVals,yVals,exVals,eyVals);
  g->SetName(name.Data());  
  g->SetMarkerStyle(particleInfo.GetParticleMarker(pid,charge));
  g->SetMarkerColor(particleInfo.GetParticleColor(pid));

  return g;
}


//____________________________________________________________________
TMultiGraph *SpectraClass::GetScaledCopy(double factor, int color){

  //Return a copy of the spectrum which has been scaled by factor
  TMultiGraph *multiG = new TMultiGraph();
  TGraphErrors *gr1 = NULL;
  TGraphErrors *gr2 = NULL;
  
  if (tpcSpectrum){
    gr1 = new TGraphErrors();
    gr1->SetMarkerStyle(tpcSpectrum->GetMarkerStyle());
    if (color < 0)
      gr1->SetMarkerColor(tpcSpectrum->GetMarkerColor());
    else
      gr1->SetMarkerColor(color);

    for (int iPoint=0; iPoint<tpcSpectrum->GetN(); iPoint++){
      gr1->SetPoint(gr1->GetN(),tpcSpectrum->GetX()[iPoint],
		    tpcSpectrum->GetY()[iPoint] * factor);
      gr1->SetPointError(gr1->GetN()-1,tpcSpectrum->GetEX()[iPoint],
			 tpcSpectrum->GetEY()[iPoint] * factor);
    }
    multiG->Add(gr1);
  }//End Copy tpcSpectrum

  if (tofSpectrum){
    gr2 = new TGraphErrors();
    gr2->SetMarkerStyle(tofSpectrum->GetMarkerStyle());
    if (color < 0)
      gr2->SetMarkerColor(tofSpectrum->GetMarkerColor());
    else
      gr2->SetMarkerColor(color);

    for (int iPoint=0; iPoint<tofSpectrum->GetN(); iPoint++){
      gr2->SetPoint(gr2->GetN(),tofSpectrum->GetX()[iPoint],
		    tofSpectrum->GetY()[iPoint] * factor);
      gr2->SetPointError(gr2->GetN()-1,tofSpectrum->GetEX()[iPoint],
			 tofSpectrum->GetEY()[iPoint] * factor);
    }
    multiG->Add(gr2);
  }//End Copy TOF Spectru

  return multiG;

}

//_________________________________________________________________
TF1 *SpectraClass::GetScaledFunction(double factor, unsigned int which, int color){

  //Return a pointer to the scaled version of the fit

  TF1 *orig = GetSpectrumFit(which);
  Double_t min(0), max(0);
  orig->GetRange(min,max);

  TF1 *scaled = new TF1(Form("%s_Scaled",orig->GetName()),
			this,&SpectraClass::scaleFunction,min,max,2);
  scaled->FixParameter(0,factor);
  scaled->FixParameter(1,which);
  scaled->SetLineWidth(3);
 
  if (color > 0)
    scaled->SetLineColor(color);
  else
    scaled->SetLineColor(kBlack);

  return scaled;
  
}

//_________________________________________________________________
Double_t SpectraClass::scaleFunction(double *x, double *par){

  Double_t xVal = x[0];

  Double_t scaleFactor = par[0];
  int which = TMath::Nint(par[1]);
  
  //Since we just converted a double to a integrer
  //there is some chance that if par[1] is close to zero
  //it ould be converted to a negative number. If this was the
  //case then assume it should be zero. 
  if (which < 0)
    which = 0;

  TF1 *orig = spectraFits.at(which);
  if (!orig)
    return 1;
    
  return orig->Eval(xVal) * scaleFactor; 

}

//________________________________________________________________
void SpectraClass::SetMarkerColor(Color_t colorTPC, Color_t colorTOF){

  if (colorTOF < 0)
    colorTOF = colorTPC;

  if (tpcSpectrum)
    tpcSpectrum->SetMarkerColor(colorTPC);
  if (tofSpectrum)
    tofSpectrum->SetMarkerColor(colorTOF);

  
}

//_________________________________________________________________
void SpectraClass::SetMarkerSize(Size_t sizeTPC, Size_t sizeTOF){

  if (sizeTOF < 0)
    sizeTOF = sizeTPC;

  if (tpcSpectrum)
    tpcSpectrum->SetMarkerSize(sizeTPC);
  if (tofSpectrum)
    tofSpectrum->SetMarkerSize(sizeTOF);
  
}

//__________________________________________________________________
void SpectraClass::SetMarkerStyle(Style_t styleTPC, Style_t styleTOF){

  if (styleTOF < 0)
    styleTOF = styleTPC;

  if (tpcSpectrum)
    tpcSpectrum->SetMarkerStyle(styleTPC);
  if (tofSpectrum)
    tofSpectrum->SetMarkerStyle(styleTOF);
  
}
