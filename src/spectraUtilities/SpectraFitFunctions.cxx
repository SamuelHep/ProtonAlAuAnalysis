//Implementation of spectra and ratio fit functions.
//Unless otherwise noted below, all functions were implemented
//by Christopher Flores (chrflores@ucdavis.edu)

#include <utility>
#include <vector>

#include <TString.h>
#include <TMath.h>
#include <TF1.h>

#include "utilityFunctions.h"
#include "ParticleInfo.h"
#include "SpectraClass.h"

#include "SpectraFitFunctions.h"

/******************************************************************
   PION RATIO COULOMB FIT AND FIT RESULT PRINT FUNCTIONS
/******************************************************************/
//_________________________________________________________
Double_t PionRatioFit(Double_t *x, Double_t *par){

  //NOTE: Parameters 2,3,and 4 should be FIXED in the following ways:
  //      Par2 = FIXED to the PION slope parameter from spectra fit
  //      Par3 = FIXED to the PROTON slope parameter from spectra fit
  //      Par4 = FIXED to the recentered rapidity, y-yMid

  ParticleInfo particleInfo;
  
  Double_t mTm0 = x[0];

  //Fit Parameters
  Double_t initialRatio = par[0];
  Double_t coulombV = par[1]/1000.0;
  Double_t pionTemp = par[2];  //This Should be Fixed the pion Spectra Fit
  Double_t protonTemp = par[3];//This Should be Fixed from the proton Spectra Fit
  Double_t y = par[4];         //This Should be Fixed (y-yMid)

  //Additional Variables for -->Pion<--
  Double_t m0 = particleInfo.GetParticleMass(PION);
  Double_t mT = mTm0 + m0;
  Double_t energy = mT*cosh(y);
  Double_t p = sqrt(energy*energy - m0*m0);
  Double_t betaGamma = p / m0;

  //Additional Variables for -->Proton<--
  Double_t mP = particleInfo.GetParticleMass(PROTON);

  //Compute the Effective Coulomb Potential
  Double_t eMax = sqrt(pow(mP*betaGamma,2) + pow(mP,2)) - mP;
  Double_t effCoulombV = coulombV * (1.0 - exp(-eMax/protonTemp));


  //Compute the Transformation Jacobian
  Double_t num = sqrt(pow(energy - effCoulombV,2) - pow(m0,2));
  Double_t denom = sqrt(pow(energy + effCoulombV,2) - pow(m0,2));
  Double_t Jacobian = (energy - effCoulombV) / (energy + effCoulombV);
  Jacobian = Jacobian * ( num/denom );

  //Compute the Emmision Factor for Bose Einstein
  Double_t emmisionFactor = exp( (energy + effCoulombV)/pionTemp ) - 1.0;
  emmisionFactor = emmisionFactor / (exp( (energy - effCoulombV)/pionTemp) - 1.0);

  //Assemble it all Together
  return initialRatio * emmisionFactor * Jacobian;


}

//________________________________________________________
TString PrintPionRatioFitResults(TF1 *fitFunc){
  
  Double_t initialPionRatio     = fitFunc->GetParameter(0);
  Double_t initialPionRatioErr = fitFunc->GetParError(0);

  Double_t coulombPotential    = fitFunc->GetParameter(1);
  Double_t coulombPotentialErr = fitFunc->GetParError(1);

  Double_t chi2 = fitFunc->GetChisquare();
  Double_t ndf = fitFunc->GetNDF();

  return TString::Format("%.04f +- %.04f\t%.04f +- %.04f\t%.04f",
                         coulombPotential,coulombPotentialErr,
			 initialPionRatio,initialPionRatioErr,
			 chi2/ndf);

}

/******************************************************************
   INTEGRATE THE UNINVERTED FIT CURVE TO FIND THE DNDY
 ******************************************************************/
std::pair<double, double> ComputedNdyFromFit(TF1 *fitFunc,
					     TGraphErrors *confInterval,
					     Double_t pMass,
					     Double_t minRange, Double_t maxRange){

  //This function returns the dNdy of the fitFunc in the range
  //[minRange,maxRange). By Default minRange=0 and maxRange=10.0

  //NOTE: The TF1 passed into this function must be the UNINVERTED
  //      form of the fit. That is, the zeroth parameter must be a
  //      generic amplitude and NOT the dNdy already.

  //NOTE: This function will only return the correct dNdy if your
  //      spectra have a (1/2pi)*(1/mT) phase factor and are
  //      a function of mT-m0.

  Double_t functionRangeMin(0), functionRangeMax(0);
  fitFunc->GetRange(functionRangeMin,functionRangeMax);
  Double_t stepWidth = fabs(functionRangeMax-functionRangeMin)/
    (double) fitFunc->GetNpx();

  //Integrate using simple rectangles at three different points
  // 0. Begining of rectangle
  // 1. Middle of rectangle
  // 2. End of rectangle

  std::vector<double> sum(3,0);
  std::vector<double> sumErr(3,0);
  
  // 0. Beginning of Rectangle
  Double_t mTm0 = minRange;
  while (mTm0 < maxRange){

    Double_t mT = mTm0 + pMass;
    sum.at(0) += mT * fitFunc->Eval(mTm0);
    sumErr.at(0) += TMath::Power(mT*GetErrorAtValue(confInterval,mTm0),2);
    mTm0 += stepWidth;
    
  }

  // 1. Middle of Rectangle
  mTm0 = minRange + stepWidth/2.0;
  while (mTm0 < maxRange){
    Double_t mT = mTm0 + pMass;
    sum.at(1) += mT * fitFunc->Eval(mTm0);
    sumErr.at(1) += TMath::Power(mT*GetErrorAtValue(confInterval,mTm0),2);
    mTm0 += stepWidth;
  }

  // 2. End of Rectangle
  mTm0 = minRange + stepWidth;
  while (mTm0 < maxRange){
    Double_t mT = mTm0 + pMass;
    sum.at(2) += mT * fitFunc->Eval(mTm0);
    sumErr.at(2) += TMath::Power(mT*GetErrorAtValue(confInterval,mTm0),2);
    mTm0 += stepWidth;
  }

  Double_t avgSum = TMath::Mean(3,&sum.at(0));
  Double_t rmsSum = TMath::RMS(3,&sum.at(0));
  for (int i=0; i<3; i++){
    sumErr.at(i) = sqrt(sumErr.at(i));
  }

  Double_t preFactors = 2.0*TMath::Pi() * stepWidth;
  
  Double_t dNdy = preFactors * avgSum;
  Double_t dNdyErr = preFactors *
    ( sqrt(TMath::MaxElement(3,&sumErr.at(0))) + rmsSum);

  return std::make_pair(dNdy,dNdyErr); 

}

/******************************************************************
   COUNT THE DNDY IN THE MEASURED REGION OF THE SPECTRUM
/******************************************************************/
//__________________________________________________________
std::pair<double,double> CountdNdyOfSpectrum(TGraphErrors *spectrum, Int_t pid,
					     Double_t minRange, Double_t maxRange){

  //Count the dNdy of a measured spectrum by numerically integrating
  //its points in the range from minRange to maxRange. Return a pair
  //corresponding to the dNdy and its error.
  
  Double_t dNdy(0);
  Double_t dNdyErr(0);
  Double_t errSum(0);
    
  ParticleInfo particleInfo;
  Double_t pMass = particleInfo.GetParticleMass(pid);

  for (Int_t iPoint=0; iPoint<spectrum->GetN(); iPoint++){

    Double_t mTm0 = spectrum->GetX()[iPoint];
    if (mTm0 < minRange || mTm0 > maxRange)
      continue;

    Double_t mT = mTm0 + pMass;
    Double_t yield = spectrum->GetY()[iPoint];
    Double_t yieldErr = spectrum->GetEY()[iPoint];
    Double_t binWidth = spectrum->GetEX()[iPoint]*2.0; //Error is half of bin width
    
    dNdy += yield * binWidth * mT;
    errSum += pow(yieldErr * binWidth * mT,2);
    
  }

  //Compute the Error Before the 2Pi Normalization
  dNdyErr = sqrt(errSum);

  //Apply the Normalization  
  dNdy = dNdy * 2.0 * TMath::Pi();
  dNdyErr = dNdyErr * 2.0 * TMath::Pi();

  return std::make_pair(dNdy,dNdyErr);

}

/******************************************************************
   BOSE-EINSTEIN SPECTRAL FIT FUNCTION AND NUMERICAL INTEGRATION
/******************************************************************/
//__________________________________________________________
Double_t BoseEinsteinKernel(Double_t *x, Double_t *par){

  //This is not a fit function, it is the functional
  //dependence of the Bose Einstein function and used below.

  Double_t mTm0 = x[0];

  Double_t slope = par[1];
  Double_t pMass = par[2];

  Double_t mT = mTm0 + pMass;

  return 1.0 / (exp(mT/slope) - 1.0);
  
}

//__________________________________________________________
Double_t InvertBoseEinstein(Double_t *par,
			    Double_t minRange, Double_t maxRange){
  
  //Returns the integral of mT*BoseEinstein*dmT.
  //Here we integrate the fitting function over dmT
  //from minRange to maxRange GeV in .01 GeV Steps
  
  Double_t stepWidth = .01;
  Int_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum = 0.0;
  Double_t mTm0Arr[] = {minRange}; //mT-m0
  Double_t pMass = par[2];
  Double_t mT(0);
  
  for (Int_t stepIndex=0; stepIndex < nSteps; stepIndex++){
    
    mT = mTm0Arr[0]+pMass;

    sum += mT * BoseEinsteinKernel(mTm0Arr,par) * stepWidth;
    
    mTm0Arr[0] += stepWidth;

  } 

  return sum;
}

//_______________________________________________________
Double_t BoseEinsteinFitFunc(Double_t *x, Double_t *par){

  //This is the Bose Einstein Function which is used to fit the spectra
  //It consists of the normalized part which contains the dNdy and the 
  //BE Functional Form

  //Parameter Definitions
  //par0 - the normalization
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!

  //Get the dNdy Parameter
  Double_t norm = par[0];

  return (1.0/(2.0*TMath::Pi())) * norm
    * BoseEinsteinKernel(x,par);

}

//_______________________________________________________
Double_t BoseEinsteinFitFuncInRange(Double_t *x, Double_t *par){

  //In this case the dNdy (parameter 0) is the dNdy in
  //the range of [minRange,maxRange].
  
  //Parameter Definitions
  //par0 - the dNdy
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!
  //par3 - minRange      - SHOULD BE FIXED!
  //par4 - maxRange      - SHOULD BE FIXED!

  Double_t dNdy     = par[0];
  Double_t tSlope   = par[1];
  Double_t pMass    = par[2]; //This Should be FIXED!
  Double_t minRange = par[3]; //This Should be FIXED!
  Double_t maxRange = par[4]; //This Should be FIXED!

  return (1.0/(2.0*TMath::Pi()))
    * dNdy
    * (1.0 / InvertBoseEinstein(par,minRange,maxRange))
    * BoseEinsteinKernel(x,par);

}

//_______________________________________________________
Double_t BoseEinsteinFitFuncWithInversion(Double_t *x, Double_t *par){

  //In this case the zeroth parameter is a generic amplidutde
  //parameter and NOT the dNdy

  //Parameter Definitions
  //par0 - the dNdy
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!
  
  return (1.0/(2.0*TMath::Pi()))
    * par[0]
    * (1.0 / InvertBoseEinstein(par))
    * BoseEinsteinKernel(x,par);

}

/******************************************************************
   MT-EXPONENTIAL SPECTRAL FIT FUNCTION AND NUMERICAL INTEGRATION
/******************************************************************/
//___________________________________________________________
Double_t pTExponentialKernel(Double_t *x, Double_t *par){

  //This is not a fit function, it is the functional
  //dependence of the pT Exponential function and used below

  Double_t mTm0 = x[0];

  Double_t slope = par[1];
  Double_t pMass = par[2];

  Double_t pT = ConvertmTm0ToPt(mTm0,pMass);

  return exp(-pT/slope);
  
}

//___________________________________________________________
Double_t pTExponentialFitFunc(Double_t *x, Double_t *par){

  //Parameter Definitions
  //par0 = Normalization
  //par1 = Slope
  //par2 = Particle Mass (Should be Fixed)

  return (1.0/(2.0*TMath::Pi())) * par[0] * pTExponentialKernel(x,par);
  
}

/******************************************************************
   MT-EXPONENTIAL SPECTRAL FIT FUNCTION AND NUMERICAL INTEGRATION
/******************************************************************/
//___________________________________________________________
Double_t mTExponentialKernel(Double_t *x, Double_t *par){

  //This is not a fit function, it is the functional
  //dependence of the mT Exponential function and used below

  Double_t mTm0 = x[0];

  Double_t slope = par[1];
  Double_t pMass = par[2];

  Double_t mT = mTm0 + pMass;

  return 1.0 / exp(mTm0/slope);
  
}

//___________________________________________________________
Double_t InvertmTExponential(Double_t *par,
			     Double_t minRange, Double_t maxRange){
  
  //Returns the integral of the mT*mTExponential*dmT fit by numerically
  //inverting the function in the specified range.

  Double_t stepWidth = .01;
  Int_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum = 0.0;
  Double_t mTm0Arr[] = {minRange}; //mT-m0
  Double_t pMass = par[2];
  Double_t mT(0);
  
  for (Int_t stepIndex=0; stepIndex < nSteps; stepIndex++){
    
    mT = mTm0Arr[0]+pMass;
    
    sum += mT * mTExponentialKernel(mTm0Arr,par) * stepWidth;
    
    mTm0Arr[0] += stepWidth;
   
  } 

  return sum;
}

//____________________________________________________________
Double_t mTExponentialFitFuncInRange(Double_t *x, Double_t *par){

  //The mTExponential funtion. The zeroth parameter is the dNdy,
  //from mT-m0=[minRange,maxRange] . See inversion routine above

  //Parameter Definitions
  //par0 - the dNdy
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!
  //par3 - minRange      - SHOULD BE FIXED!
  //par4 - maxRange      - SHOULD BE FIXED!
  
  //Get the Parameters
  Double_t dNdy     = par[0];
  Double_t tSlope   = par[1]; 
  Double_t pMass    = par[2]; //This should be FIXED
  Double_t minRange = par[3]; //This Should be FIXED!
  Double_t maxRange = par[4]; //This Should be FIXED!

  return (1.0/(2.0*TMath::Pi()))
    * dNdy
    * (1.0 / InvertmTExponential(par,minRange,maxRange) )
    * mTExponentialKernel(x,par);
    
}

//____________________________________________________________
Double_t mTExponentialFitFuncWithInversion(Double_t *x, Double_t *par){

  //The mTExponential funtion. The zeroth parameter is the dNdy,
  //from mT-m0=0 to "infinity". See inversion routine above

  //Parameter Definitions
  //par0 - the dNdy
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!
  
  //Get the Parameters
  Double_t dNdy = par[0];
  Double_t tSlope = par[1];   
  Double_t pMass = par[2];    //This should be FIXED

  return (1.0/(2.0*TMath::Pi()))
    * dNdy
    * (1.0 / InvertmTExponential(par) )
    * mTExponentialKernel(x,par);
  
}

//____________________________________________________________
Double_t mTExponentialFitFunc(Double_t *x, Double_t *par){

  //This is the mTExponential fit function which is used to fit the spectra.
  //It consists of the normalized part which contains the dNdy and
  //the exponential part.

  //Parameter Definitions
  //par0 - Generic Amplitude
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!
  
  return (1.0/(2.0*TMath::Pi())) * par[0] * mTExponentialKernel(x,par);
  
}


/******************************************************************
   MAXWWELL-BOLTZMANN FIT FUNCTION
/******************************************************************/
Double_t BoltzmannKernel(Double_t *x, Double_t *par){

  //This is not a fit function, it is the functional
  //dependence of the Boltzmann Function and used below.

  Double_t mTm0 = x[0];

  Double_t slope = par[1];
  Double_t pMass = par[2];

  Double_t mT = mTm0 + pMass;

  return mT * (1.0 / exp(mT/slope) );  

}

//_________________________________________________________________________
Double_t InvertBoltzmann(Double_t *par,
			 Double_t minRange, Double_t maxRange){

  //Returns the integral of mT*Boltzmann*dmT.
  //Here we integrate the fitting function of dmT
  //from minRange to maxRange.
  
  Double_t stepWidth = .01;
  Int_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum = 0.0;
  Double_t mTm0Arr[] = {minRange}; //mT-m0
  Double_t pMass = par[2];
  Double_t mT(0);

  for (Int_t stepIndex=0; stepIndex < nSteps; stepIndex++){

    mT = mTm0Arr[0] + pMass;

    sum += mT * BoltzmannKernel(mTm0Arr,par) * stepWidth;

    mTm0Arr[0] += stepWidth;

  }

  return sum;
  
}

//_________________________________________________________________________
Double_t BoltzmannFitFunc(Double_t *x, Double_t *par){

  //Boltzmann Distribution fit function. The zeroth parameter is the dNdy,
  //obtained by doing the integral from mT-m0=0 to "infinity"

  //Parameter Definitions
  //par0 - normalization
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!

  //Get the Parameters
  Double_t dNdy = par[0];

  return (1.0/(2.0*TMath::Pi()))
    * par[0]
    * BoltzmannKernel(x,par);
}

//_________________________________________________________________________
Double_t BoltzmannFitFuncInRange(Double_t *x, Double_t *par){

  //Boltzmann Distribution fit function. The zeroth parameter is the dNdy,
  //obtained by doing the integral from mT-m0=minRange to maxRange

  //Parameter Definitions
  //par0 - the dNdy
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!
  //par3 - minRange      - SHOULD BE FIXED!
  //par4 - maxRange      - SHOULD BE FIXED!
  
  //Get the Parameters
  Double_t dNdy     = par[0];
  Double_t tSlope   = par[1];
  Double_t pMass    = par[2];    
  Double_t minRange = par[3];  
  Double_t maxRange = par[4]; 
  
  return (1.0/(2.0*TMath::Pi()))
    * dNdy
    * (1.0 / InvertBoltzmann(par,minRange,maxRange) )
    * BoltzmannKernel(x,par);

}

//_________________________________________________________________________
Double_t BoltzmannFitFuncWithInversion(Double_t *x, Double_t *par){

  //Boltzmann Distribution Fit Function with generic amplitude
  //parameter that is the dNdy.

  //Parameter Definitions
  //par0 - the dNdy
  //par1 - the temperature slope
  //par2 - particle mass - SHOULD BE FIXED!

  //Get the Parameters
  Double_t dNdy = par[0];

  return (1.0/(2.0*TMath::Pi()))
    * dNdy
    * (1.0 / InvertBoltzmann(par))
    * BoltzmannKernel(x,par);

}

/******************************************************************
   THERMAL SPECTRAL FIT FUNCTION (MAXWELL-BOLTZMANN)
   Note: See Jen Klay's Thesis Eqn 7.5
/******************************************************************/
//___________________________________________________________________________
Double_t ThermalFitFunc(Double_t *x, Double_t *par){

  //NOTE: par2 corresponding to the particle mass should be fixed
  
  Double_t xx =x[0];

  Double_t dNdy = par[0];
  Double_t tSlope = par[1];
  Double_t pMass = par[2];  //This Should be Fixed

  return (1.0/(2.0*TMath::Pi())) * (dNdy/tSlope) *(xx+pMass) *(1.0/(2.0*tSlope*tSlope + 2.0*tSlope*pMass + pMass*pMass)) * exp(-(xx)/tSlope);

}

/******************************************************************
   DOUBLE EXPONENTIAL FIT FUNCTION 
   Note: See Lokesh's Mid Rapidity Spectra Paper
/******************************************************************/
//___________________________________________________________________
Double_t InvertDoubleExponential(Double_t A1, Double_t T1, Double_t A2, Double_t T2,
				 Double_t pMass,Double_t minRange, Double_t maxRange){

  //Returns the integral of the Double Exponential fit by numerically
  //inverting the function in the specified range.

  Double_t stepWidth = .01;
  Int_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum = 0.0;
  Double_t x = minRange; //mT-m0

  //Integrals for each term

  //Term 1
  Double_t term1(0);
  for (Int_t stepIndex=0; stepIndex < nSteps; stepIndex++){
    
    term1 += ( (x+pMass)*exp(-(x+pMass)/T1) ) * stepWidth;
    
    //Increment x
    x += stepWidth;
    
  } 

  term1 = A1 * term1;

  //Term 2
  x=minRange;
  Double_t term2(0);
  for (Int_t stepIndex=0; stepIndex < nSteps; stepIndex++){
    
    term2 += ( (x+pMass)*exp(-(x+pMass)/T2) ) * stepWidth;
    
    //Increment x
    x += stepWidth;
    
  } 

  term2 = A2 * term2;

  sum = term1 + term2;
  
  return sum;

}

//_________________________________________________________________
Double_t DoubleExponential(Double_t *x, Double_t *par){

  Double_t xx=x[0]; //mT-m0

  Double_t dNdy  = par[0];
  Double_t A1    = par[1];
  Double_t T1    = par[2];
  Double_t A2    = par[3];
  Double_t T2    = par[4];
  Double_t pMass = par[5]; //This should be fixed

  return (1.0/(2.0*TMath::Pi()))*(dNdy)*( (A1*exp(-(xx+pMass)/T1)+A2*exp(-(xx+pMass)/T2)) / InvertDoubleExponential(A1,T1,A2,T2,pMass));

}


/******************************************************************
   SIEMENS AND RASMUSSEN MODEL
   See: Phys rev 42,880 (1979)

   Original Implementation in Fortran by Mike Heffner
   Modified by Daniel Cebra to use mT
   Converted into ROOT/C++ here
/******************************************************************/
Double_t SiemensRasmussenKernel(Double_t *x, Double_t *par){

  //This is not a fit function. It is the unnormalized functional dependence
  //of the Siemens-Rasmussen Fit Function and is used in the fit functions below

  Double_t mTm0 = x[0];

  Double_t T     = par[1];
  Double_t beta  = par[2];
  Double_t pMass = par[3];

  Double_t mt = mTm0+pMass;
  Double_t E  = mt;
  Double_t p  = sqrt(pow(E,2)-pow(pMass,2));
  Double_t gamma = 1.0/sqrt(1.0-pow(beta,2));
  Double_t alpha = (gamma*beta*p) / T;

  Double_t A = exp(-gamma*E/T);
  Double_t B = (sinh(alpha)/alpha)*(gamma+(T/E));
  Double_t C = (T/E)*cosh(alpha);
  
  return (1.0/(2.0*TMath::Pi())) * A * (B - C);
}

//___________________________________________________________________
Double_t InvertSiemensRasmussen(Double_t T, Double_t beta, 
				Double_t pMass,Double_t minRange, Double_t maxRange){

  //Returns the integral of the Siemens Rasmussen Function in
  //the specified range
  
  Double_t stepWidth = .01;
  Int_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum = 0.0;
  Double_t x = minRange+0.00001; //mT-m0

  for (Int_t stepIndex=0; stepIndex < nSteps; stepIndex++){

    Double_t mt = x+pMass;
    Double_t E = mt;
    Double_t p = sqrt(pow(E,2)-pow(pMass,2));
    Double_t gamma = 1.0/sqrt(1.0-pow(beta,2));
    Double_t alpha = (gamma*beta*p) / T;

    Double_t A = exp(-gamma*E/T);
    Double_t B = (sinh(alpha)/alpha)*(gamma+(T/E));
    Double_t C = (T/E)*cosh(alpha);

    sum += ( E * A * (B-C) ) * stepWidth;
    
    //Increment x
    x += stepWidth;
    
  } 

  return sum;

  
}

//___________________________________________________________________
Double_t SiemensRasmussen(Double_t *x, Double_t *par){

  Double_t xx=x[0]; //mT-m0

  Double_t A  = par[0];
  
  return A * SiemensRasmussenKernel(x,par);

}

/******************************************************************
   BLAST WAVE MODEL
  


/******************************************************************/
//___________________________________________________________________
Double_t RadialIntegralBW(Double_t betaS, Double_t n, Double_t T, Double_t mTm0,
			  Double_t mass, Double_t R){

  //Radial Integral for the Blast Wave Model

  //Integration Variable
  Double_t r(0.0); //The lower limit of integration is 0.0

  //Integration Steps
  Double_t sum(0.0);
  Double_t stepWidth(0.01);
  Int_t    nSteps = (Int_t)(fabs(R)/stepWidth);

  //Compute mT and pT
  Double_t mT = mTm0 + mass;
  Double_t pT = ConvertmTm0ToPt(mTm0,mass);

  //Do the numerical integral
  for (int iStep=0; iStep<nSteps; iStep++){

    Double_t betaT = betaS * pow(r/R,n);
    Double_t rho   = TMath::ATanH(betaT);

    Double_t BesselI0 = TMath::BesselI0((pT/T)*TMath::SinH(rho));
    Double_t BesselK1 = TMath::BesselK1((mT/T)*TMath::CosH(rho));

    //Increment the Sum
    sum += r * mT * BesselI0 * BesselK1 * stepWidth;

    //Increment the integration variable
    r += stepWidth;
    
  }//End Loop Over iStep

  return sum;
  
}

//___________________________________________________________________
Double_t InvertBlastWaveFit(Double_t betaS, Double_t n, Double_t T,
			    Double_t mass, Double_t R, Double_t minRange, Double_t maxRange){

  //Integrate the Blast Wave Fit in range from minRange to maxRange

  Double_t stepWidth(0.01);
  Double_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum(0.0);
  Double_t x(minRange); //Integration variable (mT-m0)

  for (int iStep=0; iStep<nSteps; iStep++){

    Double_t mT = x+mass;
    
    //Increment the Sum
    sum += mT * RadialIntegralBW(betaS,n,T,x,mass,R) * stepWidth;
    
    //Increment the Integration Variable
    x += stepWidth;
    
  }
  
  return sum;
  
}

//___________________________________________________________________
Double_t BlastWaveModelFit(Double_t *x, Double_t *par){

  Double_t xx=x[0]; //mT-m0

  Double_t norm  = par[0]; //Yield
  Double_t T     = par[1]; //Tkin
  Double_t betaS = par[2]; //Surface Velocity
  Double_t n     = par[3]; //Power in the Tranverse Velocity factor (betaT)
  Double_t R     = par[4]; //Upper Limit of Radial Integral (SHOULD BE FIXED)
  Double_t mass  = par[5]; //Should be FIXED

  return (1.0/(2.0*TMath::Pi())) * norm * RadialIntegralBW(betaS,n,T,xx,mass,R);// / InvertBlastWaveFit(betaS,n,T,mass,R);

}

/*****************************************************************
Bose-Einstein + Blast Wave Model Fit
*****************************************************************/
Double_t BoseBlastWaveFit(Double_t *x, Double_t *par){

  Double_t blastWaveAmp = par[0];
  Double_t blastWaveTemp = par[1];
  Double_t blastWaveBetaS = par[2];
  Double_t blastWaveN     = par[3];
  Double_t blastWaveR     = par[4];

  Double_t particleMass = par[5];

  Double_t boseAmp = par[6];
  Double_t boseSlope = par[7];

  
  Double_t boseParams[3] = {boseAmp,boseSlope,particleMass};
  
  return (1.0/(2.0*TMath::Pi())) * (boseAmp*BoseEinsteinKernel(x,boseParams)
				    +
				    blastWaveAmp*RadialIntegralBW(blastWaveBetaS,blastWaveN,blastWaveTemp,
								  x[0],particleMass,blastWaveR));

  
}


/******************************************************************
Tsallis Distribution
*******************************************************************/
Double_t TsallisKernel(Double_t *x, Double_t *par){

  Double_t mTm0 = x[0];

  Double_t Temp  = par[1];
  Double_t n     = par[2];
  Double_t pMass = par[3];

  return TMath::Power(1.0 + ( (mTm0)/(n*Temp) ),-n);
  
}

//_________________________________________________________
Double_t TsallisFitNoInversion(Double_t *x, Double_t *par){


  return par[0] * TsallisKernel(x,par);

}



/******************************************************************
Piecewise function consisting of Bose-Einstein 
for low mT-m0 and mT Exponential for high mT-m0.
*******************************************************************/
//__________________________________________________________________
Double_t PiecewiseBosemTExpoFitKernel(Double_t *x, Double_t *par){

  //This Function is not a fit function, it is the kernel
  //of the piecewise function that is used in the fits.

  //IMPORTANT!!!!
  //If the stitchVal (par4) is set to something greater
  //than 500 (which is a ridiculously high value for a spectrum fit
  //and must have therefore been set purposfully) this kerenel
  //will result to a standard (non-piecewise) Bose-Einstein Function
  
  Double_t mTm0 = x[0];

  Double_t boseSlope = par[1];
  Double_t expoSlope = par[2];
  Double_t pMass = par[3];
  Double_t stitchVal = par[4];

  Double_t mT = mTm0 + pMass;

  Double_t boseParams[] = {1,boseSlope,pMass};
  Double_t expoParams[] = {1,expoSlope,pMass};

  Double_t stitchArr[] = {stitchVal};

  if (mTm0 <= stitchVal){

    if (stitchVal > 500)
      return BoseEinsteinKernel(x,boseParams);
    
    return mTExponentialKernel(stitchArr,expoParams)
      * (1.0/BoseEinsteinKernel(stitchArr,boseParams))
      * BoseEinsteinKernel(x,boseParams);
  }
  
  return mTExponentialKernel(x,expoParams);
  

}

//________________________________________________________________________
Double_t InvertPiecewiseBosemTExpoFit(Double_t *par,
				 Double_t minRange, Double_t maxRange){
  
  //Integrate the Default Pion Fit
  Double_t stepWidth(.01);
  Double_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum(0.0);
  Double_t mTm0Arr[] = {minRange};
  Double_t pMass = par[3];
  Double_t mT(0);
  
  for (int iStep=0; iStep<nSteps; iStep++){
 
    mT = mTm0Arr[0]+pMass;
    
    sum += mT * PiecewiseBosemTExpoFitKernel(mTm0Arr,par) * stepWidth;
    
    mTm0Arr[0] += stepWidth;
    
  }//End Loop Over Sum

  return sum;
}

//____________________________________________________________________
Double_t PiecewiseBosemTExpoFitNoInversion(Double_t *x, Double_t *par){

  //Parameter Definitions
  //par0 - Generic Amplitude
  //par1 - the Bose-Einstein temperature slope
  //par2 - the mTExponential temperature slope
  //par3 - particle mass - SHOULD BE FIXED!
  //par4 - point where Bose Function turns into mTExponential
  
  Double_t amplitude = par[0];
  
  return amplitude * PiecewiseBosemTExpoFitKernel(x,par);
  
}

//____________________________________________________________________
Double_t PiecewiseBosemTExpoFit(Double_t *x, Double_t *par){

  //Piecewise function combining Bose-Einstein for the low mTm0
  //portion of the spectrum and an mTExponential for the high
  
  //Parameter Definitions
  //par0 - dNdy
  //par1 - the Bose temperature slope
  //par2 - the mTExponential temperature slope
  //par3 - particle mass - SHOULD BE FIXED!
  //par4 - point where Bose Function turns into mTExponential

  Double_t dNdy = par[0];

  return (1.0/(2.0 * TMath::Pi())) * dNdy
    * PiecewiseBosemTExpoFitKernel(x,par)
    * (1.0/InvertPiecewiseBosemTExpoFit(par));  
  
}

/******************************************************************
Piecewise function consisting of two mT Exponential 
functions, one for the low and high mTm0 regions.
*******************************************************************/
//___________________________________________________________
Double_t PiecewiseDoublemTExponentialKernel(Double_t *x, Double_t *par){

  //This function is not a fit function, it is the kernel
  //of the piecewise double exponential fit

  Double_t mTm0 = x[0];

  Double_t slope1 = par[1];
  Double_t slope2 = par[2];
  Double_t pMass = par[3];
  Double_t stitchVal = par[4];

  Double_t mT = mTm0 + pMass;

  Double_t expo1Params[] = {1,slope1,pMass};
  Double_t expo2Params[] = {1,slope2,pMass};

  Double_t stitchArr[] = {stitchVal};

  if (mTm0 <= stitchVal)
    return mTExponentialKernel(stitchArr,expo2Params)
      * (1.0/mTExponentialKernel(stitchArr,expo1Params))
      * mTExponentialKernel(x,expo1Params);
  
  return mTExponentialKernel(x,expo2Params);
  
}

//____________________________________________________________
Double_t InvertPiecewiseDoublemTExponential(Double_t *par,
					    Double_t minRange,
					    Double_t maxRange){

  //Integrate mT*PiecewiseDoublemTExponential*dmT
  Double_t stepWidth(.01);
  Double_t nSteps = (Int_t)(fabs(maxRange-minRange)/stepWidth);
  Double_t sum(0.0);
  Double_t mTm0Arr[] = {minRange};
  Double_t pMass = par[3];
  Double_t mT(0);
  
  for (int iStep=0; iStep<nSteps; iStep++){
 
    mT = mTm0Arr[0]+pMass;
    
    sum += mT * PiecewiseDoublemTExponentialKernel(mTm0Arr,par) * stepWidth;
    
    mTm0Arr[0] += stepWidth;
    
  }//End Loop Over Sum

  return sum;
}

//_____________________________________________________________
Double_t PiecewiseDoublemTExponentialFitNoInversion(Double_t *x,
						    Double_t *par){

  //Parameter Definitions
  //par0 - Generic Amplitude
  //par1 - the mTExponential1 temperature slope
  //par2 - the mTExponential2 temperature slope
  //par3 - particle mass - SHOULD BE FIXED!
  //par4 - point where the first mTExponential turns into the second
  
  Double_t amplitude = par[0];
  
  return amplitude * PiecewiseDoublemTExponentialKernel(x,par);
  
  
}

//_____________________________________________________________
Double_t PiecewiseDoublemTExponentialFit(Double_t *x, Double_t *par){


  //Parameter Definitions
  //par0 - dNdy
  //par1 - the mTExponential1 temperature slope
  //par2 - the mTExponential2 temperature slope
  //par3 - particle mass - SHOULD BE FIXED!
  //par4 - point where the first mTExponential turns into the second

  Double_t dNdy = par[0];

  return (1.0/(2.0 * TMath::Pi())) * dNdy
    * PiecewiseDoublemTExponentialKernel(x,par)
    * (1.0/InvertPiecewiseDoublemTExponential(par)); 
  
}

