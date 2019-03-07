//Header file for spectra and ratio fit functions and routines
//See src/analysis/SpectraFitFunctions.cxx for implementations

class SpectraClass;

//Pion Ratio Coulomb Fit___________________________________________
Double_t PionRatioFit(Double_t *x, Double_t *par);
TString PrintPionRatioFitResults(TF1 *fitFunc);

//Compute the dNdy in the specified range by integrating___________
std::pair<double,double> ComputedNdyFromFit(TF1 *fitFunc, TGraphErrors *confInterval, Double_t pMass,
					    Double_t minRange=0.0, Double_t maxRange=10.0);

//Compute the dNdy in the measured range by counting_______________
std::pair<double,double> CountdNdyOfSpectrum(TGraphErrors *spectrum, Int_t pid,
					     Double_t minRange=0.0, Double_t maxRange=100.0);

//Bose-Einstein fit_________________________________________________
Double_t BoseEinsteinKernel(Double_t *x, Double_t *par);
Double_t InvertBoseEinstein(Double_t *par, Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t BoseEinsteinFitFunc(Double_t *x, Double_t *par); 
Double_t BoseEinsteinFitFuncInRange(Double_t *x, Double_t *par); 
Double_t BoseEinsteinFitFuncWithInversion(Double_t *x, Double_t *par); 

//MT-Exponential fit_________________________________________________
Double_t mTExponentialKernel(Double_t *x, Double_t *par);
Double_t InvertmTExponential(Double_t *par, Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t mTExponentialFitFuncWithInversion(Double_t *x, Double_t *par); 
Double_t mTExponentialFitFuncInRange(Double_t *x, Double_t *par); 
Double_t mTExponentialFitFunc(Double_t *x, Double_t *par); 

//PT-Exponential Fit_______________________________________________
Double_t pTExponentialKernel(Double_t *x, Double_t *par);
Double_t pTExponentialFitFunc(Double_t *x, Double_t *par);

//Boltzmann Fit____________________________________________________
Double_t BoltzmannKernel(Double_t *x, Double_t *par);
Double_t InvertBoltzmann(Double_t *par, Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t BoltzmannFitFunc(Double_t *x, Double_t *par);
Double_t BoltzmannFitFuncInRange(Double_t *x, Double_t *par);
Double_t BoltzmannFitFuncWithInversion(Double_t *x, Double_t *par);

//Thermal (Maxwell-Boltzmann) fit__________________________________
Double_t ThermalFitFunc(Double_t *x, Double_t *par);

//Double Exponential
Double_t InvertDoubleExponential(Double_t A1, Double_t T1, Double_t A2, Double_t T2,
				 Double_t pMass,Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t DoubleExponential(Double_t *x, Double_t *par);

//Siemens and Rasmussen Model______________________________________
Double_t SiemensRasmussenKernel(Double_t *x, Double_t *par);
Double_t InvertSiemensRasmussen(Double_t T, Double_t beta,
				Double_t pMass,Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t SiemensRasmussen(Double_t *x, Double_t *par);

//Blast Wave Model_________________________________________________
Double_t RadialIntegralBW(Double_t betaS, Double_t n, Double_t T, Double_t mTm0,
			  Double_t mass, Double_t R);
Double_t InvertBlastWaveFit(Double_t betaS, Double_t n, Double_t T, Double_t mass,
			    Double_t R, Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t BlastWaveModelFit(Double_t *x, Double_t *par);

//Tsallis Fit Func_________________________________________________
Double_t TsallisKernel(Double_t *x, Double_t *par);
Double_t TsallisFitNoInversion(Double_t *x, Double_t *par);


//Piece-wise Bose & mTExpo Fit_____________________________________
Double_t PiecewiseBosemTExpoFitKernel(Double_t *x, Double_t *par);
Double_t InvertPiecewiseBosemTExpoFit(Double_t *par, Double_t minRange=0.0, Double_t maxRange=10.0);
Double_t PiecewiseBosemTExpoFitNoInversion(Double_t *x, Double_t *par);
Double_t PiecewiseBosemTExpoFit(Double_t *x, Double_t *par);

//Piece-wise mTExpo & mTExpo Fit __________________________________
Double_t PiecewiseDoublemTExponentialKernel(Double_t *x, Double_t *par);
Double_t InvertPiecewiseDoublemTExponential(Double_t *par, Double_t minRange=0, Double_t maxRange=10.0);
Double_t PiecewiseDoublemTExponentialFitNoInversion(Double_t *x, Double_t *par);
Double_t PiecewiseDoublemTExponentialFit(Double_t *x, Double_t *par);

Double_t BoseBlastWaveFit(Double_t *x, Double_t *par);
