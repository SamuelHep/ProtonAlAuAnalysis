/********************************************
             !!!!WARNING!!!!
This class is constructed assumming that the
spectra it holds are ploted as a function
of mT-m0 and with normalization factors
(1/2pi)x(1/mT)x(d^2N/dmTdy). All other types
of spectra will result in erroneous values!
*********************************************/

#include <TMatrixDSym.h>
class TGraphErrors;
class TMultiGraph;
class TDirectory;
class TF1;
class TString;
class ParticleInfo;
class TFile;

class SpectraClass : public TObject {

 private:
  int pid;
  int charge;
  int centIndex;
  int yIndex;
  int nPoints;
  TString name;
  unsigned int nSpectraFits;
  std::vector<double> binWidths;

  //Values that should be reported as being
  //the results for this spectrum
  //dNdy using the first fit integrated from [0,10]
  //quad sum of the two stat and sys errors below
  //stat error on dNdy of first fit
  //mean difference between dNdy's of all fits and first fit
  //mean pT using the first fit
  //mean difference between pTs of all fits and first fit
  //mean mT using the first fit
  //mean difference between mTs of all fits and first fit
  
  double dNdyResult;        
  double dNdyResultErr;     
  double dNdyResultStatErr; 
  double dNdyResultSysErrHigh;
  double dNdyResultSysErrLow;
  double meanpTResult;      
  double meanpTResultErr;   
  double meanmTResult;      
  double meanmTResultErr;

  //Other Quantities
  std::pair<double,double> dNdyFiducial;
  std::vector< std::pair<double,double> > dNdyLowExtrap;
  std::vector< std::pair<double,double> > dNdyHighExtrap;
  std::vector< std::pair<double,double> > dNdyVec;
  std::vector< std::pair<double,double> > meanpTVec;
  std::vector< std::pair<double,double> > meanmTVec;

  //These quantities are not saved to the root file
  TFile *spectraFile = NULL;      //! The File in which the following directories exist
  TString spectraDir;      //! Directory in which the spectra exist
  TString spectraClassDir; //! Directory in which the SpectraClass Object Exists
  TString spectraFitDir;   //! Directory in which the Spectra Fits Exist
  TString namePrefix;      //! Start of spectra name
  double *xVals = NULL; //!
  double *yVals = NULL; //!
  double *exVals = NULL; //!
  double *eyVals = NULL; //!
  TGraphErrors *tpcSpectrum = NULL; //!
  TGraphErrors *tofSpectrum = NULL; //!
  TMultiGraph *totalSpectrum = NULL; //!
  std::vector<TF1 *> spectraFits; //!
  std::vector<TMatrixDSym *>  spectraFitCovarianceMatrix; //!
  ParticleInfo particleInfo; //!

  void MakeMultiGraph();
  void ComputeInternalVectorQuantities(int which);
  void ComputedNdyHighExtrapolation(int which);
  void ComputedNdyLowExtrapolation(int which);

  Double_t scaleFunction(double *x, double *par);

  
 public:
  SpectraClass();
  SpectraClass(int pidIndex, int pidCharge, int centBin, int yBin,
	       TGraphErrors *tpc, TGraphErrors *tof=NULL);
  virtual ~SpectraClass();
  int GetN(){return nPoints;}
  int GetPid(){return pid;}
  int GetCharge(){return charge;}
  int GetRapidityIndex(){return yIndex;}
  int GetCentralityIndex(){return centIndex;}
  int GetNSpectrumFits(){return nSpectraFits;}
  TString GetSpectrumName(){return name;}
  double *GetX(){return xVals;}
  double *GetY(){return yVals;}
  double *GetEX(){return exVals;}
  double *GetEY(){return eyVals;}
  std::vector<double> *GetBinWidthVector(){return &binWidths;}
  TGraphErrors *GetTPCSpectrum(){return tpcSpectrum;}
  TGraphErrors *GetTOFSpectrum(){return tofSpectrum;}
  
  //Get Result Values
  double GetdNdyResult(){return dNdyResult;}
  double GetdNdyResultErr(){return dNdyResultErr;}
  double GetdNdyResultStatErr(){return dNdyResultStatErr;}
  double GetdNdyResultSysErrHigh(){return dNdyResultSysErrHigh;}
  double GetdNdyResultSysErrLow(){return dNdyResultSysErrLow;}
  double GetMeanPtResult(){return meanpTResult;}
  double GetMeanPtResultErr(){return meanpTResultErr;}
  double GetMeanmTResult(){return meanmTResult;}
  double GetMeanmTResultErr(){return meanmTResultErr;}

  //Get Values from specific fit
  TF1 *GetSpectrumFit(unsigned int which=0);
  double GetdNdy(unsigned int which=0);
  double GetdNdyErr(unsigned int which=0);
  double GetMeanmT(unsigned int which=0);
  double GetMeanmTErr(unsigned int which=0);
  double GetMeanpT(unsigned int which=0);
  double GetMeanpTErr(unsigned int which=0);
  double GetdNdyLowExtrapolation(unsigned int which=0);
  double GetdNdyLowExtrapolationErr(unsigned int which=0);
  double GetdNdyHighExtrapolation(unsigned int which=0);
  double GetdNdyHighExtrapolationErr(unsigned int which=0);

  //Get Value from Spectrum
  double GetdNdyFiducial();
  double GetdNdyFiducialErr();
  
  //Sets
  void SetSpectrumName(TString nameVal){name = nameVal;}
  void SetBinWidthVector(std::vector<double> binWidthVec){binWidths = binWidthVec;}
  void SetSpectraFile(TFile *file){spectraFile = file;}
  void SetSpectraDir(TString dir){spectraDir = dir;}
  void SetSpectraClassDir(TString dir){spectraClassDir = dir;}
  void SetSpectraNamePrefix(TString prefix){namePrefix = prefix;}
  void SetSpectraFitDir(TString dir){spectraFitDir = dir;}
  void SetdNdy(double val, unsigned int which=0);
  void SetdNdyStatErr(double val, unsigned int which=0);

  void SetdNdyResult(double val){dNdyResult = val;}
  void SetdNdyResultTotalErr(double val){dNdyResultErr = val;}
  void SetdNdyResultStatErr(double val){dNdyResultStatErr = val;}
  void SetdNdyResultSysErrHigh(double val){dNdyResultSysErrHigh = val;}
  void SetdNdyResultSysErrLow(double val){dNdyResultSysErrLow = val;}
  void SetMeanpTResult(double val){meanpTResult = val;}
  void SetMeanpTResultTotalErr(double val){meanpTResultErr = val;}
  void SetMeanmTResult(double val){meanmTResult = val;}
  void SetMeanmTResultTotalErr(double val){meanmTResultErr = val;}
  
  //Administrative Functions
  int LoadSpectra();
  void Fit(TF1 *fitFunc, TString opts="",TString gopts="", Axis_t minRange=0, Axis_t maxRange=0);
  void WriteSpectrum(bool overwriteSpectra=true, bool overwriteSpectraClass=true);
  void DrawSpectrum(TString opts="");
  void AddSpectraFit(TF1 *fit, TMatrixDSym *matrix=NULL/*TGraphErrors *confInterval=NULL*/, int which=-1);
  void AddSpectraFitCovMatrix(TMatrixDSym *matrix, int which=-1);
  void PrintSummary(int which =0,bool printPoints=false);

  //Functions which actually compute quantities.
  std::pair<double,double> ComputeTotaldNdy(int which=-1);
  std::pair<double,double> ComputedNdyInRange(int which=-1, double min=0.00001,double max=3.0);
  std::pair<double,double> ComputeFiducialdNdy(double min=-1, double max=-1);
  std::pair<double,double> ComputeMeanmTm0(int whichFunc=-1);
  //std::pair<double,double> ComputeTotaldNdy(TF1 *thisFunc, TGraphErrors *thisConf=NULL,
  //					    double min=0.0, double max=10.0);
  std::pair<double,double> ComputeTotaldNdy(TF1 *thisFunc, TMatrixDSym *thisCovMatrix=NULL,
					    double min=0.0, double max=10.0);
  //Functions Useful for Graphing
  TGraphErrors *GetTGraph();
  TMultiGraph *GetScaledCopy(double factor, int color=-1);
  TF1 *GetScaledFunction(double factor, unsigned int which=0, int color=-1);
  void SetMarkerColor(Color_t colorTPC,Color_t colorTOF=-1);
  void SetMarkerSize(Size_t sizeTPC, Size_t sizeTOF=-1);
  void SetMarkerStyle(Style_t styleTPC, Style_t styleTOF=-1);
  ClassDef(SpectraClass,3);
};
