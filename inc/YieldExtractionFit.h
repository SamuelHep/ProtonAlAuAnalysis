class particleInfo;
class TF1;
class TH1D;

class YieldExtractionFit{

private:
  ParticleInfo *particleInfo;
  int minEntries;

  std::map<int,int> pidToFitParMap;

  //Mode TPC or TOF
  bool isTPCMode;

  //Bin Info 
  int pid, yIndex, mTm0Index;

  //Kinematic Info
  double mass, rapidity, mTm0, pT, pZ, mT, pTotal;

  //Yields
  double yield, yieldErr;

  //Yield Graphs
  TGraph *yieldGraph, *notYieldGraph;

  
  //Objects needed for the fits
  std::vector<TF1 *> tofOptFits;
  std::vector<TString> fitParNames;
  TH1D *dataHisto;

  float tofWidth;

  //Objects needed for yield extraction fit
  TF1 *fitFunc;
  TGraphErrors *confInterval;
  std::vector<TF1 *> singleDistributions;
  double minFitRange;
  double maxFitRange;
  TLine *minFitLine, *maxFitLine;

  //Useful Parameters for the zTPC Fits
  std::vector<double> predictedMean;
  std::vector<std::pair<double,bool> > finalFitAmpPars; //(val, isFixed)
  double electronAmp, electronMax;
  double slushParam;

  //Useful Parameters for the zTOF Fits
  std::vector<std::pair<double,bool> > tofFitPars; //(val,isFixed)

  double singleDistFunc(Double_t *x, Double_t *par);
  double singleDistFuncTOF(Double_t *x, Double_t *par);
  Double_t singleGaussianFunc(Double_t *x, Double_t *par);
  double yieldExtractFitTPC(Double_t *x, Double_t *par);
  double yieldExtractFitTOF(Double_t *x, Double_t *par);
  int FitTPC();
  int FitTOF();
  double ComputeTotalError(TGraphErrors *confInterval=NULL);


public:
  YieldExtractionFit();
  YieldExtractionFit(int particleSpecies, int yBin, int mTm0Bin,
		     TH1D *dHisto, ParticleInfo *pInfo, std::vector<TF1 *> tofOptVec);
  YieldExtractionFit(int particleSpecies, int yBin, int mTm0Bin,
		     TH1D *dHisto, ParticleInfo *pInfo, std::vector<TF1 *> tofFits, bool blah);
  YieldExtractionFit(int particleSpecies, int yBin, int mTm0Bin,
		     TH1D *dHisto, ParticleInfo *pInfo);
  ~YieldExtractionFit();

  int Fit();

  void DrawIndividualDistributions();
  void DrawFit(TString opt);
  std::pair<double,double> ComputeYield(double nEvents,
					TGraphErrors *altConfInterval=NULL, TF1 *altFitFunc = NULL);
  void SetElectronAmpMax(double eAmpMax){electronMax = eAmpMax;}
  void FixAmpParameter(int particle, double val);
  void FixHorizontalSlushParameter(double val){slushParam = val;}

  TF1 *GetFit(){return fitFunc;}
  std::vector<TF1 *> *GetIndividualFits(){return &singleDistributions;}
  TGraphErrors *GetFitConfidenceInterval(){return confInterval;}
  TGraph *GetYieldGraph(){return yieldGraph;}
  TGraph *GetNotYieldGraph(){return notYieldGraph;}
  double GetMinFitRange(){return minFitRange;}
  double GetMaxFitRange(){return maxFitRange;}
  double GetmTm0(){return mTm0;}
  double GetRapidity(){return rapidity;}
  double GetmT(){return mT;}
  double GetAmpParameter(int particle);
  double GetSlushParameter();

  void SetTofFitParameters(std::vector<std::pair<double,bool> > vals){tofFitPars = vals;}
  std::vector<std::pair<double,bool> > GetTofFitParameters(){return tofFitPars;}
};
