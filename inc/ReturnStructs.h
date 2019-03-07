//These are various structures used as return trypes for functions

class TH1D;
class TGraph2D;
class TCanvas;

//____________________________________________
struct NegBinomialSearchResults{

  //This Structure is returned from FindBestFitNegativeBinomialParameters()

  //These are variables used as part of FindBestFitNegativeBinomialParameters()
  Double_t npp;
  Double_t k;
  Double_t InverseChi2;
  TH1D *bestFitHisto;
  TGraph2D *chi2Graph;

  //These are variables used as part of FindCentralityBinCuts() and FindNpartNcollDistributions()
  Int_t nCentralityBins;
  std::vector<double> centralityBinDefinitions;
  std::vector<double> centralityBinCuts;
  TH1D *nPartTotalHisto;
  TH1D *nCollTotalHisto;
  TH1D *impactParamTotalHisto;
  TH1D *nPartHistos;
  TH1D *nCollHistos;
  TH1D *impactParamHistos;

  std::vector<double> nPartMeans;
  std::vector<double> nPartStatErrors;
  std::vector<double> nPartSysErrors;
  std::vector<double> nCollMeans;
  std::vector<double> nCollStatErrors;
  std::vector<double> nCollSysErrors;
  std::vector<double> impactParamMeans;
  std::vector<double> impactParamStatErrors;
  std::vector<double> impactParamSysErrors;

  //Canvas that could be useful for future drawing
  TCanvas *nPartTotalSysErrCanvas;
  TCanvas *nCollTotalSysErrCanvas;
  TCanvas *impactParamTotalSysErrCanvas;
  TCanvas *resultsTableCanvas;
};

//______________________________________________
//This class get written to a seperate "Correction Info File" that is designed to contain
//the minimal amount of things necessary to perform the correction operations during
//looping over data.
class RefMultCorrInfo {

 public:
  Long64_t triggerID;
  Double_t refMultCentBinCuts[16];

  RefMultCorrInfo();
  ~RefMultCorrInfo();

  ClassDef(RefMultCorrInfo,1);
};
