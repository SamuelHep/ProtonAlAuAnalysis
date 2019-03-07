#ifndef SIMFITTER_H
#define SIMFITTER_H

#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

struct GlobalChi2;


//_____________________________________________________________________
class SimFitter : public TObject{
  
private:
  bool isInitialized;
  unsigned int nObjects;
  std::vector<TGraphErrors *> objectsToFit;
  ROOT::Fit::Fitter fitter;
  std::vector<double> totalFitPars;
  unsigned int dataSize;
  int printLevel;
  ROOT::Fit::FitResult fitResult;

  std::vector<TF1 *> fitFuncResults;
  std::vector<ROOT::Math::WrappedMultiTF1 *> fitFuncs;
  std::vector<std::vector<double> >          fitPars;
  std::vector<ROOT::Fit::DataOptions *>      dataOption;
  std::vector<ROOT::Fit::DataRange *>        dataRange;
  std::vector<ROOT::Fit::BinData *>          binData;
  std::vector<ROOT::Fit::Chi2Function *>     chi2Func;
  std::vector<std::map<int,int> >            parameterMap; //From iFunc (iFuncParN,totalFitParN)
  std::vector<std::vector<int> >             parameterRelationships;
  std::vector<TMatrixDSym *>                 covarianceMatrix;
  
  GlobalChi2 *globalChi2;
  void ConstructCovarianceMatrix();
  TMatrixD GetRandomDeviateVector(int iFit);
  
public:
  SimFitter(std::vector<TGraphErrors *> objs, std::vector<TF1 *> fitFuncs);
  ~SimFitter();

  void Init(std::vector<double> parVals,
	    std::vector<std::pair<float,float> > *range = NULL);
  void FixParameter(int iPar){fitter.Config().ParSettings(iPar).Fix();}
  void SetParLimits(int iPar, double valLow, double valHigh){fitter.Config().ParSettings(iPar).SetLimits(valLow,valHigh);}
  void SetPrintLevel(int val){printLevel = val;}
  void SetParameterRelationships(std::vector<std::vector<int> > parRelations);
  Double_t ComputeGlobalChi2(const Double_t *pars);
  ROOT::Fit::FitResult Fit();
  std::vector<TF1 *> GetFitFunctions(){return fitFuncResults;}
  std::vector<TGraphErrors *> GetConfidenceIntervals(float cl=.638);
  TGraphErrors *GetConfidenceInterval(int iFit, float cl=.683);
  ROOT::Fit::FitResult GetFitResult(){return fitResult;}
  TMatrixDSym *GetCovarianceMatrix(int iFit){return covarianceMatrix.at(iFit);}
  
  
  ClassDef(SimFitter,1);
};

//********************************************************************
struct GlobalChi2 {

  SimFitter *fitter = NULL;

  GlobalChi2 (SimFitter *simFit){fitter = simFit;}

  double operator() (const double *par) const {
    return fitter->ComputeGlobalChi2(par);
  }

};

#endif
