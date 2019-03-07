#ifndef PARTICLE_INFO_H
#define PARTICLE_INFO_H

#include <TString.h>

//Forward Declarations
class TF1;
//class TString;
class TGraph;
class TFile;

enum dEdxMethod {TRUNCATED, MOSTPROBABLE};
enum particleSpecies {PION, KAON, PROTON, ELECTRON, DEUTERON, TRITON, HELION, ALPHA, MUON, LAMBDA};

//ParticleInfo Class Declaration
class ParticleInfo : public TObject {

 private:
  Int_t type;
  TString libraryVersion;
  TFile *bichselRootFile;
  TFile *pidCalibrationFile;

  Double_t minBetaGamma;
  Double_t midBetaGamma;
  Double_t maxBetaGamma;

  Bool_t pidFunctionsCreated;
  Bool_t useEmpiricalFits;
  Int_t empiricalSpeciesInterest;
  Int_t empiricalRapidityIndex;

  Double_t canonicalTPCWidth;
  Double_t canonicalTOFWidth;
  
  static const int nDefinedParticles = 10;
  Double_t BichselStitch(Double_t *x, Double_t *par=NULL);
  Double_t BichselPrototype(Double_t *x,Double_t *par);
  Double_t InverseBetaPrototype(Double_t *x, Double_t *par);
  //TGraph *bichselGraph;
  TF1 *bichselFunction;
  TF1 *bichselFunctionPart[2];
  TF1 *bichselPIDFunctions[nDefinedParticles];
  TF1 *inverseBetaPIDFunctions[nDefinedParticles];

  TF1 *gausFunc; //Used in Weighted Z Calcs

  std::vector<std::vector<TF1 *> > bichselEmpiricalFit;
  std::vector<std::vector<std::vector<TF1 *> > > bichselEmpiricalPIDFunctions;
  
  Bool_t FindBichselFile(TString fileName, TString path, TString &foudnPath);
  void ReadBichselFile();

 public:
  ParticleInfo();
  ParticleInfo(TString starLibaryVersion, Bool_t createAll=false, const Int_t dEdxType=TRUNCATED);
  ParticleInfo(TString starLibaryVersion, Bool_t createAll, const int rapidityBins, TString pidCalibFile);
  virtual ~ParticleInfo();

  //Generic ParticleInfo
  Int_t GetNumberOfDefinedParticles(){return nDefinedParticles;}
  Double_t GetTPCWidth(){return canonicalTPCWidth;}
  Double_t GetTOFWidth(){return canonicalTOFWidth;}
  Double_t GetParticleMass(Int_t speciesIndex);
  Double_t GetParticleChargeMag(Int_t speciesIndex);
  TString GetParticleName(Int_t speciesIndex, Int_t charge=0);

  TString GetParticleSymbol(Int_t speciesIndex, Int_t charge=0);
  Style_t GetParticleMarker(Int_t speciesIndex, Int_t charge=0);
  Color_t GetParticleColor(Int_t speciesIndex);
  //TString GetParticleLine(Int_t speciesIndex, Int_t charge=0);
  //TString GetParticleSysColor(Int_t speciesIndex, Int_t charge=0);

  void CreateAllDefinedParticlePidFunctions();

  //Bichsel related functions
  Double_t BichselEmpiricalFitPrototype(Double_t *x, Double_t *par);
  Bool_t SetEmpiricalBichselFunction(Int_t speciesInterest, Int_t yIndex);
  TF1 *GetBichselFunction();
  TF1 *GetBichselFunction(Double_t mass, Double_t charge=1.0);
  TF1 *GetBichselFunction(Int_t speciesIndex);
  Double_t ComputeZTPC(Double_t dEdxMeasured, Double_t ptMeasured, Double_t pzMeasured, Int_t speciesIndex);
  Double_t ComputeWeightedZTPC(Double_t dEdxMeasured, Double_t ptMeasured, Double_t pzMeasured, Int_t speciesIndex,
			       Double_t avgTotalMomInBin); 
  Double_t PredictZTPC(Double_t totalMomentum, Int_t speciesIndexInterest, Int_t speciesIndexConfound);
  Double_t ConfoundMeanFitTPC(Double_t *x, Double_t *par);  
  
  //Inverse Beta related functions
  TF1 *GetInverseBetaFunction(Double_t mass, Double_t charge=1.0);
  TF1 *GetInverseBetaFunction(Int_t speciesIndex);
  Double_t ComputeZTOF(Double_t inverseBetaMeasured, Double_t ptMeasured, Double_t pzMeasured, Int_t speciesIndex);
  Double_t ComputeWeightedZTOF(Double_t inverseBetaMeasured, Double_t ptMeasured, Double_t pzMeasured, Int_t speciesIndex, Double_t avgTotalMomInBin);
  Double_t PredictZTOF(Double_t totalMomentum, Int_t speciesIndexInterest, Int_t speciesIndexConfound);
  Double_t ConfoundMeanFitTOF(Double_t *x, Double_t *par);

  //Other Useful Functions
  Int_t    GetNearestSpeciesTOF(Double_t inverseBetaMeasured, Double_t totalMomentum);
  Int_t    GetNearestSpeciesTPC(Double_t dEdxMeasured, Double_t totalMomentum);
  Int_t    GetNearestSpeciesTPCTOF(Double_t dEdxMeasured, Double_t inverseBetaMeasured, Double_t totalMomentum);
  Double_t GetPercentErrorTOF(Double_t inverseBetaMeasured, Double_t totalMomentum, Int_t speciesIndexInterest);
  Double_t GetPercentErrorTPC(Double_t dEdxMeasured, Double_t totalMomentum, Int_t speciesIndexInterest);
  Double_t GetDistanceToSpeciesTOF(Double_t inverseBetaMeasured, Double_t totalMomentum, Int_t speciesIndexInterest);
  Double_t GetDistanceToSpeciesTPC(Double_t dEdxMeasured, Double_t totalMomentum, Int_t speciesIndexInterest);  
  Bool_t   DoTPCandTOFAgreeOnNearestSpecies(Double_t dEdxMeasured, Double_t inverseBetaMeasured, Double_t totalMomentum);
  Bool_t   DoTPCandTOFAgreeOnThisSpecies(Double_t dEdxMeasured, Double_t inverseBetaMeasured, Double_t totalMomentum, Int_t speciesIndexInterest);
  Bool_t   IsExclusivelyThisSpeciesTPC(Double_t dEdxMeasured, Int_t speciesIndexInterest, Double_t totalMomentum,
				       Double_t percent=7.0, Bool_t electron=true);
  Bool_t   IsExclusivelyThisSpeciesTOF(Double_t inverseBetaMeasured, Int_t speciesIndexInterest, Double_t totalMomentum,
				       Double_t percent=5.0, Bool_t electron=true);
  Bool_t   IsExclusivelyThisSpeciesTPCTOF(Double_t dEdxMeasured, Double_t inverseBetaMeasured, Int_t speciesIndexInterest,
					  Double_t totalMomentum, Double_t tpcPercent=7.0,
					  Double_t tofPercent=5.0, Bool_t electron=true);
  
  ClassDef(ParticleInfo,1);

};

#endif
