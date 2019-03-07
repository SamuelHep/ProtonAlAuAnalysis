#ifndef STREFMULTEXTENDEDCORR_H
#define STREFMULTEXTENDEDCORR_H


//This class is responsible for applying the refMult, 
//refMultPosY, and refMultNegY Correction and then determining
//the centrality bin of the event. 

//To work the class needs to know the event's:
//triggerID, refMult(PosY/NegY), and zVertex

class TFile;
class TTree;
class TGraphErrors;
class TF1;
#include <TVectorD.h>


class RefMultCorrInfo{

 public:
  RefMultCorrInfo();
  virtual ~RefMultCorrInfo();

  //public:
  Long64_t triggerID;
  Double_t refMultCentBinCuts[16];

  ClassDef(RefMultCorrInfo,2);
};


class StRefMultExtendedCorr : public TObject{


 private:
  Bool_t isInitialized;
  double cmEnergy;
  TFile *corrFile;
  TTree *corrTree;
  RefMultCorrInfo *treeEntry;

  int nTriggers;
  TGraphErrors ***refMultCorrGraph;
  TF1 ***convertToRefMultFunc;
  TVectorD *badRunList;

  int GetTriggerIndex(long triggerID);
  int GetRefMultType(double zVertex);
  Bool_t CanProcessEvent(long triggerID, double zVertex);
  Bool_t FindCorrectionFile(TString fileName, TString path, TString &foundPath);


 public:
  StRefMultExtendedCorr();
  virtual ~StRefMultExtendedCorr();

  void Initialize(double energy);
  int GetNTriggers(){return nTriggers;}
  double GetRefMultCorrection(long triggerID, double zVertex, int refMultType=-1);
  double GetCorrectedRefMult(long triggerID, double zVertex, double refMult, int refMultType=-1);
  double ConvertToRefMult(long triggerID, double zVertex, double refMultPosNegY);
  int GetCentralityBin16(long triggerID, double zVertex, double refMult, int refMultType=-1);
  int GetCentralityBin9(long triggerID, double zVertex, double refMult, int refMultType=-1);
  Bool_t IsBadRun(long runNumber);
 
  ClassDef(StRefMultExtendedCorr,1);

};

#endif
