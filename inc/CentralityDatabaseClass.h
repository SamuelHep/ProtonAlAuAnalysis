#ifndef CENTRALITYDATABASECLASS_H
#define CENTRALITYDATABASECLASS_H

//This is a simple class to hold basic run/event info
//so that it can be stored to classify the centrality
//of embedded events.

#include <vector>

class TClonesArray;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CentralityDatabase : public TObject {

 private:
  Long64_t runNumber;
  Double_t minVz;
  Double_t maxVz;
  Double_t minVr;
  Double_t maxVr;
  Double_t beamSpotX;
  Double_t beamSpotY;
  
  TClonesArray *eventArray;

 public:
  CentralityDatabase();
  CentralityDatabase(Long64_t rNumber);
  ~CentralityDatabase();

  void SetRunInfo(Long64_t rNumber);
  void SetVzRange(Double_t valMin, Double_t valMax){ minVz=valMin; maxVz=valMax; }
  void SetVrRange(Double_t valMin, Double_t valMax){ minVr=valMin; maxVr=valMax; }
  void SetBeamSpot(Double_t valX, Double_t valY){ beamSpotX=valX, beamSpotY=valY; }
  
  void AddEvent(Long64_t eNumber, Int_t centBin, Double_t xVertex, Double_t yVertex, Double_t zVertex);
  Int_t GetCentralityBin(Long64_t eNumber);
  Long64_t GetEventNumber(Int_t iEvent);
  Long64_t GetRunNumber(){ return runNumber; }
  Int_t GetNEvents();
  Double_t GetMinVz(){ return minVz; }
  Double_t GetMaxVz(){ return maxVz; }
  Double_t GetMinVr(){ return minVr; }
  Double_t GetMaxVr(){ return maxVr; }
  Double_t GetBeamSpotX(){ return beamSpotX; }
  Double_t GetBeamSpotY(){ return beamSpotY; }
  
  ClassDef(CentralityDatabase,1);

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CentralityDatabaseEventInfo : public TObject {

 private:
  Long64_t eventNumber;
  Int_t centralityIndex;
  Double_t xVertex;
  Double_t yVertex;
  Double_t zVertex;

 public:
  CentralityDatabaseEventInfo();
  virtual ~CentralityDatabaseEventInfo();

  void SetEventInfo(Long64_t eNumber,Int_t centBin);
  void SetEventVertex(Double_t xVal, Double_t yVal, Double_t zVal){ xVertex=xVal; yVertex=yVal; zVertex=zVal; };
  Long64_t GetEventNumber(){ return eventNumber; }
  Int_t GetCentralityBin(){ return centralityIndex; }
  Double_t GetXVertex(){ return xVertex; }
  Double_t GetYVertex(){ return yVertex; }
  Double_t GetZVertex(){ return zVertex; }

  
  ClassDef(CentralityDatabaseEventInfo,1);

};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Other Useful functions
Bool_t IsNewRunNumber(std::vector<CentralityDatabase *> *vec,
		      Long64_t rNumber);
CentralityDatabase *GetRunNumberEntry(std::vector<CentralityDatabase *> *vec,
					Long64_t rNumber);
#endif
