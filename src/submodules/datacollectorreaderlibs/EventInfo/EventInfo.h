#ifndef EVENT_INFO_H
#define EVENT_INFO_H

#ifdef private
#undef private
#endif

//Forward Declarations
class TClonesArray;
class PrimaryVertexInfo;

//EventInfo Class Definition
class EventInfo : public TObject {

 private:

  //WARNING: The data members must be in the same order
  //         and have the same type as the data members
  //         in stdatacollector/eventinfo

  //Data Members
  Int_t nPrimaryVertices;  //Number of stored in davisdst
  Int_t runNumber;         //Run Number
  Int_t eventNumber;       //Event Number
  Int_t nTotalVertices;    //Number of primary vertices reported in mudst
  Int_t nPileUpVertices;   //Number of pile up vertices reported in mudst
  Int_t adcSumBBCWest;     //Sum of ADC values from BBC West
  Int_t adcSumBBCEast;     //Sum of ADC values from BBC East
  Int_t tofMultiplicity;   //Total number of ToF modules with a hit
  Float_t meanPt;          //Average Pt for all global tracks in this event
  Float_t meanEta;         //Average Eta for all global tracks in this event
  std::vector <unsigned int> triggerIDs;

  PrimaryVertexInfo *vertex; //!
  TClonesArray *primaryVertexArray; //!
  TClonesArray *primaryTrackArray;  //!

 public:

  //Constructors and Destructors
  EventInfo();
  ~EventInfo();

  //Member Functions to get the stored quantities
  Int_t GetNPrimaryVertices(){return nPrimaryVertices;}
  Int_t GetRunNumber(){return runNumber;}
  Int_t GetEventNumber(){return eventNumber;}
  Int_t GetNTotalVertices(){return nTotalVertices;}
  Int_t GetNPileUpVertices(){return nPileUpVertices;}
  Int_t GetADCSumBBCWest(){return adcSumBBCWest;}
  Int_t GetADCSumBBCEast(){return adcSumBBCEast;}
  Float_t GetMeanPt(){return meanPt;}
  Float_t GetMeanEta(){return meanEta;}
  std::vector <unsigned int> *GetTriggerIDVector(){return &triggerIDs;}
  Int_t GetTofMultiplicity(){return tofMultiplicity;}
  Bool_t ContainsTriggerID(unsigned int trigID);
  UInt_t GetTriggerID(std::vector <unsigned int> *allowedTriggers);

  //Sets
  void SetVertexArrayPtr(TClonesArray *ptr);
  void SetTrackArrayPtr(TClonesArray *ptr);

  //Member Functions to get other objects
  PrimaryVertexInfo *GetPrimaryVertex(Int_t index);

  void ResetEventInfo();
  void PrintEventInfo(Bool_t printVertexInfo=false, Bool_t printTrackInfo=false, Int_t wichVertex=-1);

  ClassDef(EventInfo,1);

};

#endif
