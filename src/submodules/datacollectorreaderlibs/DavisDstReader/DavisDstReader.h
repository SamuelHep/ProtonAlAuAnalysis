#ifndef DAVISDSTREADER_H
#define DAVISDSTREADER_H

#ifdef private
#undef private
#endif

//Forward Declarations
class TrackInfo;
class PrimaryVertexInfo;
class EventInfo;
class TFile;
class TTree;
class TBranch;
class TString;
class TChain;
class TEntryList;
class TEventList;

class DavisDstReader : public TObject {

 private:
  TFile *file;
  TChain *tree;
  TEventList *eventList;
  TBranch *eventBranch;
  TBranch *vertexBranch;
  TBranch *trackBranch;

  Long64_t totalEntries;
  Long64_t currentEntry;
  
  EventInfo *eventInfo;
  TClonesArray *vertexArray;
  TClonesArray *trackArray;

  void ConstructTree(TString fileNameWithFullPath, Double_t minVz=-999, Double_t maxVz=-999);
  

 public:
  DavisDstReader();
  DavisDstReader(TString fileNameWithFullPath);
  DavisDstReader(TString fileNameWithFullPath, Double_t minVz, Double_t maxVz);
  virtual ~DavisDstReader();

  //Sets
  void SetBranchStatus(TString branchName,Int_t active);

  //Gets
  TTree *GetTree();
  Long64_t GetEntries(){return totalEntries;}
  EventInfo *GetEntry(Long64_t iEntry);
  EventInfo *FindEvent(Long64_t runNumber, Long64_t eventNumber);
  std::vector<long> GetEventsWithRunNumber(Long64_t runNumber);
  std::vector<long> GetEventsWithTriggerID(unsigned int triggerID);
  TClonesArray *GetPrimaryVertexArray(){return vertexArray;}

  ClassDef(DavisDstReader,1);
};


#endif
