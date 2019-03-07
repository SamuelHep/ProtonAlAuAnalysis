#include <iostream>
#include <vector>
#include <TObject.h>
#include <TClonesArray.h>

#include "CentralityDatabaseClass.h"

using namespace std;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ClassImp(CentralityDatabase);

//_____________________________________________________________
CentralityDatabase::CentralityDatabase(){
  runNumber = -999;

  eventArray = new TClonesArray("CentralityDatabaseEventInfo",0);
}

//_____________________________________________________________
CentralityDatabase::CentralityDatabase(Long64_t rNumber){

  runNumber = rNumber;

  eventArray = new TClonesArray("CentralityDatabaseEventInfo",0);

}

//_____________________________________________________________
CentralityDatabase::~CentralityDatabase(){
  //Destructor
  delete eventArray;
}

//_____________________________________________________________
void CentralityDatabase::SetRunInfo(Long64_t rNumber){

  runNumber = rNumber;

}

//_____________________________________________________________
Int_t CentralityDatabase::GetNEvents(){
  return eventArray->GetEntriesFast();
}

//_____________________________________________________________
Long64_t CentralityDatabase::GetEventNumber(Int_t iEvent){

  return ((CentralityDatabaseEventInfo *)eventArray->At(iEvent))->GetEventNumber();
}

//_____________________________________________________________
void CentralityDatabase::AddEvent(Long64_t eNumber, Int_t centBin, Double_t xVertex,
				  Double_t yVertex, Double_t zVertex){
  
  Int_t back = eventArray->GetLast() + 1;
  CentralityDatabaseEventInfo *tempEvent = (CentralityDatabaseEventInfo *)eventArray->ConstructedAt(back);
  tempEvent->SetEventInfo(eNumber,centBin);
  tempEvent->SetEventVertex(xVertex,yVertex,zVertex);
}

//______________________________________________________________
Int_t CentralityDatabase::GetCentralityBin(Long64_t eNumber){

  //Loop Over the Elements of the Event Array and 
  //return the centrality Bin of the event matching eNumber
  for (Int_t iEvent=0; iEvent<eventArray->GetEntriesFast(); iEvent++){

    CentralityDatabaseEventInfo *event = (CentralityDatabaseEventInfo *)eventArray->At(iEvent);

    if ( event->GetEventNumber() != eNumber )
      continue;

    return event->GetCentralityBin();

  }

  //Failure Mode
  return -1;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ClassImp(CentralityDatabaseEventInfo);

//______________________________________________________________
CentralityDatabaseEventInfo::CentralityDatabaseEventInfo(){

  eventNumber = -999;
  centralityIndex = -1;

}

//______________________________________________________________
CentralityDatabaseEventInfo::~CentralityDatabaseEventInfo(){
  //Destructor
}

//______________________________________________________________
void CentralityDatabaseEventInfo::SetEventInfo(Long64_t eNumber,Int_t centBin){

  eventNumber = eNumber;
  centralityIndex = centBin;

}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Some Functions for using these classes

//_______________________________________________________________
Bool_t IsNewRunNumber(std::vector<CentralityDatabase *> *vec, 
		      Long64_t rNumber){

  //Loops over the vector vec containing pointers to 
  //CentralityDatabase objects and compared rNumber to the
  //runNumbers stored in the objects of vec.  If rNumber
  //matches a runNumber of the object the result is false,
  //otherwise the result is true - e.g. rNumber is a new runNumber

  for (unsigned int i=0; i<vec->size(); i++){

    if (rNumber == vec->at(i)->GetRunNumber())
      return false;

  }

  return true;
}

//________________________________________________________________
CentralityDatabase *GetRunNumberEntry(std::vector<CentralityDatabase *> *vec,
					Long64_t rNumber){

  //Loop Over the entries of vec and return a pointer to the one
  //that has a runNumber that Matches rNumber
  
  for (unsigned int i=0; i<vec->size(); i++){

    if (rNumber == vec->at(i)->GetRunNumber())
      return vec->at(i);

  }

  return NULL;
}
