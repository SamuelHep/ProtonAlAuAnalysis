#include <iostream>
#include <TObject.h>
#include "TClonesArray.h"

#include "../PrimaryVertexInfo/PrimaryVertexInfo.h"
#include "EventInfo.h"

ClassImp(EventInfo);

//___________________________________________________________
EventInfo::EventInfo(){

  nPrimaryVertices   = -999;
  runNumber          = -999;
  eventNumber        = -999;
  nTotalVertices     = -999;
  nPileUpVertices    = -999;
  adcSumBBCWest      = -999;
  adcSumBBCEast      = -999;
  meanPt             = -999;
  meanEta            = -999;
  tofMultiplicity    = -999;
  triggerIDs.clear();

  primaryVertexArray = NULL;
  primaryTrackArray  = NULL;

}

//__________________________________________________________
EventInfo::~EventInfo(){

}

//__________________________________________________________
void EventInfo::ResetEventInfo(){

  nPrimaryVertices    = -999;
  runNumber          = -999;
  eventNumber        = -999;
  nTotalVertices     = -999;
  nPileUpVertices    = -999;
  adcSumBBCWest      = -999;
  adcSumBBCEast      = -999;
  meanPt             = -999;
  meanEta            = -999;
  tofMultiplicity    = -999;

  primaryVertexArray = NULL;
  primaryTrackArray  = NULL;
}

//__________________________________________________________
void EventInfo::PrintEventInfo(Bool_t printVertexInfo, Bool_t printTrackInfo, Int_t whichVertex){

  //Print information about the event.
  //To print information about the event's vertices set printVertexInfo to true
  //To print information about the event's vertices and tracks also set printTrackInfo to true
  //If you would like to choose a particular vertex to print use the whichVertex Variable

  cout <<Form("RunNumber: %d  EventNumber: %d PrimaryVertices: %d \n",runNumber, eventNumber,nPrimaryVertices);

  if (printVertexInfo){
    
    if (whichVertex < 0){
      for (Int_t iVertex=0; iVertex < nPrimaryVertices; iVertex++){
	cout <<"    ";
	GetPrimaryVertex(iVertex)->PrintPrimaryVertexInfo(printTrackInfo,"        ");
      }
    }
    else
      GetPrimaryVertex(whichVertex)->PrintPrimaryVertexInfo(printTrackInfo);
  }

}

//__________________________________________________________
void EventInfo::SetVertexArrayPtr(TClonesArray *ptr){

  //Pass in the pointer to the vertex array
  if (!ptr){
    fputs("ERROR: EventInfo::SetVertexArrayPtr() - The value of the pointer is NULL!",stderr);
    exit (EXIT_FAILURE);
  }

  primaryVertexArray = ptr;
}

//__________________________________________________________
void EventInfo::SetTrackArrayPtr(TClonesArray *ptr){

  //Pass in the pointer to the Events Track Array
  if (!ptr){
    fputs("ERROR: EventInfo::SetTrackArrayPtr() - The value of the pointer is NULL!",stderr);
    exit (EXIT_FAILURE);
  }

  primaryTrackArray = ptr;
}


//__________________________________________________________
PrimaryVertexInfo *EventInfo::GetPrimaryVertex(Int_t index){

  //Return a Pointer to the primary vertex at index index

  //Check to make sure the index is valid
  if (index < 0){
    fputs("ERROR: EventInfo::GetPrimaryVertex() - Vertex index must be positive!",stderr);
    exit (EXIT_FAILURE);
  }
  if (index >= primaryVertexArray->GetEntries()){
    fputs("WARNING: EventInfo::GetPrimaryVertex() - Vertex index requested is out of range! Returning NULL.",stderr);
    return NULL;
  }

  vertex = (PrimaryVertexInfo *)primaryVertexArray->At(index);
  vertex->SetEventTrackArrayPtr(primaryTrackArray);

  return vertex;

}

//___________________________________________________________
Bool_t EventInfo::ContainsTriggerID(unsigned int trigID){

  //Returns true if the trigger ID vector contains an entry
  //that is equal to the argument trigID

  for (unsigned int i=0; i<triggerIDs.size(); i++){
    if (trigID == triggerIDs.at(i))
      return true;
  }

  //If the function makes it this far then all the entries
  //have been checked and none match trigID
  return false;
  
}

//___________________________________________________________
UInt_t EventInfo::GetTriggerID(std::vector <unsigned int> *allowedTriggers){

  //Return the trigger ID of this event from the list of allowed triggers
  //If the event does not have one of the allowed tiggers then
  //return 999

  for (unsigned int i=0; i<allowedTriggers->size(); i++){

    if (this->ContainsTriggerID(allowedTriggers->at(i)))
      return allowedTriggers->at(i);
    
  }
  
  return 999;

}
