#include <iostream>
#include <TObject.h>
#include <TClonesArray.h>
#include <TString.h>

#include "../TrackInfo/TrackInfo.h"
#include "PrimaryVertexInfo.h"

ClassImp(PrimaryVertexInfo);

//____________________________________________________________________________
PrimaryVertexInfo::PrimaryVertexInfo(){

  trackStartIndex    = -999;
  trackStopIndex     = -999;
  vertexIndexMuDst   = -999;
  nPrimaryTracks     = -999;
  ntofMatches        = -999;
  nPions             = -999;
  refMult            = -999;
  refMultUser        = -999;
  refMultPosY        = -999;
  refMultNegY        = -999;
  xVertex            = -999;
  yVertex            = -999;
  zVertex            = -999;

  trackArrayVertex = NULL;
  trackArrayEvent  = NULL;
}

//___________________________________________________________________________
PrimaryVertexInfo::~PrimaryVertexInfo(){

}

//____________________________________________________________________________
void PrimaryVertexInfo::ResetPrimaryVertexInfo(){
  
  trackStartIndex   = -999;
  trackStopIndex    = -999;
  vertexIndexMuDst  = -999;
  nPrimaryTracks    = -999;
  ntofMatches       = -999;
  nPions            = -999;
  refMult           = -999;
  refMultUser       = -999;
  refMultPosY       = -999;
  refMultNegY       = -999;
  xVertex           = -999;
  yVertex           = -999;
  zVertex           = -999;

  trackArrayVertex = NULL;
  trackArrayEvent  = NULL;
}

//____________________________________________________________________________
TClonesArray *PrimaryVertexInfo::GetPrimaryTrackArray(){

  //Return a pointer to an array of tracks belonging to vertexIndex

  //If trackArrayVertex exists then clear it
  if (trackArrayVertex)
    trackArrayVertex->Clear();
  else
    trackArrayVertex = new TClonesArray("TrackInfo",0);

  TClonesArray tempCopy(*trackArrayEvent);
  trackArrayVertex->AbsorbObjects(&tempCopy,trackStartIndex,trackStopIndex-1);

  return trackArrayVertex;

}

//____________________________________________________________________________
TrackInfo *PrimaryVertexInfo::GetPrimaryTrack(Int_t index){

  //Return a pointer to the primary track associated with this event at
  //index index

  //Make sure the Index is valid
  if (index < 0 || index >= nPrimaryTracks){
    fputs("ERROR: PrimaryVertexInfo::GetPrimaryTrack() - Requested track index is out of range!",stderr);
    exit (EXIT_FAILURE);
  }

  return (TrackInfo *)trackArrayEvent->At(trackStartIndex+index);

}

//____________________________________________________________________________
Float_t PrimaryVertexInfo::GetRefMultPosY(){

  //This needs to be recomputed (e.g. can not simply return from tree)
  //because the correct track quality cuts were not implemented at the
  //time of data collection.

  //The Track quality cuts for refMultPosY should be the same as
  //the ones used for refMult

  Float_t refMultPosY(0);
  TrackInfo *track;
  for (Int_t iTrack=0; iTrack<this->GetNPrimaryTracks(); iTrack++){

    track = this->GetPrimaryTrack(iTrack);

    if (track->GetEta() < 0.0 || track->GetEta() > 1.0)
      continue;
    
    if (track->GetFlag() < 0 || track->GetHitsFit() < 10)
      continue;

    if (track->GetPTotal() < 1.0e-10)
      continue;

    if (track->GetGlobalDCA()->Mag() > 3.0)
      continue;

    refMultPosY++;
  }

  return refMultPosY;
}

//____________________________________________________________________________
Float_t PrimaryVertexInfo::GetRefMultNegY(){

  //See comments in GetRefMultPosY
  
  Float_t refMultNegY(0);
  TrackInfo *track;
  for (Int_t iTrack=0; iTrack<this->GetNPrimaryTracks(); iTrack++){

    track = this->GetPrimaryTrack(iTrack);

    if (track->GetEta() < -1.0 || track->GetEta() > 0.0)
      continue;

    if (track->GetFlag() < 0 || track->GetHitsFit() < 10)
      continue;

    if (track->GetPTotal() < 1.0e-10)
      continue;

    if (track->GetGlobalDCA()->Mag() > 3.0)
      continue;

    refMultNegY++;
  }

  return refMultNegY;
}
  
//____________________________________________________________________________
void PrimaryVertexInfo::PrintPrimaryVertexInfo(Bool_t printTrackInfo,TString ind){

  cout <<"nPrimaryTracks: "   <<nPrimaryTracks  <<" "
       <<"nTofMatches: "      <<ntofMatches     <<" "
       <<Form("VertexLocation(X,Y,Z): (%g,%g,%g)\n",xVertex,yVertex,zVertex);

  if (!printTrackInfo)
    return;
  
  for (Int_t iTrack=0; iTrack<nPrimaryTracks; iTrack++){
    cout <<ind.Data();
    GetPrimaryTrack(iTrack)->PrintTrackInfo();
  }

}

