#ifndef PRIMARYVERTEX_INFO_H
#define PRIMARYVERTEX_INFO_H

#ifdef private
#undef private
#endif

//Forward Declarations
class TClonesArray;
class TrackInfo;

//PrimaryVertexInfo Class Definition
class PrimaryVertexInfo : public TObject {

  //WARNING: The data members here must be in the same
  //         order and have the same type as in 
  //         datacollector/primaryvertexinfo.

 private:
  Int_t trackStartIndex;     //Index in the track branch
  Int_t trackStopIndex;      //Index in the track branch
  Int_t vertexIndexMuDst;    //Index of vertex in the MuDst
  Int_t nPrimaryTracks;      //Number of Primary Tracks
  Int_t ntofMatches;         //Number of Primary Tracks with TOF Matches
  Int_t nPions;              //Rough estimate of number of pions
  Float_t refMult;           //Standard RefMult |dN/dEta| < .5
  Float_t refMultUser;       //RefMult in User defined Range
  Float_t refMultPosY;       //RefMult in [0,1)
  Float_t refMultNegY;       //RefMult in [-1,0)
  Float_t xVertex;           //X vertex location (cm)
  Float_t yVertex;           //Y vertex location (cm)
  Float_t zVertex;           //Z vertex location (cm)

  TClonesArray *trackArrayEvent;  //!
  TClonesArray *trackArrayVertex; //!

 public:
  
  //Constructors and Destructors
  PrimaryVertexInfo();
  ~PrimaryVertexInfo();
  
  //Member Functions to get the stored values
  Int_t GetVertexIndex(){return vertexIndexMuDst;}
  Int_t GetNPrimaryTracks(){return nPrimaryTracks;}
  Int_t GetNTofMatches(){return ntofMatches;}
  Int_t GetNPions(){return nPions;}
  Float_t GetRefMult(){return refMult;}
  Float_t GetRefMultUser(){return refMultUser;}
  Float_t GetRefMultPosY();
  Float_t GetRefMultNegY();
  Float_t GetXVertex(){return xVertex;}
  Float_t GetYVertex(){return yVertex;}
  Float_t GetZVertex(){return zVertex;}

  //Member Functions which do something
  void SetEventTrackArrayPtr(TClonesArray *arr){trackArrayEvent = arr;}
  TClonesArray *GetPrimaryTrackArray();
  TrackInfo *GetPrimaryTrack(Int_t index);
  Float_t GetRVertex(){return sqrt(pow(xVertex,2)+pow(yVertex,2));}
  void ResetPrimaryVertexInfo();
  void PrintPrimaryVertexInfo(Bool_t printTrackInfo = false,TString ind="");

  ClassDef(PrimaryVertexInfo,1);

};

#endif
