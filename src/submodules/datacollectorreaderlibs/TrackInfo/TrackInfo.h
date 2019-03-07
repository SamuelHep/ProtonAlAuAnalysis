#ifndef TRACK_INFO_H
#define TRACK_INFO_H

#include "TVector3.h"

#ifdef private
#undef private
#endif

//Forward Declarations
class TVector3;

//StTrackInfo Class Definition
class TrackInfo : public TObject {

 private:
  //NOTE: All variables are for primary tracks
  //      unless prefixed with global.

  //WARNING: The list of variables must be in exactly the
  //         same order and have the same types as those in
  //	     StDataCollector/TrackInfo.h

  //Data Members
  Int_t parentVertex;
  Int_t nHits;
  Int_t nHitsFit;
  Int_t nHitsPoss;
  Int_t dEdxHits;
  Int_t charge;
  short flag;
  Float_t pZ;
  Float_t pT;
  Float_t dEdxTruncated;
  Float_t dEdxFit;
  Float_t tof;
  Float_t tofBeta;
  Float_t eta;
  Float_t phi;
  Float_t nSigmaElectron;
  Float_t nSigmaPion;
  Float_t nSigmaKaon;
  Float_t nSigmaProton;
  Float_t length;
  Float_t lengthMeasured;
  Float_t chi2;
  Float_t nSigmaElectronTOF;
  Float_t nSigmaPionTOF;
  Float_t nSigmaKaonTOF;
  Float_t nSigmaProtonTOF;
  Float_t yLocalTOF;
  Float_t zLocalTOF;
  Float_t pathLengthTOF;

  TVector3 globalDCA;
  TVector3 firstPoint;
  TVector3 lastPoint;

  unsigned char tofMatchFlag;


 public:
  TrackInfo();
  virtual ~TrackInfo();

  //Member Functions
  Int_t GetParentVertex(){return parentVertex;}
  Int_t GetHits(){return nHits;}
  Int_t GetHitsFit(){return nHitsFit;}
  Int_t GetHitsPoss(){return nHitsPoss;}
  Int_t GetdEdxHits(){return dEdxHits;}
  Int_t GetCharge(){return charge;}
  short GetFlag(){return flag;}
  Float_t GetPz(){return pZ;}
  Float_t GetPt(){return pT;}
  Float_t GetdEdx(Int_t method=0){ return method==0 ? dEdxTruncated : dEdxFit;}
  Float_t GetTof(){return tof;}
  Float_t GetTofBeta(){return tofBeta;}
  Float_t GetEta(){return eta;}
  Float_t GetPhi(){return phi;}
  Float_t GetNSigmaElectron(){return nSigmaElectron;}
  Float_t GetNSigmaPion(){return nSigmaPion;}
  Float_t GetNSigmaKaon(){return nSigmaKaon;}
  Float_t GetNSigmaProton(){return nSigmaProton;}
  Float_t GetLength(){return length;}
  Float_t GetMeasuredLength(){return lengthMeasured;}
  Float_t GetChi2(){return chi2;}
  Float_t GetTofNSigmaElectron(){return nSigmaElectronTOF;}
  Float_t GetTofNSigmaPion(){return nSigmaPionTOF;}
  Float_t GetTofNSigmaKaon(){return nSigmaKaonTOF;}
  Float_t GetTofNSigmaProton(){return nSigmaProtonTOF;}
  Float_t GetTofYLocal(){return yLocalTOF;}
  Float_t GetTofZLocal(){return zLocalTOF;}
  Float_t GetTofPathLength(){return pathLengthTOF;}
  Int_t   GetTofMatchFlag(){return (Int_t)tofMatchFlag;}

  TVector3 *GetGlobalDCA(){return &globalDCA;}
  TVector3 *GetFirstPoint(){return &firstPoint;}
  TVector3 *GetLastPoint(){return &lastPoint;}

  //Member Functions that compute quantities
  Float_t GetNSigmaOfSpecies(Int_t particleSpecies);
  Float_t GetTofNSigmaOfSpecies(Int_t particleSpecies);
  Float_t GetRapidity(Float_t massAssumption);
  Float_t GetmTm0(Float_t massAssumption);
  Float_t GetBetaGamma(Float_t massAssumption);
  Float_t GetPTotal();
  Float_t GetNFrac();

  //Print
  void PrintTrackInfo();

  ClassDef(TrackInfo,1);

};

#endif
