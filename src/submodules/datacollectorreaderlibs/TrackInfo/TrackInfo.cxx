#include <iostream>
#include "TObject.h"
#include "TVector3.h"
#include "TString.h"

#include "TrackInfo.h"

ClassImp(TrackInfo);

//__________________________________________________________
TrackInfo::TrackInfo(){

  nHits = -999;
  nHitsFit = -999;
  nHitsPoss = -999;
  dEdxHits = -999;
  flag = -999;
  pZ = -999;
  pT = -999;
  charge = -999;
  dEdxTruncated = -999;
  dEdxFit = -999;
  tof = -999;
  tofBeta = -999;
  eta = -999;
  phi =-999;
  nSigmaElectron = -999;
  nSigmaPion = -999;
  nSigmaKaon = -999;
  nSigmaProton = -999;
  length = -999;
  lengthMeasured = -999;
  chi2 = -999;
  nSigmaElectronTOF = -999;
  nSigmaPionTOF = -999;
  nSigmaKaonTOF = -999;
  nSigmaProtonTOF = -999;
  yLocalTOF = -999;
  zLocalTOF = -999;
  globalDCA.SetXYZ(-999,-999,-999);
  firstPoint.SetXYZ(-999,-999,-999);
  lastPoint.SetXYZ(-999,-999,-999);
  tofMatchFlag = 0;
}

//_________________________________________________________
TrackInfo::~TrackInfo(){

  //Nothing to do here!

}

//_________________________________________________________
Float_t TrackInfo::GetNSigmaOfSpecies(Int_t particleSpecies){

  //Returns the NSigma of track for particleSpecies
  //NOTE: This uses the same enumeration as ParticleInfo,
  //but be careful because there is not an nSigma value for
  //all particle spcies!

  if (particleSpecies == 0)
    return nSigmaPion;
  else if (particleSpecies == 1)
    return nSigmaKaon;
  else if (particleSpecies == 2)
    return nSigmaProton;
  else if (particleSpecies == 3)
    return nSigmaElectron;
  else 
    cout <<"WARNING - TrackInfo::GetNSigmaOfSpecies(): Species does not have an NSigma "
	 <<"value to return. RETURNING -999!\n";

  return -999;
}

//_________________________________________________________
Float_t TrackInfo::GetTofNSigmaOfSpecies(Int_t particleSpecies){

  //Returns the TofNSigma of track for particleSpecies
  //NOTE: This uses the same enumeration as ParticleInfo,
  //but be careful because there is not an a TofNSigma value
  //for all particle species!

  if (particleSpecies == 0)
    return nSigmaPionTOF;
  else if (particleSpecies == 1)
    return nSigmaKaonTOF;
  else if (particleSpecies == 2)
    return nSigmaProtonTOF;
  else if (particleSpecies == 3)
    return nSigmaElectronTOF;
  else
    cout <<"WARNING - TrackInfo::GetTofNSigmaOfSpecies(): Species does not have an TofNSigma "
	 <<"value to return. RETURNING -999!\n";

  return -999;
}

//_________________________________________________________
Float_t TrackInfo::GetRapidity(Float_t massAssumption){

  Float_t energy = sqrt( pow(pT,2) + pow(pZ,2) + pow(massAssumption,2) );

  return 0.5 * log( (energy + pZ) / (energy - pZ) );

}

//_________________________________________________________
Float_t TrackInfo::GetmTm0(Float_t massAssumption){
  
  return sqrt( pow(pT,2) + pow(massAssumption,2) ) - massAssumption;

}

//_________________________________________________________
Float_t TrackInfo::GetBetaGamma(Float_t massAssumption){

  return GetPTotal() / massAssumption;
}

//_________________________________________________________
Float_t TrackInfo::GetNFrac(){

  //Return the fraction of hits used in the fit
  return (Float_t)nHitsFit / (Float_t)nHitsPoss;

}

//_________________________________________________________
Float_t TrackInfo::GetPTotal(){

  return sqrt( pow(pT, 2) + pow(pZ, 2) );
  
}

//_________________________________________________________
void TrackInfo::PrintTrackInfo(){

  cout <<"ParentVertex: "        <<parentVertex <<" "
       <<Form("Momentum (pt,pz): (%g,%g)\n",pT,pZ);

}
