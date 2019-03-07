//A header file for the user accessible cuts

#include "StRefMultExtendedCorr.h"
#include "TrackInfo.h"
#include "PrimaryVertexInfo.h"
#include "EventInfo.h"


void SetVariableUserCuts(Double_t energy=0, TString eventConfig="",TString starLibrary="");
void CheckCutStatus();

void SetStarLibraryVersion(TString starLibVer="");
TString GetStarLibraryVersion();

void SetEventConfiguration(TString eventConfig="");
TString GetEventConfiguration();

void SetCollisionEnergy(Double_t energy = 0.0);
Double_t GetCollisionEnergy();

Bool_t IsGoodEvent(EventInfo *event);
Bool_t IsGoodVertex(PrimaryVertexInfo *vertex);
Bool_t IsGoodTrack(TrackInfo *track, Bool_t useThisCut=true);
Bool_t IsGoodTofTrack(TrackInfo *track);

void SetCentralityPercents(Double_t *centPercentArr=NULL, const int nCentBins=-1);
std::vector<double> GetCentralityPercents();

void SetCentralityCuts(Int_t *centCutsArr=NULL, const int nCentBins=-1);
std::vector<int> GetCentralityCuts();

Int_t GetCentralityVariable(PrimaryVertexInfo *vertex);
Int_t GetUserCentralityVariable(PrimaryVertexInfo *vertex);

Int_t GetCentralityBin(Int_t centralityVariable, Long64_t triggerID=999, Double_t zVertex=999);
UInt_t GetNCentralityBins();

void SetZVertexCuts(Double_t zMin=-999, Double_t zMax=999);
std::pair<double,double> GetZVertexCuts();

void SetRadialVertexCuts(Double_t rMin=0, Double_t rMax=999);
std::pair<double,double> GetRadialVertexCuts();

void SetBeamSpotLocation(Double_t beamX=0, Double_t beamY=0);
std::pair<double,double> GetBeamSpotLocation();

void SetAllowedTriggers(unsigned int triggerID);
void SetAllowedTriggers();
std::vector<unsigned int> GetAllowedTriggers();

void SetPileUpCut(Int_t pileUpCut=99999);
Int_t GetPileUpCut();

Int_t GetMinRapidityIndexOfInterest(Int_t pidIndex = -1);
Int_t GetMaxRapidityIndexOfInterest(Int_t pidIndex = -1);

Int_t GetEnergyIndex(Double_t energy);
std::pair<double,double> GetAverageNpart(Double_t energy, Int_t centBin);
std::pair<double,double> GetAverageNcoll(Double_t energy, Int_t centBin);
std::pair<double,double> GetAverageImpactParameter(Double_t energy, Int_t centBin);
