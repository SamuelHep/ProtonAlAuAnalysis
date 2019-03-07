/**************************************************************
 **  This is an example file showing how you will implement your event, 
 **  vertex, and track cuts in each of the following functions.
 **
 **  IMPORTANT!: 
 **  You should modify the cuts within the functions below to
 **  make the cuts whatever you need/want them to be.
 **
 **  Since track quality cuts are very particular to each analysis
 **  You must specify your own cuts in the IsGoodTrack and 
 **  IsGoodTofTrack functions.
 **
 **  Some event/vertex cuts can be specified by using a combination of 
 **  energy, event configuration, and star library version
 **
 **************************************************************/

#include <iostream>
#include <utility>
#include <vector>

#include <TClonesArray.h>

#include <globalDefinitions.h>
#include <TrackInfo.h>
#include <PrimaryVertexInfo.h>
#include <EventInfo.h>
#include <ParticleInfo.h>
#include <StRefMultExtendedCorr.h>
#include <UserCuts.h>

/**********GLOBAL VARIABLES *************/
/***** DO NOT CHANGE ANYTHING HERE ******/
bool areCutsSet = false;
double energy_GLOBAL = 0;
TString eventConfig_GLOBAL = "";
TString starLibVersion_GLOBAL = "";
std::pair<double,double>  beamSpotLocation_GLOBAL(0.0,0.0);    //(x,y)                                                   
std::pair<double,double>  zVertexCuts_GLOBAL(-999,999);        //(zMin,zMax)                                             
std::pair<double,double>  radialVertexCuts_GLOBAL(0,999);      //(rMin,rMax)                                             
std::vector<unsigned int> allowedTriggers_GLOBAL;              //(min bias trigs)                                        
std::vector<double>       centralityPercents_GLOBAL;           //vector of percent cent bins                             
std::vector<int>          centralityCuts_GLOBAL;               //vector of centrality cuts                               
Int_t                     pileUpCut_GLOBAL(99999);              //cut to eliminate pileup events in most central bin
StRefMultExtendedCorr *refMultExtendedCorr_GLOBAL = NULL;      //Only used if centrality cuts are not manual
/****************************************/
/**** END GLOBAL VARIABLES **************/

//_______________________________________________________________
void SetVariableUserCuts(Double_t energy, TString eventConfig, TString starLibrary){

  //Set the collision Energy, Event Configuration, and STARVER
  //NOTE: These values can either be passed in as arguments to this
  //      function, or they can be set by hand directly below.
  SetCollisionEnergy(energy);
  SetEventConfiguration(eventConfig);
  SetStarLibraryVersion(starLibrary);

  //If no arguments are given in the following functions, the
  //prefconfigured default values are used for the given
  //energy and event configuration.

  //WARNING - These are just some basic, frequently used cuts in
  //          all analyses. You can specify more cuts in the
  //          boolean functions below! It is your responsibility
  //          to make sure the boolean functions below return
  //          true or false based on YOUR cuts. 
  
  SetZVertexCuts(210,212);        //You can also pass in zMin and zMax
  SetRadialVertexCuts();   //You can also pass in rMin and rMax
  SetBeamSpotLocation();   //You can also pass in xBeam and yBeam
  SetAllowedTriggers(1);    //You can also pass in a trigger ID
  SetPileUpCut(125);          //You can also pass in pileUpCut

  //NOTE: You must manually specify your centrality bin percentiles
  //      and cuts unless you are using StRefMultExtnededCorr (collider)
  const Int_t NCENTBINS = 6;
  Int_t centCutsArr[NCENTBINS] = {85,74,65,57,49,41};
  Double_t centPercentArr[NCENTBINS] = {5,10,15,20,25,30};
  SetCentralityCuts(centCutsArr,NCENTBINS);     //You can also pass in an array of ints
  SetCentralityPercents(centPercentArr,NCENTBINS); //You can also pass in an array of doubles

  //Set this boolean to true to keep track of whether the cuts
  //were set or not.
  areCutsSet = true;
  
}

//_______________________________________________________________
Bool_t IsGoodEvent(EventInfo *event){
  //Returns true if this event passes all your cuts
  
  CheckCutStatus();

  if (!event)
    return false;
  
  if (event->GetNPrimaryVertices() <= 0)
    return false;

  //Loop Over the allowed triggers and determine whether
  //this event contains one of them
  Bool_t hasAllowedTrigID(false);
  for (unsigned int iTrig=0; iTrig<allowedTriggers_GLOBAL.size(); iTrig++){
    hasAllowedTrigID = event->ContainsTriggerID(allowedTriggers_GLOBAL.at(iTrig));
    if (hasAllowedTrigID)
      break;      
  }
  
  if (!hasAllowedTrigID && areCutsSet)
    return false;

  return true;
  
}

//_______________________________________________________________
Bool_t IsGoodVertex(PrimaryVertexInfo *vertex){
  
  //Returns true if this vertex passes all your cuts
  CheckCutStatus();
  
  if (vertex->GetVertexIndex() != 0) 
    return false;
  
  //Z Vertex Selection
  if (vertex->GetZVertex() < zVertexCuts_GLOBAL.first ||
      vertex->GetZVertex() > zVertexCuts_GLOBAL.second)
    return false;
  
  //Radial Vertex Selection
  double rVertex = sqrt(pow(vertex->GetXVertex() - beamSpotLocation_GLOBAL.first,2) +
			pow(vertex->GetYVertex() - beamSpotLocation_GLOBAL.second,2));
  if (rVertex < radialVertexCuts_GLOBAL.first ||
      rVertex > radialVertexCuts_GLOBAL.second)
    return false;
  /*
  //Number of ToF Mathces
  if (vertex->GetNTofMatches() < 3)
    return false;
  */
  return true;
}

//_______________________________________________________________
Bool_t IsGoodTrack(TrackInfo *track, Bool_t useThisCut){
  
  //Returns true if this track passes all your cuts
  CheckCutStatus();
  /*
  if (track->GetFlag() < 0)
    return false;

  if (track->GetFlag() >= 1000)
    return false;

  if (track->GetHitsFit() < 15)
    return false;
  */
  if (track->GetdEdxHits() < 1)
    return false;
  
  if (track->GetNFrac() < 0.52)
    return false;

    
  if (useThisCut){
    if (fabs(track->GetGlobalDCA()->Mag()) > 3.0)
      return false;
  }
  
  return true;
}

//_______________________________________________________________
Bool_t IsGoodTofTrack(TrackInfo *track){

  //Returns true if this track passes all your cuts
  CheckCutStatus();
  
  if (track->GetTofMatchFlag() < 1)
    return false;

  if (1.0/track->GetTofBeta() < 0.0)
    return false;

  if (fabs(track->GetTofYLocal()) > 1.6)
    return false;

  if (fabs(track->GetTofZLocal()) > 3.0)
    return false;

  return true;
  
}

//_______________________________________________________________
Int_t GetUserCentralityVariable(PrimaryVertexInfo *vertex){

  //Use this function if you would like to compute your own
  //centrality variable.

  //!!!!!!!!VERY IMPORTANT NOTE!!!!!!!
  //When you are writing code you should call the function
  //GetCentralityVariable instead of this function. It is more
  //general and will call this function if needed.

  Int_t totalMultiplicity(0);
  TrackInfo *track;
  
  for (Int_t iTrack=0; iTrack<vertex->GetNPrimaryTracks(); iTrack++){
    
    track = vertex->GetPrimaryTrack(iTrack);
    
    if (IsGoodTrack(track, false))
      totalMultiplicity++;
    
  }
  
  return totalMultiplicity;
}

//________________________________________________________________
Int_t GetEnergyIndex(Double_t energy){
  
  //Get the index of the requested energy.
  //Note: This function can be useful if you
  //are using multiple energies and would like to
  //switch between them.
  
  Int_t energyIndex(-1);
  const int nEnergies(7);
  Double_t energies[nEnergies] = {7.7,11.5,14.5,19.6,27.0,39.0,62.4};
  for (int i=0; i<nEnergies; i++){
    if (fabs(energy - energies[i]) < .2){
      energyIndex = i;
      break;
    }
  }
  
  //Check that the energy index was assigned
  if (energyIndex < 0){
    cout <<"ERROR: UserCuts::GetEnergyIndex - Requested energy is not found! EXITING" <<endl;
    exit (EXIT_FAILURE);
  }

  return energyIndex;
}

//________________________________________________________________
std::pair<double,double> GetAverageNpart(Double_t energy, Int_t centBin){

  //Return the average Number of participant nucleons
  //for the requested energy and centrality bin
  
  Double_t avgNpart[7][16] =
    { {12,16,22,29,38,50,63,79,98, 119,144,172,205,244,289,341}, //AuAu 7.7
      {11,16,21,29,38,50,63,79,98, 119,144,173,205,244,289,341}, //AuAu 11.5
      {11,16,21,29,39,50,64,80,98, 120,144,173,206,244,290,342}, //AuAu 14.5
      {12,16,21,29,39,50,64,79,98, 120,145,174,208,246,291,342}, //AuAu 19.6
      {11,15,21,29,39,50,63,79,98, 121,146,175,208,246,292,343}, //AuAu 27.0
      {11,16,21,29,39,50,64,80,100,121,147,175,209,248,293,345}, //AuAu 39.0
      {11,16,22,30,39,51,65,81,100,122,148,177,211,250,284,346}  //AuAu 62.4
    };

  Double_t avgNpartErr[7][16] =
    { {2,2,2,2,2,2,2,3,3,4,4, 5, 6, 8, 9,5}, //AuAu 7.7
      {2,2,2,2,2,2,3,3,4,5,6, 7, 8, 9,11,6}, //AuAu 11.5
      {1,2,2,2,2,2,3,3,4,5,6, 7, 8,10,12,7}, //AuAu 14.5
      {1,2,2,2,2,2,3,3,4,5,5, 6, 8, 9,11,6}, //AuAu 19.6
      {1,1,2,2,2,2,3,3,4,5,6, 7, 8,10,11,6}, //AuAu 27.0
      {1,1,1,2,2,2,3,4,4,5,6, 7, 8,10,12,6}, //AuAu 39.0
      {1,1,2,2,3,3,4,5,7,7,8,11,13,16,19,9} //AuAu 62.4
    };

  Int_t energyIndex = GetEnergyIndex(energy);
  
  //Nothing needs to be converted if we are using all 16 bins
  if (GetNCentralityBins() == 16)
    return std::make_pair(avgNpart[energyIndex][centBin],avgNpartErr[energyIndex][centBin]);

  //**The rest of this function assumes 9 centrality bins are being used
  if (centBin == 8)
    return std::make_pair(avgNpart[energyIndex][15],avgNpartErr[energyIndex][15]);
  else if (centBin == 7)
    return std::make_pair(avgNpart[energyIndex][14],avgNpartErr[energyIndex][14]);

  Int_t binA, binB;
  if (centBin == 6){
    binA = 13; binB = 12;
  }
  else if (centBin == 5){
    binA = 11; binB = 10;
  }
  else if (centBin == 4){
    binA = 9; binB = 8;
  }
  else if (centBin == 3){
    binA = 7; binB = 6;
  }
  else if (centBin == 2){
    binA = 5; binB =4;
  }
  else if (centBin == 1){
    binA = 3; binB = 2;
  }
  else if (centBin == 0){
    binA = 1; binB =0;
  }

  Double_t binMean =    (avgNpart[energyIndex][binA]+avgNpart[energyIndex][binB])/2.0;
  Double_t binMeanErr = sqrt(pow(avgNpartErr[energyIndex][binA],2) +
			     pow(avgNpartErr[energyIndex][binB],2))/2.0;

  return std::make_pair(binMean, binMeanErr);
  
}

//________________________________________________________________
std::pair<double,double> GetAverageNcoll(Double_t energy, Int_t centBin){
  
  //Return the average Number of participant nucleons
  //for the requested energy and centrality bin
  /*
  Double_t avgNcoll[7][9] =
    { {}, //AuAu 7.7
      {}, //AuAu 11.5
      {}. //AuAu 14.5
      {}, //AuAu 19.6
      {}, //AuAu 27.0
      {}, //AuAu 39.0
      {}, //AuAu 62.4
    };
  
  Double_t avgNcollErr[7][9] =
    { {}, //AuAu 7.7
      {}, //AuAu 11.5
      {}. //AuAu 14.5
      {}, //AuAu 19.6
      {}, //AuAu 27.0
      {}, //AuAu 39.0
      {}, //AuAu 62.4
    };
  
  Int_t energyIndex = GetEnergyIndex(energy);
  return std::make_pair(avgNcoll[energyIndex][centBin],avgNcollErr[energyIndex][centBin]);
  */
}

//________________________________________________________________
std::pair<double,double> GetAverageImpactParameter(Double_t energy, Int_t centBin){

  //Return the average Number of participant nucleons
  //for the requested energy and centrality bin
  /*
  Double_t avgImp[7][9] =
    { {}, //AuAu 7.7
      {}, //AuAu 11.5
      {}. //AuAu 14.5
      {}, //AuAu 19.6
      {}, //AuAu 27.0
      {}, //AuAu 39.0
      {}, //AuAu 62.4
    };
  
  Double_t avgImpErr[7][9] =
    { {}, //AuAu 7.7
      {}, //AuAu 11.5
      {}. //AuAu 14.5
      {}, //AuAu 19.6
      {}, //AuAu 27.0
      {}, //AuAu 39.0
      {}, //AuAu 62.4
    };
  
  Int_t energyIndex = GetEnergyIndex(energy);
  return std::make_pair(avgImp[energyIndex][centBin],avgImpErr[energyIndex][centBin]);
  */
}


//________________________________________________________________
/*****************************************************************
 **                                                             **
 **    YOU SHOULD NOT NEED TO MODIFY ANY CODE BELOW, UNLESS     **
 **         YOU ARE DOING SOMETHING NEW/EXPERIMENTAL.           **
 **                                                             **
 *****************************************************************/

//_______________________________________________________________
void SetCollisionEnergy(Double_t energy){
  energy_GLOBAL = energy;
}

//_______________________________________________________________
Double_t GetCollisionEnergy(){
  return energy_GLOBAL;
}

//_______________________________________________________________
void SetEventConfiguration(TString eventConfig){
  eventConfig_GLOBAL = eventConfig;
}

//_______________________________________________________________
TString GetEventConfiguration(){
  return eventConfig_GLOBAL;
}

//_______________________________________________________________
void SetStarLibraryVersion(TString starLibVer){
  starLibVersion_GLOBAL = starLibVer;
}

//_______________________________________________________________
TString GetStarLibraryVersion(){
  return starLibVersion_GLOBAL;
}

//_______________________________________________________________
void SetCentralityPercents(Double_t *centPercentArr, const int nCentBins){

  //Check that the cuts haven't already been set
  if (areCutsSet){
    cout <<"WARNING - UseCuts::SetCentralityPercents should only be called once! No Action taken.\n";
    return;
  }

  //If the user has specified values
  if (centPercentArr != NULL && nCentBins != -1){
    for (unsigned int i=0; i<(unsigned int)nCentBins; i++){
      centralityPercents_GLOBAL.push_back(centPercentArr[i]);
    }
    return;
  }

  //Otherwise assume that StRefMultExtendedCorr is being used
  centralityPercents_GLOBAL.push_back(80);
  centralityPercents_GLOBAL.push_back(70);
  centralityPercents_GLOBAL.push_back(60);
  centralityPercents_GLOBAL.push_back(50);
  centralityPercents_GLOBAL.push_back(40);
  centralityPercents_GLOBAL.push_back(30);
  centralityPercents_GLOBAL.push_back(20);
  centralityPercents_GLOBAL.push_back(10);
  centralityPercents_GLOBAL.push_back(5);
    

}

//_______________________________________________________________
std::vector<double> GetCentralityPercents(){

  //Returns a vector containing the matching percent bins of
  //the centrality cuts in GetCentralityCuts() below
  CheckCutStatus();
  
  return centralityPercents_GLOBAL;

}

//_______________________________________________________________
std::vector<int> GetCentralityCuts(){

  //Returns a the vector of centrality Cuts
  CheckCutStatus();
  
  return centralityCuts_GLOBAL;

}

//_______________________________________________________________
UInt_t GetNCentralityBins(){

  //Return the number of centrality bins
  CheckCutStatus();
  
  //If the cuts have been speficied manually return the
  //size of the centrality cut vector
  if (areCutsSet && centralityCuts_GLOBAL.size() > 0)
    return centralityCuts_GLOBAL.size();

  //Otherwise 9 bins will be used from StRefMultExtendedCorr
  return 9;

}

//_______________________________________________________________
Int_t GetCentralityVariable(PrimaryVertexInfo *vertex){

  //Returns the centrality variable you would like to cut on

  //Make sure the cuts are set
  CheckCutStatus();
  
  //Compute the Centrality Variable you want to cut on if
  //you are not using a value already stored in PrimaryVertexInfo


  //First Check if StRefMultExtendedCorr is being used
  if (refMultExtendedCorr_GLOBAL){

    //For Vertices with |Vz| <= 95 cm use the standard refmult
    if (fabs(vertex->GetZVertex()) <= 95)
      return vertex->GetRefMult();

    //For Vertices with Vz < -95 use RefMultPosY
    else if (vertex->GetZVertex() < -95)
      return vertex->GetRefMultPosY();

    //For Vertices with Vz > 95 use RefMultNegY
    else if (vertex->GetZVertex() > 95)
      return vertex->GetRefMultNegY();
    
  }//End StRefMultExtendedCorr

  //If you are not using StRefMultExtendedCorr for your centrality
  //determination then you can compute your own variable here
  return GetUserCentralityVariable(vertex);
}

//______________________________________________________________
void SetCentralityCuts(Int_t *centCutsArr, const int nCentBins){

  //Set the centrality cuts and percentage bins here.

  //Check that the cuts haven't already been set
  if (areCutsSet){
    cout <<"WARNING - UseCuts::SetCentralityCuts should only be called once! No Action taken.\n";
    return;
  }

  //IMPORTANT!!!! There are two methods in use depending on the analysis.
  //              If you specify a manual cut here it will override StRefMultExtendedCorr.

  //------------ METHOD 1 - MANUAL DEFINITIONS --------------
  //*** Add as many calls to push_back as you have bins ***
  //If the user has specified values
  if (centCutsArr != NULL && nCentBins != -1){
    for (unsigned int i=0; i<(unsigned int)nCentBins; i++){
      centralityCuts_GLOBAL.push_back(centCutsArr[i]);
    }
  }

  if (centralityCuts_GLOBAL.size() > 0)
    return;

  //------------- METHOD 2 - STREFMULTEXTENDEDCORR ------------
  refMultExtendedCorr_GLOBAL = new StRefMultExtendedCorr();
  refMultExtendedCorr_GLOBAL->Initialize(energy_GLOBAL);
 
}

//______________________________________________________________
Int_t GetCentralityBin(Int_t centralityVariable, Long64_t triggerID, Double_t zVertex){

  //Returns the centrality Bin of a primary vertex with
  //quantity centralityVariable from GetCentralityVariable
  CheckCutStatus();
  
  //NOTE: triggerID and zVertex are set to 999 by default. You only need to
  //      use them if you are suing StRefMultExtendedCorr

  //If StRefMultExtendedCorr is being used then let it determine the bin
  if (areCutsSet && centralityCuts_GLOBAL.size() == 0 && triggerID != 999)
    return refMultExtendedCorr_GLOBAL->GetCentralityBin9(triggerID,zVertex,centralityVariable);

  //Check to make sure there is at least one centrality bin, if there is none do not bin in centrality
  if (centralityCuts_GLOBAL.size() == 0 ){
    return 0;
  }

  //Check that the centrality variable does not exceed the pile up cut
  if (pileUpCut_GLOBAL != 99999 && centralityVariable > pileUpCut_GLOBAL)
    return -2;

  //Check the high bin
  if (centralityVariable >= centralityCuts_GLOBAL.at(0))
    return 0;
  
  //Loop Over the remaining bins and find where this vertex should be binned
  for (int iCentBin=1; iCentBin<(int)centralityCuts_GLOBAL.size(); iCentBin++){
    
    if (centralityVariable >= centralityCuts_GLOBAL.at(iCentBin) &&
	centralityVariable < centralityCuts_GLOBAL.at(iCentBin-1))
      return iCentBin;
    
  }

  //If no centrality bin was found then return -1. In this case
  //the vertex will always be skipped;
  return -1;

}

//______________________________________________________________
std::pair<double,double> GetZVertexCuts(){

  if (!areCutsSet)
    SetZVertexCuts();

  return zVertexCuts_GLOBAL;
  
}

//______________________________________________________________
void SetZVertexCuts(Double_t zMin, Double_t zMax){

  //Here you can specify the zVertex cuts for your event configuration

  //You can either pass in values for zMin and zMax
  if (zMin != -999 && zMax != 999){
    zVertexCuts_GLOBAL.first = zMin;
    zVertexCuts_GLOBAL.second = zMax;
    return;
  }
    
  //Or you can rely on the event configuration
  if (eventConfig_GLOBAL.Contains("FixedTarget2015")){
    zVertexCuts_GLOBAL.first  = 210;
    zVertexCuts_GLOBAL.second = 212;
  }

  else if (eventConfig_GLOBAL.Contains("FixedTarget2014")){
    zVertexCuts_GLOBAL.first  = 209;
    zVertexCuts_GLOBAL.second = 211;
  }

  else if (eventConfig_GLOBAL.Contains("BeamPipePosY")){
    zVertexCuts_GLOBAL.first  = -200;
    zVertexCuts_GLOBAL.second = -150;
  }

  else if (eventConfig_GLOBAL.Contains("BeamPipeNegY")){
    zVertexCuts_GLOBAL.first  = 150;
    zVertexCuts_GLOBAL.second = 200;
  }

  else if (eventConfig_GLOBAL.Contains("ColliderCenter")){
    zVertexCuts_GLOBAL.first  = -30;
    zVertexCuts_GLOBAL.second = 30;

    //Use a large window for 7.7 Gev for more statistics
    if (fabs(energy_GLOBAL-7.7) < 0.1){
      zVertexCuts_GLOBAL.first = -50;
      zVertexCuts_GLOBAL.second = 50;
    }
  }

  else if (eventConfig_GLOBAL.Contains("ColliderPosY")){
    zVertexCuts_GLOBAL.first  = -200;
    zVertexCuts_GLOBAL.second = -100;
  }

  else if (eventConfig_GLOBAL.Contains("ColliderNegY")){
    zVertexCuts_GLOBAL.first  = 100;
    zVertexCuts_GLOBAL.second = 200;
  }

  
}

//______________________________________________________________
std::pair<double,double> GetRadialVertexCuts(){

  if (!areCutsSet)
    SetRadialVertexCuts();

  return radialVertexCuts_GLOBAL;
  
}

//______________________________________________________________
void SetRadialVertexCuts(Double_t rMin, Double_t rMax){

  //Set the radial Vertex Cuts here

  //You can either pass in values
  if (rMin != 0 && rMax != 999){
    radialVertexCuts_GLOBAL.first = rMin;
    radialVertexCuts_GLOBAL.second = rMax;
  }
  
  //Or you can rely on the preconfiguration for the energy
  if (fabs(energy_GLOBAL-7.7) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 2.0;
  }

  else if (fabs(energy_GLOBAL-11.5) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 2.0;
  }
  
  else if (fabs(energy_GLOBAL-14.5) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 0.5;
  }
  
  else if (fabs(energy_GLOBAL-19.6) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 2.0;
  }
  
  else if (fabs(energy_GLOBAL-27.0) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 2.0;
  }
  
  else if (fabs(energy_GLOBAL-39.0) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 2.0;
  }
  
  else if (fabs(energy_GLOBAL-62.4) < 0.1){
    radialVertexCuts_GLOBAL.first  = 0.0;
    radialVertexCuts_GLOBAL.second = 2.0;
  }
  
}

//_______________________________________________________________
std::pair<double,double> GetBeamSpotLocation(){

  if (!areCutsSet)
    SetBeamSpotLocation();
  
  return beamSpotLocation_GLOBAL;
  
}
  
//______________________________________________________________
void SetBeamSpotLocation(Double_t beamX, Double_t beamY){

  //Set the location of the beams spot, this is the xy coordiated
  //around which the radial vertex selection will be made

  //You can either pass in values
  if (beamX != 0 && beamY != 0){
    beamSpotLocation_GLOBAL.first = beamX;
    beamSpotLocation_GLOBAL.second = beamY;
  }

  //Or you can rely on the preconfigured values for the
  //given energy you are using
  if (fabs(energy_GLOBAL-7.7) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.20;
    beamSpotLocation_GLOBAL.second = -0.25;
  }

  else if (fabs(energy_GLOBAL-11.5) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.10;
    beamSpotLocation_GLOBAL.second = -0.23;
  }

  else if (fabs(energy_GLOBAL-14.5) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.0;
    beamSpotLocation_GLOBAL.second = -0.9;
  }
  
  else if (fabs(energy_GLOBAL-19.6) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.37;
    beamSpotLocation_GLOBAL.second = -0.03;
  }
  
  else if (fabs(energy_GLOBAL-27.0) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.13;
    beamSpotLocation_GLOBAL.second = -0.07;
  }
  
  else if (fabs(energy_GLOBAL-39.0) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.30;
    beamSpotLocation_GLOBAL.second = -0.07;
  }
  
  else if (fabs(energy_GLOBAL-62.4) < 0.1){
    beamSpotLocation_GLOBAL.first  = 0.39;
    beamSpotLocation_GLOBAL.second = -0.06;
  }
  

}

//__________________________________________________________________________
std::vector<unsigned int> GetAllowedTriggers(){

  if (!areCutsSet)
    SetAllowedTriggers();
  
  return allowedTriggers_GLOBAL;
  
}

//__________________________________________________________________________
void SetAllowedTriggers(unsigned int triggerID){

  //Overload for the case in which the user has just one trigger ID

  //Check to make sure this has not already been called
  if (areCutsSet){
    cout <<"WARNING - UserCuts:: SetAllowedTriggers() has already been called! It should only be called once. No action taken!\n";
    return;
  }
  
  allowedTriggers_GLOBAL.push_back(triggerID);
  
}

//__________________________________________________________________________
void SetAllowedTriggers(){

  //Check to make sure this has not already been called
  if (areCutsSet){
    cout <<"WARNING - UserCuts:: SetAllowedTriggers() has already been called! It should only be called once. No action taken!\n";
    return;
  }

  if (eventConfig_GLOBAL.Contains("FixedTarget2015")){
    allowedTriggers_GLOBAL.push_back(1);
  }
  
  else if (fabs(energy_GLOBAL-7.7) < 0.1){
    allowedTriggers_GLOBAL.push_back(290001);
    allowedTriggers_GLOBAL.push_back(290004);
  }

  else if (fabs(energy_GLOBAL-11.5) < 0.1){
    allowedTriggers_GLOBAL.push_back(310014);
  }

  else if (fabs(energy_GLOBAL-14.5) < 0.1){
    allowedTriggers_GLOBAL.push_back(440015);
  }
  
  else if (fabs(energy_GLOBAL-19.6) < 0.1){
    allowedTriggers_GLOBAL.push_back(340001);
    allowedTriggers_GLOBAL.push_back(340011);
    allowedTriggers_GLOBAL.push_back(340021);
  }
  
  else if (fabs(energy_GLOBAL-27.0) < 0.1){
    allowedTriggers_GLOBAL.push_back(360001);
  }
  
  else if (fabs(energy_GLOBAL-39.0) < 0.1){
    allowedTriggers_GLOBAL.push_back(280001);
  }
  
  else if (fabs(energy_GLOBAL-62.4) < 0.1){
    allowedTriggers_GLOBAL.push_back(270001);
    allowedTriggers_GLOBAL.push_back(270011);
    allowedTriggers_GLOBAL.push_back(270021);
  }

  if (allowedTriggers_GLOBAL.size() == 0){
    cout <<"WARNING - UserCuts::SetAllowedTriggers() - No triggers found for this energy/eventConfig combination! No trigger cut will be applied.\n";
  }
  
}

//_______________________________________________________________
Int_t GetPileUpCut(){

	if (!areCutsSet)
		SetPileUpCut();

	return pileUpCut_GLOBAL;

}

//_______________________________________________________________
void SetPileUpCut(Int_t pileUpCut){
	
  //Here you can specify the pile up cut for your event configuration
  
  //You can pass in a value for pileUpCut
  if (pileUpCut != 99999){
    pileUpCut_GLOBAL = pileUpCut;
    return;
  }
  
}

//_______________________________________________________________
Int_t GetMinRapidityIndexOfInterest(Int_t pidIndex){

  //By default (pidIndex < 0) the min rapidity index for the
  //config is returned without regard to species.
  
  if (eventConfig_GLOBAL.Contains("FixedTarget2015")){
    if (pidIndex < 0)
      return 0;
    else if (pidIndex == 0) //Pion
      return 0;// nRapidityBins-1;
    else if (pidIndex == 1) //Kaon
      return 0;//nRapidityBins-1;
    else if (pidIndex == 2) //Proton
      return 0;//nRapidityBins-1;
  }
  
  else if (eventConfig_GLOBAL.Contains("ColliderCenter")){
    if (pidIndex < 0)
      return 10;
    else if (pidIndex == 0) //Pion
      return 11;
    else if (pidIndex == 1) //Kaon
      return 14;
    else if (pidIndex == 2) //Proton
      return 14;
  }

  else if (eventConfig_GLOBAL.Contains("ColliderPosY")){
    if (pidIndex < 0)
      return 18;
    else if (pidIndex == 0) //Pion
      return 20;
    else if (pidIndex == 1) //Kaon
      return 20;
    else if (pidIndex == 2) //Proton
      return 20;
  }

  else if (eventConfig_GLOBAL.Contains("ColliderNegY")){
    if (pidIndex < 0)
      return 5;
    else if (pidIndex == 0) //Pion
      return 7;
    else if (pidIndex == 1) //Kaon
      return nRapidityBins-1;
    else if (pidIndex == 2) //Proton
      return 8;
    
  }

  //Return value if eventconfig or pidIndex is not recognized
  return 0;     

}

//_______________________________________________________________
Int_t GetMaxRapidityIndexOfInterest(Int_t pidIndex){

  //By default if pidIndex < 0 then return the max rapidity index
  //of interest for the event config with out regard to species

  if (eventConfig_GLOBAL.Contains("FixedTarget2015")){
    if (pidIndex < 0)
      return nRapidityBins-1;
    else if (pidIndex == 0) //Pion
      return 20;//nRapidityBins-1;
    else if (pidIndex == 1) //Kaon
      return 20;//nRapidityBins-1;
    else if (pidIndex == 2) //Proton
      return 20;//nRapidityBins-1;
  }
  
  else if (eventConfig_GLOBAL.Contains("ColliderCenter")){
    if (pidIndex < 0)
      return 30;
    else if (pidIndex == 0) //Pion
      return 29;
    else if (pidIndex == 1) //Kaon
      return 26;
    else if (pidIndex == 2) //Proton
      return 26;
  }
  
  else if (eventConfig_GLOBAL.Contains("ColliderPosY")){
    if (pidIndex < 0)
      return 35;
    else if (pidIndex == 0) //Pion
      return 33;
    else if (pidIndex == 1) //Kaon
      return nRapidityBins-1;
    else if (pidIndex == 2) //Proton
      return 33;
    
  }
  
  else if (eventConfig_GLOBAL.Contains("ColliderNegY")){
    if (pidIndex < 0)
      return 20;
    else if (pidIndex == 0) //Pion
      return 20;
    else if (pidIndex == 1) //Kaon
      return 20;
    else if (pidIndex == 2) //Proton
      return 20;
  }

  //Return mode if eventConfig or pidIndex is not recognized
  return nRapidityBins-1;

}


//_______________________________________________________________
void CheckCutStatus(){

  //Check that the function SetVariableUserCuts has been called
  if (!areCutsSet){
    cout <<"WARNING - UserCuts::SetVariableUserCuts() - User cuts have not been defined!!!\n"
	 <<"          Have you called SetVariableUserCuts(energy,eventConfig)?\n";
  }

}
