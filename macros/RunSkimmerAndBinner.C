//This macro runs the compilded skimmerAndBinner binary.
//The source code for skimmerAndBinner.cxx can be found at
//    src/analysis/skimmerAndBinner.cxx

void RunSkimmerAndBinner(TString inputFile, TString starLibrary, Long64_t nEvents=-1, TString outputFile="", TString tofFile="", Int_t tofCut = 0, Double_t energy=0, TString eventConfig=""){

  //Load the Necessary Libraries
  gSystem->Load("../bin/TrackInfo_cxx.so");
  gSystem->Load("../bin/PrimaryVertexInfo_cxx.so");
  gSystem->Load("../bin/EventInfo_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");

  gSystem->Load("../bin/DavisDstReader_cxx.so");
  gSystem->Load("../bin/StRefMultExtendedCorr_cxx.so");
  gSystem->Load("../bin/UserCuts_cxx.so");
  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/skimmerAndBinner_cxx.so");

  //Set the User Cuts
  SetVariableUserCuts(energy,eventConfig,starLibrary);
  
  skimmerAndBinner(inputFile,starLibrary,
		   nEvents,outputFile,tofFile,
		   energy,tofCut,eventConfig);

}
