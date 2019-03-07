
void RunFeedDownBackgroundAnalysis(TString inputFile, TString outputFile){
  
  //Load the Necessary Libraries
  gSystem->Load("../bin/TrackInfo_cxx.so");
  gSystem->Load("../bin/PrimaryVertexInfo_cxx.so");
  gSystem->Load("../bin/EventInfo_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  
  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/usefulDataStructs_cxx.so");
  gSystem->Load("../bin/DavisDstReader_cxx.so");
  gSystem->Load("../bin/StRefMultExtendedCorr_cxx.so");
  gSystem->Load("../bin/UserCuts_cxx.so");
  gSystem->Load("../bin/feeddownBackground_cxx.so");

  //Set the User Cuts
  SetVariableUserCuts(4.5,"FixedTarget2015","SL16a");
  
  feeddownBackground(inputFile,outputFile);

}
