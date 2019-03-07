void RunFitdNdy(TString blastwaveInputFile, TString outputFile, TString imageDir){

  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/fitZTPCUtilities_cxx.so");
  gSystem->Load("../bin/fitdNdy_cxx.so");

  fitdNdy(blastwaveInputFile,outputFile,imageDir);

}
