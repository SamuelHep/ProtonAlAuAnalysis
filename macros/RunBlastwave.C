
void RunBlastwave(){


  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/TrackInfo_cxx.so");
  gSystem->Load("../bin/PrimaryVertexInfo_cxx.so");
  gSystem->Load("../bin/EventInfo_cxx.so");
  gSystem->Load("../bin/DavisDstReader_cxx.so");
  gSystem->Load("../bin/StRefMultExtendedCorr_cxx.so");
  gSystem->Load("../bin/UserCuts_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/SpectraFitFunctions_cxx.so");
  gSystem->Load("../bin/SpectraClass_cxx.so");
  gSystem->Load("../bin/SimFitter_cxx.so");
  gSystem->Load("../bin/blastwaveFit_cxx.so");


  blastwaveFit();



}
