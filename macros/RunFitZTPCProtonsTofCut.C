void RunFitZTPCProtonsTofCut(TString yieldHistoFile,TString outRoot, TString starLibrary, TString imageDir=""){

  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/fitZTPCUtilities_cxx.so");
  gSystem->Load("../bin/fitZTPCProtonsTofCut_cxx.so");
 
  fitZTPCProtonsTofCut(yieldHistoFile,outRoot,starLibrary,imageDir);

}
