void RunFitThermalProtonSpectra(TString spectraFile,TString outFile,TString imageDir, int centBin){

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
  gSystem->Load("../bin/fitThermalProtonSpectra_cxx.so");

  double energy =4.5;
  TString starLibrary = "SL16a";
  TString eventConfig = "FixedTarget2015"; 
  SetVariableUserCuts(energy,eventConfig,starLibrary);

  fitThermalProtonSpectra(spectraFile,outFile,imageDir,centBin);

}
