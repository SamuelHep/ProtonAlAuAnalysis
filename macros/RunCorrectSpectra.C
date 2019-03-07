//Runs the source code responsible for applying all of the corrections to the spectra
//see src/analysis/correctSpectra.cxx

void RunCorrectSpectra(TString spectraFileName, TString correctionFileName, Int_t pid, Int_t charge,
		       TString starLibrary, Double_t energy, TString eventConfig){
  
  //Load the necessary Librarires
  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/TrackInfo_cxx.so");
  gSystem->Load("../bin/PrimaryVertexInfo_cxx.so");
  gSystem->Load("../bin/EventInfo_cxx.so");
  gSystem->Load("../bin/DavisDstReader_cxx.so");
  gSystem->Load("../bin/StRefMultExtendedCorr_cxx.so");
  gSystem->Load("../bin/UserCuts_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/SpectraCorrectionCurves_cxx.so");
  gSystem->Load("../bin/SpectraCorrectionFunctions_cxx.so");
  gSystem->Load("../bin/SpectraClass_cxx.so");
  gSystem->Load("../bin/correctSpectra_cxx.so");

  SetVariableUserCuts(energy,eventConfig,starLibrary);
  
  CorrectSpectra(spectraFileName,correctionFileName,pid,charge);

}
