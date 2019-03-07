void RunFitZTPCProtons(TString yieldHistoFile, TString outputSpectraFile, TString protonwidth, TString deuteronwidth,TString ampRatioFile, TString starLibrary,
		       Int_t centralityIndex, Double_t rapidity=-999, TString imageDir="",Int_t first=0,Int_t last=0){

  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/fitZTPCUtilities_cxx.so");
  gSystem->Load("../bin/fitZTPCProtons_cxx.so");

  //Find the rapidity index of the user specified rapidity value
  Int_t rapidityIndex = -1;
  if (rapidity != -999)
    rapidityIndex = GetRapidityIndex(rapidity);

  //Check to make sure the rapidityIndex is good before continuing
  if (rapidity != -999 && rapidityIndex < 0){
    cout <<"ERROR in RunFitZTPCProtons.C: value of rapidity is not in rapidity range! EXITING" <<endl;
    return;
  }
  

  fitZTPCProtons(yieldHistoFile,outputSpectraFile,protonwidth,deuteronwidth,ampRatioFile,starLibrary,centralityIndex,rapidityIndex,imageDir,first,last);

}
