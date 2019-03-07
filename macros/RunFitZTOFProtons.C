void RunFitZTOFProtons(TString yieldHistoFile, TString outputSpectraFile, TString starLibrary,
		       TString outputGausFile, Double_t rapidity=-999, TString imageDir=""){

  gSystem->Load("../bin/utilityFunctions_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/fitZTPCUtilities_cxx.so");
  gSystem->Load("../bin/fitZTOFProtons_cxx.so");

  //Find the rapidity index of the user specified rapidity value
  Int_t rapidityIndex = -1;
  if (rapidity != -999)
    rapidityIndex = GetRapidityIndex(rapidity);

  //Check to make sure the rapidityIndex is good before continuing
  if (rapidity != -999 && rapidityIndex < 0){
    cout <<"ERROR in RunFitZTOFProtons.C: value of rapidity is not in rapidity range! EXITING" <<endl;
    return;
  }

  //  imageDir = "../userfiles/AuAu_4_5GeV_2015/analysisCleanup/TOF_Images/";

  cout << "rapidity: " << rapidity << endl;
  cout << "Image Dir: " << imageDir << endl;

  fitZTOFProtons(yieldHistoFile,outputSpectraFile,starLibrary,outputGausFile,rapidityIndex,imageDir);

}
