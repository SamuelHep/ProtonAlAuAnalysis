TString FindBinDir(){

  TString binDirPath("");
  if (gSystem->OpenDirectory("bin") != 0)
    binDirPath = TString::Format("%s/bin",gSystem->pwd());
  else if (gSystem->OpenDirectory("../bin") != 0)
    binDirPath = TString::Format("%s/../bin",gSystem->pwd());
  else if (gSystem->OpenDirectory("../../bin") != 0)
    binDirPath = TString::Format("%s/../../bin",gSystem->pwd());
  else {
    cout <<"Could not find bin directory. Have you compiled with make? Are you more than three subdirectories down?" <<endl;
  }

  return binDirPath;
}

void LoadDavisDstReaderModule(){

  TString binDir = FindBinDir();
  if (gSystem->OpenDirectory(Form("%s",binDir.Data())) == 0)
    return;

  gSystem->Load(Form("%s/TrackInfo_cxx.so",binDir.Data()));
  gSystem->Load(Form("%s/PrimaryVertexInfo_cxx.so",binDir.Data()));
  gSystem->Load(Form("%s/EventInfo_cxx.so",binDir.Data()));
  gSystem->Load(Form("%s/DavisDstReader_cxx.so",binDir.Data()));

}

void LoadParticleInfoModule(){

  TString binDir = FindBinDir();
  if (gSystem->OpenDirectory(Form("%s",binDir.Data())) == 0)
    return;

  gSystem->Load(Form("%s/ParticleInfo_cxx.so",binDir.Data()));

}

void LoadCentralityModule(){

  TString binDir = FindBinDir();
  if (gSystem->OpenDirectory(Form("%s",binDir.Data())) == 0)
    return;

  gSystem->Load(Form("%s/GlauberClass_cxx",binDir.Data()));
  gSystem->Load(Form("%s/GlauberModel_cxx",binDir.Data()));
  gSystem->Load(Form("%s/FindBestFitNegativeBinomialParameters_cxx",binDir.Data()));
  gSystem->Load(Form("%s/FindCentralityBinCuts_cxx",binDir.Data()));
  gSystem->Load(Form("%s/DetermineCentralityBin_cxx",binDir.Data()));
  gSystem->Load(Form("%s/ReturnStructs_cxx",binDir.Data()));
  gSystem->Load(Form("%s/FindNpartNcollDistributions_cxx",binDir.Data()));
  gSystem->Load(Form("%s/PrintResults_cxx",binDir.Data()));

}

void LoadUserCutsModule(){

  TString binDir = FindBinDir();
  if (gSystem->OpenDirectory(Form("%s",binDir.Data())) == 0)
    return;

  gSystem->Load(Form("%s/utilityFunctions_cxx",binDir.Data()));
  gSystem->Load(Form("%s/TrackInfo_cxx",binDir.Data()));
  gSystem->Load(Form("%s/PrimaryVertexInfo_cxx",binDir.Data()));
  gSystem->Load(Form("%s/EventInfo_cxx",binDir.Data()));
  gSystem->Load(Form("%s/ParticleInfo_cxx",binDir.Data()));
  gSystem->Load(Form("%s/StRefMultExtendedCorr_cxx",binDir.Data()));
  gSystem->Load(Form("%s/UserCuts_cxx",binDir.Data()));
  
}
