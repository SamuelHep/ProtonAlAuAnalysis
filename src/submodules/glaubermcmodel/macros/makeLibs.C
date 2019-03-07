//This is a simple ROOT Style "Make" file

void makeLibs(TString opt=""){

  if (opt=="clean"){
    gSystem->Exec("rm -f bin/*");
    return;
  }
  
  //Check to Make sure the bin/ and data/ directories exist...if they don't make them
  if (gSystem->OpenDirectory("bin") == 0)
    gSystem->MakeDirectory("bin");

  if (gSystem->OpenDirectory("data") == 0)
    gSystem->MakeDirectory("data");

  
  gSystem->SetBuildDir("$PWD/bin/",true);
  gSystem->SetIncludePath(TString::Format("-I$ROOTSYS/include -I%s/inc",gSystem->pwd()));
  
  gSystem->CompileMacro("src/GlauberClass.cxx","k-");
  gSystem->CompileMacro("src/ReturnStructs.cxx","k-");
  
  gSystem->CompileMacro("src/GlauberModel.cxx","k-");
  gSystem->CompileMacro("src/FindBestFitNegativeBinomialParameters.cxx","k-");
  gSystem->CompileMacro("src/FindCentralityBinCuts.cxx","k-");
  gSystem->CompileMacro("src/DetermineCentralityBin.cxx","k-");
  gSystem->CompileMacro("src/FindNpartNcollDistributions.cxx","k-");
  gSystem->CompileMacro("src/PrintResults.cxx","k-");
  gSystem->CompileMacro("src/SystematicErrors.cxx","k-");
  gSystem->CompileMacro("src/MakeResultsTable.cxx","k-");
  gSystem->CompileMacro("src/TriggerEfficiency.cxx","k-");
}
