//This gets called by the makefile to compile all your source code. Add to it as needed!

void makeLibs(TString opt=""){

  //Remove the contents of the binary directory if option "clean" is called
  if (opt == "clean"){
    gSystem->Exec("rm -f bin/*");
    return;
  }

  //Check to make sure the bin directory exists. If not create it.
  if (gSystem->OpenDirectory("bin") == 0)
    gSystem->MakeDirectory("bin");

  //Set the Build and Include Directories
  gSystem->SetBuildDir("$PWD/bin/",true);
  gSystem->SetIncludePath(TString::Format("-I$ROOTSYS/include -I%s/inc",gSystem->pwd()));

  //Compile the contents of the sub-modules
  gSystem->CompileMacro("src/submodules/datacollectorreaderlibs/TrackInfo/TrackInfo.cxx","k");
  gSystem->CompileMacro("src/submodules/datacollectorreaderlibs/PrimaryVertexInfo/PrimaryVertexInfo.cxx","k");
  gSystem->CompileMacro("src/submodules/datacollectorreaderlibs/EventInfo/EventInfo.cxx","k");
  gSystem->CompileMacro("src/submodules/datacollectorreaderlibs/DavisDstReader/DavisDstReader.cxx","k");
  gSystem->CompileMacro("src/submodules/ParticleInfo/ParticleInfo/ParticleInfo.cxx","k");
  gSystem->CompileMacro("src/submodules/strefmultextendedcorr/StRefMultExtendedCorr.cxx","k");
  gSystem->CompileMacro("src/submodules/centralitydatabaseclass/CentralityDatabaseClass/CentralityDatabaseClass.cxx","k");
  gSystem->CompileMacro("src/submodules/SimultaneousFitter/SimFitter.cxx","k");
  
  //Compile your source code
  gSystem->CompileMacro("userfiles/UserCuts.cxx","k");
  gSystem->CompileMacro("src/analysis/utilityFunctions.cxx","k");
  gSystem->CompileMacro("src/analysis/usefulDataStructs.cxx","k");
  gSystem->CompileMacro("src/analysis/skimmerAndBinner.cxx","k");
  gSystem->CompileMacro("src/analysis/dcaBinner.cxx","k");
  gSystem->CompileMacro("src/analysis/fitZTPCUtilities.cxx","k");
  gSystem->CompileMacro("src/analysis/SpectraCorrectionFunctions.cxx","k");
  gSystem->CompileMacro("src/spectraUtilities/SpectraCorrectionCurves.cxx","k");
  gSystem->CompileMacro("src/spectraUtilities/SpectraClass.cxx","k");
  gSystem->CompileMacro("src/spectraUtilities/SpectraFitFunctions.cxx","k");
  gSystem->CompileMacro("src/spectraUtilities/SpectraFitUtilities.cxx","k");
  gSystem->CompileMacro("src/analysis/correctSpectra.cxx","k");
  gSystem->CompileMacro("src/analysis/fitZTPCProtons.cxx","k");
  gSystem->CompileMacro("src/analysis/fitZTPCProtonsTofCut.cxx","k");
  gSystem->CompileMacro("src/analysis/fitZTPCDeuteronsTofCut.cxx","k");
  gSystem->CompileMacro("src/analysis/fitZTOFProtons.cxx","k");
  gSystem->CompileMacro("src/analysis/fitdNdy.cxx","k");
  gSystem->CompileMacro("src/analysis/testCorrectSpectra.cxx","k");

  gSystem->CompileMacro("src/spectraFitting/fitThermalProtonSpectra.cxx","k");
  //gSystem->CompileMacro("src/analysis/fitProtonSpectra.cxx","k");
  gSystem->CompileMacro("src/background/feeddownBackground.cxx","k");
  

  gSystem->CompileMacro("src/analysis/blastwaveFit.cxx","k");
  
  //PID Calibrator
  //  gSystem->CompileMacro("src/pidCalibration/BindEdxInRapidityMtM0.cxx","k");
  //  gSystem->CompileMacro("src/pidCalibration/FitdEdxBetaGamma.cxx","k");
  //  gSystem->CompileMacro("src/pidCalibration/ExtractMeansAndWidths.cxx","k");
  
  //If your code requires any dependencies from the glaubermcmodel submodule
  //make them under the following line
  gSystem->Load("bin/GlauberClass_cxx.so");
  
  
}
