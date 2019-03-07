#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <cmath>

#include <TObject.h>
#include <TF1.h>
#include <TFile.h>
#include <TString.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TMath.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include "ParticleInfo.h"

ClassImp(ParticleInfo);

//______________________________________________________
ParticleInfo::ParticleInfo(){

  //This constructor creates a particle info object, but
  //does not create any of the pid functions. It will
  //still let you get particle masses, names, ect... though.

  canonicalTPCWidth = 0.07;
  canonicalTOFWidth = 0.012;

  //Set All the Pid Functions to NULL
  bichselFunction = NULL;
  bichselFunctionPart[0] = NULL;
  bichselFunctionPart[1] = NULL;
  for (Int_t i=0; i<nDefinedParticles; i++){
    bichselPIDFunctions[i] = NULL;
    inverseBetaPIDFunctions[i] = NULL;
  }
  bichselRootFile = NULL;
  pidCalibrationFile = NULL;
  gausFunc = NULL;
  pidFunctionsCreated = false;
}

//______________________________________________________
ParticleInfo::ParticleInfo(TString starLibraryVersion, Bool_t createAll, const Int_t dEdxType){

  //This constructor forms the bichsel functions for each particle speices using the
  //parameterizations for starLibraryVersion

  canonicalTPCWidth = 0.07;
  canonicalTOFWidth = 0.012;
  
  libraryVersion = starLibraryVersion;
  type = dEdxType;

  pidFunctionsCreated = false;
  useEmpiricalFits = false;
  empiricalSpeciesInterest = -1;
  empiricalRapidityIndex   = -1;
  
  minBetaGamma = 0.1;
  midBetaGamma = 10.0;
  maxBetaGamma = 10000;

  //Set All the Pid Functions to NULL
  bichselFunction = NULL;
  bichselFunctionPart[0] = NULL;
  bichselFunctionPart[1] = NULL;
  for (Int_t i=0; i<nDefinedParticles; i++){
    bichselPIDFunctions[i] = NULL;
    inverseBetaPIDFunctions[i] = NULL;
  }
  bichselRootFile = NULL;
  pidCalibrationFile = NULL;

  gausFunc = new TF1("","gaus(0)",-100,100);
  
  //Read the Bichsel File
  ReadBichselFile();

  //If the option to create all pid Functions is true then create them all
  if (createAll)
    CreateAllDefinedParticlePidFunctions();
  
}

//_______________________________________________________
ParticleInfo::ParticleInfo(TString starLibraryVersion, Bool_t createAll, const int rapidityBins, TString pidCalibFile){

  //This constructor forms the bichsel functions for each particle speices using
  //the emperical 

  canonicalTPCWidth = 0.07;
  canonicalTOFWidth = 0.012;
  
  libraryVersion=starLibraryVersion;
  type = -1;

  pidFunctionsCreated = false;
  useEmpiricalFits = true;
  empiricalSpeciesInterest = -1;
  empiricalRapidityIndex   = -1;
  
  pidCalibrationFile = new TFile(pidCalibFile,"READ");
  
  minBetaGamma = 0.1;
  maxBetaGamma = 10000;

  //Set All the Pid Functions to NULL
  bichselFunction = NULL;
  bichselFunctionPart[0] = NULL;
  bichselFunctionPart[1] = NULL;
  for (Int_t i=0; i<nDefinedParticles; i++){
    bichselPIDFunctions[i] = NULL;
    inverseBetaPIDFunctions[i] = NULL;
    
    bichselEmpiricalFit.push_back(std::vector<TF1 *> (rapidityBins,(TF1 *)NULL));
    bichselEmpiricalPIDFunctions.push_back(std::vector<std::vector<TF1 *> > (rapidityBins, std::vector<TF1 *>
									     (nDefinedParticles, (TF1 *) NULL)));

  }
  bichselRootFile = NULL;

  gausFunc = new TF1("","gaus(0)",-100,100);
  
  //Read the Bichsel File
  ReadBichselFile();
  
  //If the option to create all pid Functions is true then create the pid Functions
  if (createAll)
    CreateAllDefinedParticlePidFunctions();
  
  
}

//______________________________________________________
ParticleInfo::~ParticleInfo(){

  if (bichselFunction)
    delete bichselFunction;

  for (int iParticle=0; iParticle<nDefinedParticles; iParticle++){

    if (bichselPIDFunctions[iParticle])
      delete bichselPIDFunctions[iParticle];
    if (inverseBetaPIDFunctions[iParticle])
      delete inverseBetaPIDFunctions[iParticle];

    if (bichselEmpiricalFit.size() > 0){
      for (unsigned int yIndex=0; yIndex<bichselEmpiricalPIDFunctions.at(iParticle).size(); yIndex++){

	if (bichselEmpiricalFit.at(iParticle).at(yIndex))
	  delete bichselEmpiricalFit.at(iParticle).at(yIndex);
	
	for (int iSubParticle=0; iSubParticle<nDefinedParticles; iSubParticle++){
	  if (bichselEmpiricalPIDFunctions.at(iParticle).at(yIndex).at(iSubParticle))
	    delete bichselEmpiricalPIDFunctions.at(iParticle).at(yIndex).at(iSubParticle);
	}
      }
    }
    
  }

  if(gausFunc)
    delete gausFunc;

  if (bichselRootFile){
    bichselRootFile->Close();
    delete bichselRootFile;
  }

  if (pidCalibrationFile){
    pidCalibrationFile->Close();
    delete pidCalibrationFile;
  }
  
}


//_______________________________________________________
Double_t ParticleInfo::BichselPrototype(Double_t *x, Double_t *par){
  
  //If The Bichsel Graph has not been created then create it
  if (!bichselFunction){
    std::cout <<"ERROR: ParticleInfo::BichselPrototype() bichselFunction NOT set. If you want to use the\n"
         <<"       PID Functions (Bichsel/TOF) then use the constructor:\n"
         <<"       ParticleInfo::ParticleInfo(TString starLibraryVersion, createAll=true)\n"
	 <<"       or call ParticleInfo::CreateAllDefinedParticlePidFunctions() prior to this method.\n"
         <<"       EXITING!\n";
    exit (EXIT_FAILURE);
  }

  Double_t mass = par[0];
  Double_t charge = par[1];

  Double_t betaGamma = x[0] * charge;
   
  if (par[0] <= 0)
    return bichselFunction->Eval(betaGamma);
  
  return pow(charge,2) * bichselFunction->Eval(betaGamma/par[0]);

}

//_______________________________________________________
Double_t ParticleInfo::InverseBetaPrototype(Double_t *x, Double_t *par){

  //Can be used to create a function for Inverse Beta as a function
  //of momemtum for a given mass assumption

  Double_t mass = par[0];
  Double_t charge = par[1];

  Double_t p = x[0] * charge;
  Double_t beta = sqrt(pow(p,2)/(pow(mass,2)+pow(p,2)));

  return 1.0/beta;

}

//_______________________________________________________
Double_t ParticleInfo::BichselEmpiricalFitPrototype(Double_t *x, Double_t *par){

  Double_t betaGamma = x[0];
  
  Double_t m1 = par[0];
  Double_t m2 = par[1];
  Double_t m3 = par[2];
  Double_t m4 = par[3];
  Double_t m5 = par[4];
  Double_t m6 = par[5];
  Double_t m7 = par[6];
  
  Double_t term1 = TMath::Power(1.0 + TMath::Power(1.0/betaGamma,2),m2);
  Double_t term2 = TMath::Power(TMath::Abs(TMath::Log(m3*TMath::Power(betaGamma,2))),m4);
  term2 += m5 * TMath::Power(1.0+TMath::Power(1.0/betaGamma,2),m6);
  
  return m1 * term1 * term2 - m7;
}

//_______________________________________________________
Bool_t ParticleInfo::FindBichselFile(TString fileName, TString path, TString &foundPath){

  //Recursively search the directory tree until the file
  //with the name fileName is found, return the full path
    
  //Get the Contents of path
  TSystemDirectory dir(path.Data(),path.Data());
  TList *list = dir.GetListOfFiles();
  
  //Look through the contents of pwd until we find the file
  if (list){
    
    TSystemFile *file;
    TString name;
    TIter next(list);
    
    while((file=(TSystemFile*)next())){
      name = file->GetName();

      if (name.Contains(".") && !name.Contains(".root"))
	continue;

      if (file->IsDirectory()){
	if (FindBichselFile(fileName, path+"/"+name, foundPath))
	  return true;
      }
      else if (!name.CompareTo(fileName)){
	foundPath = path;
	return true;
      }
      
    }
  }
  
  return false;
  
}

//_______________________________________________________
void ParticleInfo::ReadBichselFile(){

  //Read the BichselFile and create the bichsel curve as a 
  //function of beta*gamma from it

  //If a pidCalibration File has been used then get the bichsel Functions for it instead
  if (pidCalibrationFile){

    pidCalibrationFile->cd("dEdxBetaGammaFits");
    for (int iParticle=0; iParticle<nDefinedParticles; iParticle++){
      for (unsigned int yIndex=0; yIndex<bichselEmpiricalFit.at(iParticle).size(); yIndex++){

	//Get the Function from the File
	TF1 *tempFunc = (TF1 *)gDirectory->Get(Form("dEdxBetaGammaFit_%s_%d",
						    GetParticleName(iParticle).Data(),yIndex));

	if (!tempFunc)
	  continue;
	
	//Rebuild the function
	bichselEmpiricalFit.at(iParticle).at(yIndex) =
	  new TF1(tempFunc->GetName(),this,&ParticleInfo::BichselEmpiricalFitPrototype,minBetaGamma,maxBetaGamma*2,7);
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetNpx(10000000);
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(0,tempFunc->GetParameter(0));
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(1,tempFunc->GetParameter(1));
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(2,tempFunc->GetParameter(2));
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(3,tempFunc->GetParameter(3));
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(4,tempFunc->GetParameter(4));
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(5,tempFunc->GetParameter(5));
	bichselEmpiricalFit.at(iParticle).at(yIndex)->SetParameter(6,tempFunc->GetParameter(6));

	delete tempFunc;
	
      }//End Loop Over rapidity bins
    }//End Loop Over Particles    
    
    return;
  }

  //Use the file for this star library version
  TString bichselFileName = TString::Format("%s_BichselCurves.root",
					    libraryVersion.Data());

  TString pathToBichselFile;

  //Search for the file in the current directory
  Bool_t fileFound = FindBichselFile(bichselFileName,TString(gSystem->pwd()),pathToBichselFile);
 
  //If the file is not found in the current directory then try searching one directory up
  if (!fileFound){
    fileFound = FindBichselFile(bichselFileName,TString(gSystem->pwd())+"/../",pathToBichselFile);
  }

  if (!fileFound){
    std::cout <<"ERROR: ParticleInfo::ReadBichselFile(): BichselFile " <<bichselFileName <<" not found!\n";
    exit (EXIT_FAILURE);
  }
  
  bichselFileName = pathToBichselFile+"/"+bichselFileName;

  //Open the bichselFile
  bichselRootFile = new TFile(bichselFileName,"READ");

  if (type == TRUNCATED){
    bichselFunctionPart[0] = (TF1*)bichselRootFile->Get("bichselTruncated_0");
    bichselFunctionPart[1] = (TF1*)bichselRootFile->Get("bichselTruncated_1");
  }
  else {
    bichselFunctionPart[0] = (TF1*)bichselRootFile->Get("bichselMostProbable_0");
    bichselFunctionPart[1] = (TF1*)bichselRootFile->Get("bichselMostProbable_1");
  }

  bichselFunction = new TF1("bichselFunction",this,&ParticleInfo::BichselStitch,minBetaGamma,maxBetaGamma,0);
  bichselFunction->SetNpx(100000);
}

//____________________________________________________________________
Double_t ParticleInfo::BichselStitch(Double_t *x, Double_t *par){

  Double_t betaGamma = x[0];

  if (par != NULL){
    std::cout <<"WARNING: ParticleInfo::BichselStitch(): Function par array is not null." <<std::endl;
  }

  if (betaGamma < midBetaGamma)
    return bichselFunctionPart[0]->Eval(betaGamma);
  
  return bichselFunctionPart[1]->Eval(betaGamma);  

}

//____________________________________________________________________
Bool_t ParticleInfo::SetEmpiricalBichselFunction(Int_t speciesInterest, Int_t yIndex){

  //When the empirical Bichsel Fits are used there is a parameterization
  //for each particle species and rapidity bin. Use this function to speficy
  //the current particle and rapidity bin of interest.
  //Returns true if successful and false otherwise

  if (!useEmpiricalFits){
    std::cout <<"WARNING - ParticleInfo::SetEmpiricalBichselFunction() - You only need to make a call to this" <<std::endl
	 <<"          function if you are using the empirical fits. No action taken." <<std::endl;
    return true;
  }

  empiricalSpeciesInterest = speciesInterest;
  empiricalRapidityIndex   = yIndex;

  bichselFunction = NULL;
  bichselFunction = bichselEmpiricalFit.at(empiricalSpeciesInterest).at(empiricalRapidityIndex);

  if (!bichselFunction)
    return false;

  return true;
}

//____________________________________________________________________
TF1 *ParticleInfo::GetBichselFunction(){

  //Returns Bichsel curve as a function of beta*gamma

  return bichselFunction;

}

//____________________________________________________________________
TF1 *ParticleInfo::GetBichselFunction(Double_t mass, Double_t charge){

  //Returns Bichsel curve as a function of momentum assuming mass
  
  //If the bichselGraph does not exist then make it
  if (!bichselFunction){
    std::cout <<"ERROR: ParticleInfo::GetBichselFunction() - BichselFuntion is NOT set. If you want to use the\n"
	 <<"       PID Functions (Bichsel/TOF) then use the constructor:\n"
	 <<"       ParticleInfo::ParticleInfo(TString starLibraryVersion, createAll=true)\n"
	 <<"       or call ParticleInfo::CreateAllDefinedParticlePidFunctions() prior to this method.\n"
	 <<"       EXITING!\n";
    exit (EXIT_FAILURE);
  }

  //Get Momentum Range based on the range of bichsel curve  
  Double_t minMomentum, maxMomentum;
  GetBichselFunction()->GetRange(minMomentum,maxMomentum);
  minMomentum = mass*minMomentum;
  maxMomentum = mass*maxMomentum;
  
  TF1 *bichselCurve = new TF1("",this,&ParticleInfo::BichselPrototype,
			      minMomentum,maxMomentum,2,
			      "ParticleInfo","BichselPrototype");
  bichselCurve->SetParameter(0,mass);
  bichselCurve->SetParameter(1,charge);
  bichselCurve->SetNpx(100000);

  return bichselCurve;
}


//___________________________________________________________________
TF1 *ParticleInfo::GetBichselFunction(Int_t speciesIndex){

  //Return the Bichsel function as a function of momentum for 
  //particle of type speciesIndex

  //If using the empirical fit you also need to specify which particle of interest and
  //its rapidity index that you want

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"ERROR: ParticleInfo::GetBichselFunction() Unkown speciesIndex!\n";
    exit (EXIT_FAILURE);
  }

  //--------------------------------------------
  //If using the STAR Library
  if (!useEmpiricalFits){

    //If it already exists return it
    if (bichselPIDFunctions[speciesIndex])
      return bichselPIDFunctions[speciesIndex];

    bichselPIDFunctions[speciesIndex] = 
      GetBichselFunction(GetParticleMass(speciesIndex),GetParticleChargeMag(speciesIndex));
    
    return bichselPIDFunctions[speciesIndex];
  }

  //--------------------------------------------
  //If Using the Empirical Fits
  if (empiricalSpeciesInterest < 0 || empiricalSpeciesInterest >= nDefinedParticles){
    std::cout <<"ERROR: ParticleInfo::GetBichselFunction() - Unkown index of empiricalSpeciesInterest!\n";
    std::cout <<"       Have you called ParticleInfo::SetEmpiricalBichselFunction(speciesInterest,rapidityIndex)?\n";
    exit (EXIT_FAILURE);
  }
  if (empiricalRapidityIndex < 0 || empiricalRapidityIndex >= (int)bichselEmpiricalFit.at(speciesIndex).size()){
    std::cout <<"ERROR: ParticleInfo::GetBichselFunction() - Invalid empiricalRapidityIndex value!\n";
    std::cout <<"       Have you called ParticleInfo::SetEmpiricalBichselFunction(speciesInterest,rapidityIndex)?\n";
    exit (EXIT_FAILURE);
  }

  return bichselEmpiricalPIDFunctions.at(empiricalSpeciesInterest).at(empiricalRapidityIndex).at(speciesIndex);

}

//___________________________________________________________________
Double_t ParticleInfo::GetParticleMass(Int_t speciesIndex){

  //Returns the particle mass of type speciesIndex in GeV.

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"WARNING - ParticleInfo::GetParticleMass(): Unknown speciesIndex! Returning NAN.\n";
    return NAN;
  }

  Double_t particleSpeciesMass[] = 
    {
      0.13957,    //Pion
      0.49367,    //Kaon
      0.93827,    //Proton
      0.00051099, //Electron
      1.87561,    //Deuteron
      2.80925,    //Triton
      2.8094,     //Helion
      3.72742,    //Alpha
      0.105658,    //Muon
      1.1157      //Lambda 
    };

  return particleSpeciesMass[speciesIndex];
}

//___________________________________________________________________
Double_t ParticleInfo::GetParticleChargeMag(Int_t speciesIndex){

  //Returns the magnitude of the charge for a particle of type speciesIndex

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"WARNING - ParticleInfo::GetParticleMass(): Unknown speciesIndex! Returning NAN.\n";
    return NAN;
  }

  Double_t particleSpeciesCharge[] = 
    {
      1.0,    //Pion
      1.0,    //Kaon
      1.0,    //Proton
      1.0,    //Electron
      1.0,    //Deuteron
      1.0,    //Triton
      2.0,    //Helion
      2.0,    //Alpha
      1.0,     //Muon
      0.0     //Lambda
    };

  return particleSpeciesCharge[speciesIndex];
}

//___________________________________________________________________
TString ParticleInfo::GetParticleName(Int_t speciesIndex,Int_t charge){

  //Returns the name of a particle species. By default (charge=0) the 
  //generic name for the species is given, otherwise Plus or Minus is
  //appended to the name depending of whether charge is positive or negative.

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"WARNING - ParticleInfo::GetParticleName():Unknown speciesIndex! Returning empty string.\n";
    return TString("");
  }

  TString chargeName;
  if (charge > 0)
    chargeName = "Plus";
  else if (charge < 0)
    chargeName = "Minus";
  else 
    chargeName = "";

  TString particleSpeciesName[] = 
    {"Pion","Kaon","Proton","Electron",
     "Deuteron","Triton","Helion","Alpha","Muon","Lambda"};

  return particleSpeciesName[speciesIndex]+chargeName;

}

//___________________________________________________________________
TString ParticleInfo::GetParticleSymbol(Int_t speciesIndex, Int_t charge){

  //Returns the root version of the latex symbol for a particle. By
  //default (charge=0) so the generic symbol is returned. If the charge
  //is set to either positive or negative the particle or antiparticle 
  //symbol will be returned

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"WARNING - ParticleInfo::GetParticleSymbol(): Unknown speciesIndex! Returning empty string.\n";
    return TString("");
  }

  TString particleSymbol[] = 
    {"#pi","K","p","e","d","t","{}^{3}He","#alpha","#mu","#Lambda"};

  if (charge == 0)
    return particleSymbol[speciesIndex];

  TString positiveParticleSymbol[] = 
    {"#pi^{+}","K^{+}","p","e^{+}","d","t","{}^{3}He","#alpha","#mu^{+}","#Lambda"};

  if (charge > 0)
    return positiveParticleSymbol[speciesIndex];

  TString negativeParticleSymbol[] =
    {"#pi^{-}","K^{-}","#bar{p}","e^{-}","#bar{d}","#bar{t}","{}^{3}#bar{He}","#bar{#alpha}","#mu^{-}","#Lambda"};
  
  return negativeParticleSymbol[speciesIndex];

}

//___________________________________________________________________
Style_t ParticleInfo::GetParticleMarker(Int_t speciesIndex, Int_t charge){

  //Returns the marker style for particle of type speciesIndex.
  //Negative particles are given closed symbols, while positive particles
  //are given open symbols. By default (charge=0) and the negative particle
  //symbol will be returned.

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"WARNING - ParticleInfo::GetParticleMarker(): Unknown speciesIndex! Using kDot.\n";
    return kDot;
  }

  Style_t negativeParticleMarker[] = 
    {kFullCircle,       //Pion
     33,                //Kaon
     kFullSquare,       //Proton
     kFullCircle,       //Electron
     34,                //Deuteron
     kFullTriangleUp,   //Triton
     kFullTriangleDown, //Helion
     kPlus,             //Alpha
     kStar,            //Muon     
     kOpenCircle};     //lambda

  if (charge <= 0)
    return negativeParticleMarker[speciesIndex];

  Style_t positiveParticleMarker[] =
    {kOpenCircle,       //Pion
     kOpenDiamond,      //Kaon
     kOpenSquare,       //Proton
     kOpenCircle,       //Electron
     kOpenCross,        //Deuteron
     kOpenTriangleUp,   //Triton
     32,                //Helion
     kPlus,             //Alpha
     kStar,            //Muon
     kOpenCircle};     //lambda

  return positiveParticleMarker[speciesIndex];
  
}

//___________________________________________________________________
Color_t ParticleInfo::GetParticleColor(Int_t speciesIndex){

  //Returns the color of the particle of type speciesIndex.

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"WARNING - ParticleInfo::GetParticleColor(): Unknown speciesIndex! Using kBlack.\n";
    return kBlack;
  }

  Color_t particleColor[] = 
    {kRed-2,           //Pion
     kBlue+2,          //Kaon
     kGreen+2,         //Proton
     kViolet-1,        //Electron
     kCyan+3,          //Deuteron
     kOrange+2,        //Triton
     kMagenta+2,       //Helion
     kOrange-2,        //Alpha
     kPink-9,         //Muon
     kOrange-2};       //lambda

  return particleColor[speciesIndex];
}

//___________________________________________________________________
void ParticleInfo::CreateAllDefinedParticlePidFunctions(){

  //Creates the Bichsel and 1/Beta Curves for all the Particles defined
  for (Int_t iParticle=0; iParticle<nDefinedParticles; iParticle++){

    //The Inverse Beta Function
    inverseBetaPIDFunctions[iParticle] = GetInverseBetaFunction(iParticle);
    inverseBetaPIDFunctions[iParticle]->SetLineColor(GetParticleColor(iParticle));
    
    if (!useEmpiricalFits){
      bichselPIDFunctions[iParticle] = GetBichselFunction(iParticle);
      bichselPIDFunctions[iParticle]->SetLineColor(GetParticleColor(iParticle));
    }
    else {
      for (unsigned int yIndex=0; yIndex<bichselEmpiricalFit.at(iParticle).size(); yIndex++){
	for (unsigned int iSubParticle=0; iSubParticle<nDefinedParticles; iSubParticle++){

	  if (!SetEmpiricalBichselFunction(iParticle,yIndex))
	    continue;
	  bichselEmpiricalPIDFunctions.at(iParticle).at(yIndex).at(iSubParticle) =
	    GetBichselFunction(GetParticleMass((int)iSubParticle),GetParticleChargeMag((int)iSubParticle));
	  
	  bichselEmpiricalPIDFunctions.at(iParticle).at(yIndex).at(iSubParticle)->
	    SetLineColor(GetParticleColor(iSubParticle));
	 
	}//End Loop Over sub particle
      }//End Loop Over rapidity bin

    }//End else
    
  }//End Loop Over iParticle

  empiricalSpeciesInterest = -1;
  empiricalRapidityIndex = -1;

  pidFunctionsCreated = true;
  
}

//___________________________________________________________________
TF1 *ParticleInfo::GetInverseBetaFunction(Double_t mass, Double_t charge){

  //Return a pointer to a TF1 of inverse beta as a function of 
  //momentum with the given mass assumption

  TF1 *inverseBetaTemp = new TF1("",this,&ParticleInfo::InverseBetaPrototype,
				 .001,10,2,
				 "ParticleInfo","InverseBetaPrototype");
  inverseBetaTemp->SetParameter(0,mass);
  inverseBetaTemp->SetParameter(1,charge);
  inverseBetaTemp->SetNpx(100000);
  
  return inverseBetaTemp;
}

//____________________________________________________________________
TF1 *ParticleInfo::GetInverseBetaFunction(Int_t speciesIndex){

  //Return the Inverse beta function as a function of momentum for 
  //this particle species

  if (speciesIndex < 0 || speciesIndex >= nDefinedParticles){
    std::cout <<"ERROR: ParticleInfo::GetInverseBetaFunction() Unkown speciesIndex!\n";
    exit (EXIT_FAILURE);
  }

  //Make Sure it doesn't already exist
  if (inverseBetaPIDFunctions[speciesIndex])
    return inverseBetaPIDFunctions[speciesIndex];

  inverseBetaPIDFunctions[speciesIndex] = 
    GetInverseBetaFunction(GetParticleMass(speciesIndex),GetParticleChargeMag(speciesIndex));

  inverseBetaPIDFunctions[speciesIndex]->SetLineColor(GetParticleColor(speciesIndex));

  return inverseBetaPIDFunctions[speciesIndex];

}

//____________________________________________________________________
Double_t ParticleInfo::ComputeZTPC(Double_t dEdxMeasured, Double_t ptMeasured,
		     Double_t pzMeasured, Int_t speciesIndex){

  //Returns the Z variable of the mesured track to the expected
  //dEdx of speciesIndex

  Double_t totalMomentum = sqrt(pow(ptMeasured,2) + pow(pzMeasured,2));
  Double_t dEdxExpected  = GetBichselFunction(speciesIndex)->Eval(totalMomentum);

  return log( (dEdxMeasured*pow(10,6)) / dEdxExpected);

}

//______________________________________________________________________
Double_t ParticleInfo::ComputeWeightedZTPC(Double_t dEdxMeasured, Double_t ptMeasured,
					   Double_t pzMeasured, Int_t speciesIndex,Double_t avgTotalMomInBin){
  

  //This is the simultaneous recentering method used by Evan Sagaline in his thesis.

  //The Likelihood function for each species is a normal gaussian (function in constructor)

  //Expected dEdx value for this speciesIndex given avgTotalMomentum
  Double_t dEdxSpeciesAvg = TMath::Log(GetBichselFunction(speciesIndex)->Eval(avgTotalMomInBin));

  //Total Measured Momentum of this Track
  Double_t totalMomentum = sqrt(pow(ptMeasured,2) + pow(pzMeasured,2));
  dEdxMeasured = TMath::Log(pow(10,6) * dEdxMeasured);
  
  //First Term
  Double_t term1(0);
  Double_t term1Num(0);
  Double_t term1Denom(0);
  for (int iParticle=0; iParticle<nDefinedParticles; iParticle++){

    //Do not consider Electrons or Muons
    if (iParticle == ELECTRON || iParticle == MUON)
      continue;

    Double_t avgMean = TMath::Log(GetBichselFunction(iParticle)->Eval(avgTotalMomInBin));
    Double_t mean  = GetBichselFunction(iParticle)->Eval(totalMomentum);
    Double_t width = canonicalTPCWidth * mean;
    mean = TMath::Log(mean);

    gausFunc->SetParameters(1.0/(width*TMath::Sqrt(2.0*TMath::Pi())),mean,width);
    
    term1Num += gausFunc->Eval(dEdxMeasured) * (dEdxMeasured + avgMean - mean);

    term1Denom += gausFunc->Eval(dEdxMeasured);
  }

  term1 = term1Num / term1Denom;

  return term1 - dEdxSpeciesAvg;

}

//____________________________________________________________________
Double_t ParticleInfo::ComputeWeightedZTOF(Double_t inverseBetaMeasured, Double_t ptMeasured,
				   Double_t pzMeasured, Int_t speciesIndex, Double_t avgTotalMomInBin){

  //This is the simultaneous recentering method used by Evan Sagaline in his thesis

  //Expected 1/beta value for a particle of this species given avgTotalMomInBin
  Double_t invBetaSpeciesAvg = GetInverseBetaFunction(speciesIndex)->Eval(avgTotalMomInBin);

  //Total Momentum of This Track
  Double_t totalMomentum = sqrt(pow(ptMeasured,2) + pow(pzMeasured,2));
  
  //First Term
  Double_t term1(0);
  Double_t term1Num(0);
  Double_t term1Denom(0);
  for (int iParticle=0; iParticle<nDefinedParticles; iParticle++){

    //Do not consider Electrons or Muons
    if (iParticle == ELECTRON || iParticle == MUON)
      continue;

    Double_t avgMean = GetInverseBetaFunction(iParticle)->Eval(avgTotalMomInBin);
    Double_t mean  = GetInverseBetaFunction(iParticle)->Eval(totalMomentum);
    Double_t width = canonicalTOFWidth * mean;
    
    gausFunc->SetParameters(1.0/(width*TMath::Sqrt(2.0*TMath::Pi())),mean,width);
    
    term1Num += gausFunc->Eval(inverseBetaMeasured) *
      (inverseBetaMeasured + avgMean - mean);

    term1Denom += gausFunc->Eval(inverseBetaMeasured);
  }

  term1 = term1Num / term1Denom;

  return term1 - invBetaSpeciesAvg;
  

}

//____________________________________________________________________
Double_t ParticleInfo::PredictZTPC(Double_t totalMomentum, Int_t speciesIndexInterest, 
				   Int_t speciesIndexConfound){

  //Returns the Z TPC variable of the particle species of interest compared to
  //the confounding species.

  Double_t dEdxInterest = GetBichselFunction(speciesIndexInterest)->Eval(totalMomentum);
  Double_t dEdxConfound = GetBichselFunction(speciesIndexConfound)->Eval(totalMomentum);

  return log(dEdxConfound / dEdxInterest);
}

//____________________________________________________________________
Double_t ParticleInfo::ConfoundMeanFitTPC(Double_t *x, Double_t *par){

  //Function that can be used to fit the Means of particles as a function of mT-m0
  Double_t mTm0 = x[0];

  Int_t pidInterest = (Int_t)par[0];
  Int_t pidConfound = (Int_t)par[1];
  Double_t rapidity = par[2];

  Double_t mass = GetParticleMass(pidInterest);
  Double_t mT   = mTm0 + mass;
  
  Double_t pZ   = mT*sinh(rapidity);
  Double_t pT   = sqrt( pow(mTm0,2) + 2.0*mass*mTm0 );
  Double_t totalMomentum = sqrt(pow(pZ,2) + pow(pT,2));
  
  return par[3] + par[4]*pow(PredictZTPC(totalMomentum,pidInterest,pidConfound),par[5]);;
}

//____________________________________________________________________
Double_t ParticleInfo::ComputeZTOF(Double_t inverseBetaMeasured, Double_t ptMeasured,
				   Double_t pzMeasured, Int_t speciesIndex){

  //Returns the Z TOF variable of the measured track to the expected 
  //Inverse beta of the species Index

  Double_t totalMomentum = sqrt(pow(ptMeasured,2) + pow(pzMeasured,2));
  Double_t inverseBetaExpected = GetInverseBetaFunction(speciesIndex)->Eval(totalMomentum);

  return inverseBetaMeasured - inverseBetaExpected;
}


//_____________________________________________________________________
Double_t ParticleInfo::PredictZTOF(Double_t totalMomentum, Int_t speciesIndexInterest, 
				   Int_t speciesIndexConfound){

  //Returns the Z TOF variable of the particle species of interest compared to
  //the confounding species
  
  Double_t inverseBetaInterest = GetInverseBetaFunction(speciesIndexInterest)->Eval(totalMomentum);
  Double_t inverseBetaConfound = GetInverseBetaFunction(speciesIndexConfound)->Eval(totalMomentum);

  return inverseBetaConfound - inverseBetaInterest;
}

//____________________________________________________________________
Double_t ParticleInfo::ConfoundMeanFitTOF(Double_t *x, Double_t *par){

  //Function that can be used to fit the Means of particles as a function of mT-m0
  Double_t mTm0 = x[0];

  Int_t pidInterest = (Int_t)par[0];
  Int_t pidConfound = (Int_t)par[1];
  Double_t rapidity = par[2];

  Double_t mass = GetParticleMass(pidInterest);
  Double_t mT   = mTm0 + mass;
  
  Double_t pZ   = mT*sinh(rapidity);
  Double_t pT   = sqrt( pow(mTm0,2) + 2.0*mass*mTm0 );
  Double_t totalMomentum = sqrt(pow(pZ,2) + pow(pT,2));
  
  return par[3] + par[4]*pow(PredictZTOF(totalMomentum,pidInterest,pidConfound),par[5]);
}

//_____________________________________________________________________
Double_t ParticleInfo::GetDistanceToSpeciesTOF(Double_t inverseBetaMeasured,
					       Double_t totalMomentum, Int_t speciesIndexInterest){

  //Return the distance between the measured inverse beta and the expected inverse beta
  //for particle of type speciesIndexInterest and at totalMomentum

  //Get the Inverse Beta Function for this speciesIndexInterest
  TF1 *temp = GetInverseBetaFunction(speciesIndexInterest);

  if (!temp){
    std::cout <<"ParticleInfo::GetDistanceToSpeciesTOF() - The Inverse Beta Function Does not exist!\n";
    std::cout <<"                                          Returning -999\n";
    return -9999;
  }
  

  return fabs(inverseBetaMeasured - temp->Eval(totalMomentum));

}

//_____________________________________________________________________
Double_t ParticleInfo::GetDistanceToSpeciesTPC(Double_t dEdxMeasured,
					       Double_t totalMomentum, Int_t speciesIndexInterest){

  //Return the distance between the measured dEdx and the expected dEdx
  //for particle of type speciesIndexInterest and at totalMomentum

  //NOTE: dEdxMeasured must be in units of KeV/cm

  //Get the Bichsel Function for this speciesIndexInterest
  TF1 *temp = GetBichselFunction(speciesIndexInterest);

  if (!temp){
    std::cout <<"ParticleInfo::GetDistanceToSpeciesTPC() - The Bichsel Function Does not exist!\n";
    std::cout <<"                                          Returning -999\n";
    return -9999;
  }

  return fabs(dEdxMeasured - temp->Eval(totalMomentum));
	      
}

//_____________________________________________________________________
Double_t ParticleInfo::GetPercentErrorTOF(Double_t inverseBetaMeasured,
					  Double_t totalMomentum, Int_t speciesIndexInterest){

  //Return the percent error between the measured inverse beta and the expected
  //inverse beta for a particle of type speciesIndexInterst with totalMomentum

  return fabs(GetDistanceToSpeciesTOF(inverseBetaMeasured,totalMomentum,speciesIndexInterest)
	      / GetInverseBetaFunction(speciesIndexInterest)->Eval(totalMomentum)) * 100.0;

}

//_____________________________________________________________________
Double_t ParticleInfo::GetPercentErrorTPC(Double_t dEdxMeasured,
					  Double_t totalMomentum, Int_t speciesIndexInterest){

  //Return the percent error between the measured dEdx and the expected
  //dEdx for a particle of type speciesIndexInterst with totalMomentum

  //NOTE: dEdxMeasured must be in units of KeV/cm
  
  return fabs(GetDistanceToSpeciesTPC(dEdxMeasured,totalMomentum,speciesIndexInterest)
	      / GetBichselFunction(speciesIndexInterest)->Eval(totalMomentum)) * 100.0;
	     
}

//______________________________________________________________________
Int_t ParticleInfo::GetNearestSpeciesTOF(Double_t inverseBetaMeasured, Double_t totalMomentum){

  //Return the particle species index which minimized the distance between the meausered
  //inverse beta and the expected inverse beta of particles with totalMomentum

  //NOTE: THIS DOES NOT CONSIDER MUONS!

  Double_t minDistance = 999999999999;
  Int_t    minSpecies  = -1;
  
  for (unsigned int iParticle=PION; iParticle<MUON; iParticle++){

    Double_t distance = GetDistanceToSpeciesTOF(inverseBetaMeasured,totalMomentum,iParticle);    
    if (distance < minDistance){
      minDistance = distance;
      minSpecies  = iParticle;
    }

  }

  return minSpecies;  
}

//______________________________________________________________________
Int_t ParticleInfo::GetNearestSpeciesTPC(Double_t dEdxMeasured, Double_t totalMomentum){

  //Return the particle species index which minimized the distance between the meausered
  //dEdx and the expected dEdx of particles with totalMomentum

  //NOTE: THIS DOES NOT CONSIDER MUONS!

  //NOTE: dEdxMeasured must be in units of KeV/cm

  Double_t minDistance = 999999999999;
  Int_t    minSpecies  = -1;
  
  for (unsigned int iParticle=PION; iParticle<MUON; iParticle++){

    Double_t distance = GetDistanceToSpeciesTPC(dEdxMeasured,totalMomentum,iParticle);    
    if (distance < minDistance){
      minDistance = distance;
      minSpecies  = iParticle;
    }

  }

  return minSpecies;  
}

//_________________________________________________________________________
Int_t ParticleInfo::GetNearestSpeciesTPCTOF(Double_t dEdxMeasured,
					    Double_t inverseBetaMeasured, Double_t totalMomentum){

  //Return the particle species index which minimizes the 2D distance between the
  //measured dEdx and expected dEdx and measured inverse Beta and expected inverse Beta.

  //NOTE: THIS DOES NOT CONSIDER MUONS!

  //NOTE: dEdxMeasured must be in units of KeV/cm
  
  Double_t minDistance = 999999999999;
  Int_t    minSpecies  = -1;
  
  for (unsigned int iParticle=PION; iParticle<MUON; iParticle++){

    Double_t tpcDistance = GetDistanceToSpeciesTPC(dEdxMeasured,totalMomentum,iParticle);
    Double_t tofDistance = GetDistanceToSpeciesTOF(inverseBetaMeasured,totalMomentum,iParticle);
    Double_t distance = sqrt(tpcDistance*tpcDistance + tofDistance*tofDistance); 
    if (distance < minDistance){
      minDistance = distance;
      minSpecies  = iParticle;
    }    
  }
  
  return minSpecies;  
  
}

//__________________________________________________________________________
Bool_t ParticleInfo::DoTPCandTOFAgreeOnNearestSpecies(Double_t dEdxMeasured, Double_t inverseBetaMeasured,
						      Double_t totalMomentum){

  //Returns true if the Both the TPC and TOF aggree that the particle is of the species
  //based on the nearest species functions

  Int_t tpcSpecies = GetNearestSpeciesTPC(dEdxMeasured,totalMomentum);
  Int_t tofSpecies = GetNearestSpeciesTOF(inverseBetaMeasured,totalMomentum);

  if (tpcSpecies == tofSpecies)
    return true;

  return false;
  
}

//__________________________________________________________________________
Bool_t ParticleInfo::DoTPCandTOFAgreeOnThisSpecies(Double_t dEdxMeasured, Double_t inverseBetaMeasured,
						   Double_t totalMomentum, Int_t speciesIndexInterest){

  //Returns true if the Both the TPC and TOF aggree that the particle is of the species
  //based on the nearest species functions

  Int_t tpcSpecies = GetNearestSpeciesTPC(dEdxMeasured,totalMomentum);
  Int_t tofSpecies = GetNearestSpeciesTOF(inverseBetaMeasured,totalMomentum);
  
  if (tpcSpecies == tofSpecies && tpcSpecies == speciesIndexInterest)
    return true;

  return false;
  
}

//___________________________________________________________________________
Bool_t ParticleInfo::IsExclusivelyThisSpeciesTPC(Double_t dEdxMeasured, Int_t speciesIndexInterest,
						 Double_t totalMomentum, Double_t percent,
						 Bool_t electron){

  //Return true if the Percent error for this track and speciesIndex is less that percent
  //AND no other particles have a percent error for this species less than percent

  //NOTE: THIS DOES NOT CONSIDER MUONS!
  
  //NOTE: dEdxMeasured must be in units of KeV/cm

  //NOTE: The optional argument electron can be used to turn off consideration of electron
  //      as a confounding species. By default electron=true so it is considered.

  //Compute the percent error of species of interst and check that
  //it is within the required percentage
  Double_t percentErrorOfInterest = GetPercentErrorTPC(dEdxMeasured,totalMomentum,speciesIndexInterest);
  if (percentErrorOfInterest > percent)
    return false;  
  
  //Compute the percent error for other species and check if any
  //percentage is less that percent
  for (unsigned int iParticle=PION; iParticle<MUON; iParticle++){

    //Skip the species of interest
    if (iParticle == (unsigned int)speciesIndexInterest)
      continue;

    //Skip the electron if requested
    if (!electron && iParticle == ELECTRON)
      continue;

    Double_t percentError = GetPercentErrorTPC(dEdxMeasured,totalMomentum,iParticle);
    if (percentError < percent)
      return false;
    
  }//End Loop Over Particles

  return true;
}

//___________________________________________________________________________
Bool_t ParticleInfo::IsExclusivelyThisSpeciesTOF(Double_t inverseBetaMeasured, Int_t speciesIndexInterest,
						 Double_t totalMomentum, Double_t percent,
						 Bool_t electron){

  //Return true if the Percent error for this track and speciesIndex is less that percent
  //AND no other particles have a percent error for this species less than percent

  //NOTE: THIS DOES NOT CONSIDER MUONS!
  
  //NOTE: The optional argument electron can be used to turn off consideration of electron
  //      as a confounding species. By default electron=true so it is considered.

  //Compute the percent error of species of interst and check that
  //it is within the required percentage
  Double_t percentErrorOfInterest = GetPercentErrorTOF(inverseBetaMeasured,totalMomentum,speciesIndexInterest);
  if (percentErrorOfInterest > percent)
    return false;  
  
  //Compute the percent error for other species and check if any
  //percentage is less that percent
  for (unsigned int iParticle=PION; iParticle<MUON; iParticle++){

    //Skip the species of interest
    if (iParticle == (unsigned int)speciesIndexInterest)
      continue;

    //Skip the electron if requested
    if (!electron && iParticle == ELECTRON)
      continue;

    Double_t percentError = GetPercentErrorTOF(inverseBetaMeasured,totalMomentum,iParticle);
    if (percentError < percent)
      return false;
    
  }//End Loop Over Particles

  return true;
}

//___________________________________________________________________________
Bool_t ParticleInfo::IsExclusivelyThisSpeciesTPCTOF(Double_t dEdxMeasured, Double_t inverseBetaMeasured,
						    Int_t speciesIndexInterest, Double_t totalMomentum,
						    Double_t tpcPercent, Double_t tofPercent,
						    Bool_t electron){

  //Return true if the Percent error for this track and speciesIndex is less that percent
  //AND no other particles have a percent error for this species less than percent

  //NOTE: THIS DOES NOT CONSIDER MUONS!
  
  //NOTE: dEdxMeasured must be in units of KeV/cm

  //NOTE: The optional argument electron can be used to turn off consideration of electron
  //      as a confounding species. By default electron=true so it is considered.

  Bool_t exclusiveTPC = IsExclusivelyThisSpeciesTPC(dEdxMeasured,speciesIndexInterest,totalMomentum,
						    tpcPercent,electron);

  if (!exclusiveTPC)
    return false;

  Bool_t exclusiveTOF = IsExclusivelyThisSpeciesTOF(inverseBetaMeasured,speciesIndexInterest,totalMomentum,
						    tofPercent,electron);
  if (!exclusiveTOF)
    return false;
  
  return true;
}
