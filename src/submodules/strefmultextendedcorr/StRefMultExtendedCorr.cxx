#include <iostream>
#include <stdlib.h>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include "StRefMultExtendedCorr.h"

ClassImp(RefMultCorrInfo);

//______________________________________________________
RefMultCorrInfo::RefMultCorrInfo(){
  //Nothing to do here
}

//______________________________________________________
RefMultCorrInfo::~RefMultCorrInfo(){
  //Nothing to do here
}


ClassImp(StRefMultExtendedCorr);

//______________________________________________________
StRefMultExtendedCorr::StRefMultExtendedCorr(){

  //Default Constructor
  isInitialized = false;
  cmEnergy = -1;
  corrFile = NULL;
  corrTree = NULL;
  treeEntry = NULL;

}

//______________________________________________________
Bool_t StRefMultExtendedCorr::FindCorrectionFile(TString fileName, TString path, TString &foundPath){

  //Recursively seartch the directory tree until the file
  //with the name fileName is found and then return its path

  //Get the Contents of path
  TSystemDirectory dir(path.Data(),path.Data());
  TList *list = dir.GetListOfFiles();

  //Look through the contents of pwd until the file is found
  if (list){

    TSystemFile *file;
    TString name;
    TIter next(list);

    while ((file=(TSystemFile*)next())){
      name = file->GetName();

      if (name.Contains(".") && !name.Contains(".root"))
	continue;

      if (file->IsDirectory()){
	if (FindCorrectionFile(fileName,path+"/"+name,foundPath))
	  return true;
      }

      else if (!name.CompareTo(fileName)){
	foundPath = path;
	return true;
      }

    }//End While Loop

  }

  //If the file was not found
  return false;

}

//______________________________________________________
void StRefMultExtendedCorr::Initialize(double energy){

  //Set all the data members

  //Only Do this if it hasn't already been done
  if (isInitialized)
    return;

  cmEnergy = energy;

  TString correctionFile;
  
  if (energy > 4.0 && energy < 8.0) // Fix for now
    //    if (energy > 7.0 && energy < 8.0)
    correctionFile = "AuAu07_CorrectionInfo.root";
  else if (energy > 10.0 && energy < 12.0)
    correctionFile = "AuAu11_CorrectionInfo.root";
  else if (energy > 13.0 && energy < 16.0)
    correctionFile = "AuAu14_CorrectionInfo.root";
  else if (energy > 18.0 && energy < 20.0)
    correctionFile = "AuAu19_CorrectionInfo.root";
  else if (energy > 26.0 && energy < 29.0)
    correctionFile = "AuAu27_CorrectionInfo.root";
  else if (energy > 38.0 && energy < 40.0)
    correctionFile = "AuAu39_CorrectionInfo.root";
  else if (energy > 61.0 && energy < 63.0)
    correctionFile = "AuAu62_CorrectionInfo.root";

  TString pathToCorrectionFile;

  //Search for the Correction File in the current directory
  Bool_t fileFound = FindCorrectionFile(correctionFile,TString(gSystem->pwd()),pathToCorrectionFile);

  //If the file was not found in the current directory then try searching one directory higher
  if (!fileFound){
    fileFound = FindCorrectionFile(correctionFile,TString(gSystem->pwd())+"/../",pathToCorrectionFile);
  }

  //If the file was still not found then try searching from one more directory higher
  if (!fileFound){
    fileFound = FindCorrectionFile(correctionFile,TString(gSystem->pwd())+"/../../",pathToCorrectionFile);
  }
 
  
  //If the file was not found then throw an error
  if (!fileFound){
    std::cout <<"ERROR: StRefMultExtendedCorr::Inititialize - Correction File " <<correctionFile <<" not found!\n";
    exit (EXIT_FAILURE);
  }

  correctionFile = pathToCorrectionFile+"/"+correctionFile;

  //Open the File Containing the Correction Info
  corrFile = new TFile(correctionFile,"READ");

  //Get the Correction Tree 
  corrTree = (TTree *)corrFile->Get("CorrectionTree");
  corrTree->SetBranchAddress("CorrectionInfo",&treeEntry);

  //corrTree->Print();

  //Get the Number of Triggers which is the same as the
  //number of entries in the correction tree
  nTriggers = corrTree->GetEntries();

  //Create the Arrays of the Correction Graphs
  refMultCorrGraph = new TGraphErrors **[nTriggers];

  //Create the Array for the functions that convert 
  //RefMultPosY and RefMultNegY to RefMult
  convertToRefMultFunc = new TF1 **[nTriggers];

  //Get All the Correction Graphs
  int startTrig(0);

  //Get the Bad Run List
  badRunList = (TVectorD *)corrFile->Get("BadRunList");

  //Loop Over the TriggerIDs and get the correction/conversion objects
  for (int iTrig=startTrig; iTrig<nTriggers; iTrig++){
    
    refMultCorrGraph[iTrig] = new TGraphErrors *[3];
    convertToRefMultFunc[iTrig] = new TF1 *[3];

    corrTree->GetEntry(iTrig);

    //Get the RefMult Correction Graphs from the Correction File
    refMultCorrGraph[iTrig][0] = 
      (TGraphErrors *)corrFile->Get(TString::Format("refMultCorrectionGraph_%ld",(long)treeEntry->triggerID));
    refMultCorrGraph[iTrig][1] = 
      (TGraphErrors *)corrFile->Get(TString::Format("refMultPosYCorrectionGraph_%ld",(long)treeEntry->triggerID));
    refMultCorrGraph[iTrig][2] = 
      (TGraphErrors *)corrFile->Get(TString::Format("refMultNegYCorrectionGraph_%ld",(long)treeEntry->triggerID));

    //Get the Conversion Function from the Correction File
    convertToRefMultFunc[iTrig][0] = NULL; //THIS SHOULD NEVER BE NEEDED
    convertToRefMultFunc[iTrig][1] = 
      (TF1 *)corrFile->Get(TString::Format("refMultPosYToRefMultConversion_%ld",(long)treeEntry->triggerID));
    convertToRefMultFunc[iTrig][2] = 
      (TF1 *)corrFile->Get(TString::Format("refMultNegYToRefMultConversion_%ld",(long)treeEntry->triggerID));
    

  }//End loop to get Correction Graphs and Conversion Functions

  //Change the Initialization State
  isInitialized=true;

}

//______________________________________________________
StRefMultExtendedCorr::~StRefMultExtendedCorr(){
  //Destructor
}

//______________________________________________________
int StRefMultExtendedCorr::GetTriggerIndex(long triggerID){

  int startTrig(0);

  for (int iTrig = startTrig; iTrig<nTriggers; iTrig++){
    corrTree->GetEntry(iTrig);
    
    if (treeEntry->triggerID == triggerID)
      return iTrig;
  }

  //If a match to the trigger is not found then report error
  std::cout <<"FATAL ERROR: NO MATCH TO TRIGGER ID FOUND... " <<triggerID <<" EXITING!" <<std::endl;

  return -1;

}

//________________________________________________________
int StRefMultExtendedCorr::GetRefMultType(double zVertex){

  //The Centrality Variable that will be used to characterize the
  //event is dependent on the zVertex Location
  if (zVertex >=-95 && zVertex <= 95)
    return 0; //RefMult
  else if (zVertex <-95)
    return 1; //RefMultPosY
  else if (zVertex > 95)
    return 2; //RefMultNegY

  return -1;
}


//________________________________________________________
Bool_t StRefMultExtendedCorr::CanProcessEvent(long triggerID, double zVertex){

  //Get the Trigger Index
  int triggerIndex = GetTriggerIndex(triggerID);

  //Is the triggerID of the event a trigger there is info about?
  if (triggerIndex < 0){
    std::cout <<"WARNING!!!!!!!! StRefMultExtendedCorr: Can not find this event's triggerID....SKIPPING EVENT!\n";
    return false;
  }
  
  //Is the zVertex of the Event within the range processed of the 
  //correction factor graphs?
  Int_t nPoints = refMultCorrGraph[triggerIndex][GetRefMultType(zVertex)]->GetN();
  Double_t *xVals = refMultCorrGraph[triggerIndex][GetRefMultType(zVertex)]->GetX();
  double minRange = TMath::MinElement(nPoints, xVals);
  double maxRange = TMath::MaxElement(nPoints, xVals);
  if (zVertex > maxRange || zVertex < minRange){
    std::cout <<"WARNING!!!!!!!! StRefMultExtendedCorr: Event outside of analyzed zRange...SKIPPING EVENT!\n ";
    return false;
  }

  //If the event doesn't fail then it can be processed
  return true;
}


//________________________________________________________
double StRefMultExtendedCorr::GetRefMultCorrection(long triggerID, double zVertex, int refMultType){

  if (!CanProcessEvent(triggerID,zVertex))
    return 1.0;

  //If the User didn't specify type
  if (refMultType < 0)
    refMultType = GetRefMultType(zVertex);

  //Returns the correction factor for the refMultType 
  return refMultCorrGraph[GetTriggerIndex(triggerID)][refMultType]->Eval(zVertex);

}

//_______________________________________________________
double StRefMultExtendedCorr::GetCorrectedRefMult(long triggerID, double zVertex, double refMult, int refMultType){

  //WARNING!!!:  The refMult argument must be either the true refMult,
  //refMultPosY, or refMultNegY.  THIS IS THE USER'S RESPONSIBILITY!

  if (!CanProcessEvent(triggerID,zVertex))
    return refMult;

  //If the User didn't specify type
  if (refMultType < 0)
    refMultType = GetRefMultType(zVertex);

  return refMult / GetRefMultCorrection(triggerID,zVertex,refMultType);

}

//_______________________________________________________
double StRefMultExtendedCorr::ConvertToRefMult(long triggerID, double zVertex, double refMult){

  //Converts RefMultPosY and RefMultNegY to RefMult
  
  if (!CanProcessEvent(triggerID,zVertex))
    return 1.0;

  //Get the RefMult Type of refMultPosYNegY
  int refMultType = GetRefMultType(zVertex);

  //If the refMultType is RefMult (0) then just return refMult (it doesn't need conversion)
  if (refMultType == 0)
    return refMult;

  //Otherwise Use the conversion functions to convert RefMultPosY or RefMultNegY to RefMult
  return convertToRefMultFunc[GetTriggerIndex(triggerID)][refMultType]->GetX(refMult);

}

//_______________________________________________________
int StRefMultExtendedCorr::GetCentralityBin16(long triggerID, double zVertex, double refMult, int refMultType){

  //WARNING!!! The refMult argument must be either the true refMult,
  //refMultPosY, or refMultNegY. THIS IS THE USER'S RESPONSIBILITY!

  if (!CanProcessEvent(triggerID,zVertex))
    return -1;

  //If the User didn't specify a type
  if (refMultType < 0)
    refMultType = GetRefMultType(zVertex);

  //Get the Corrected RefMult
  double correctedRefMult = GetCorrectedRefMult(triggerID,zVertex,refMult,refMultType);
  
  //Convert RefMultPosY or RefMultNegY to RefMult
  if (refMultType != 0)
    correctedRefMult = ConvertToRefMult(triggerID,zVertex,correctedRefMult);
  
  int triggerIndex = GetTriggerIndex(triggerID);

  corrTree->GetEntry(triggerIndex);
  
  for (Int_t iCentBin=0; iCentBin<16; iCentBin++){
    
    if (iCentBin+1 < 16 &&
	correctedRefMult >= treeEntry->refMultCentBinCuts[iCentBin] &&
	correctedRefMult < treeEntry->refMultCentBinCuts[iCentBin +1]){
      return iCentBin;
    }
    else if (iCentBin+1 == 16 &&
	     correctedRefMult >= treeEntry->refMultCentBinCuts[iCentBin])
      return iCentBin;
       
  }

  //cout <<"WARNING: CENTRALITY BIN NOT FOUND!" <<endl;
  return -1;

}


//_______________________________________________________
int StRefMultExtendedCorr::GetCentralityBin9(long triggerID, double zVertex, double refMult, int refMultType){

  //WARNING!!! The refMult argument must be either the true refMult, 
  //refMultPosY, or refMultNegY.  THIS IS THE USER'S RESPONSIBILITY!

  //If the User didn't specify a type
  if (refMultType < 0)
    refMultType = GetRefMultType(zVertex);

  int centBin16 = GetCentralityBin16(triggerID,zVertex,refMult,refMultType);

  //First Consider Special Cases
  if (centBin16 == 15)
    return 8;
  else if (centBin16 == 14)
    return 7;
  else if (centBin16 < 0)
    return centBin16; //Centrality Bin Not Defined


  return (int)(0.5 * centBin16);

}

//______________________________________________________
Bool_t StRefMultExtendedCorr::IsBadRun(long runNumber){

  //Search the badRunList to see if if runNumber is contained within
  //it. It it is return true (this is bad run), otherwise return
  //false (it is a good run).
  
  for (Int_t i=0; i<badRunList->NonZeros(); i++){
    std::cout <<(long)runNumber <<" " <<(long)( (*badRunList) )[i] <<std::endl;
    if ((long)( (*badRunList) )[i] == (long)runNumber)
      return true;

  }

  return false;
}
