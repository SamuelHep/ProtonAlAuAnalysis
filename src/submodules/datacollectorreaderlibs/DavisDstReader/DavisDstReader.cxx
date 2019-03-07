#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TMath.h>
#include <TVirtualTreePlayer.h>

#include "DavisDstReader.h"
#include "../TrackInfo/TrackInfo.h"
#include "../PrimaryVertexInfo/PrimaryVertexInfo.h"
#include "../EventInfo/EventInfo.h"

ClassImp(DavisDstReader);

//____________________________________________________________
DavisDstReader::DavisDstReader(){

  //NULL constructor;

  //Set Pointers to NULL
  file = NULL;
  tree = NULL;
  eventList = NULL;
  eventInfo = NULL;
  vertexArray = NULL;
  trackArray = NULL;
}

//____________________________________________________________
DavisDstReader::DavisDstReader(TString fileNameWithFullPath){

  //Set Pointers to NULL
  file = NULL;
  tree = NULL;
  eventList = NULL;
  eventInfo = NULL;
  vertexArray = NULL;
  trackArray = NULL;
  
  ConstructTree(fileNameWithFullPath);
  
  totalEntries = tree->GetEntries();
  currentEntry = 0;

}

//____________________________________________________________
DavisDstReader::DavisDstReader(TString fileNameWithFullPath, Double_t minVz, Double_t maxVz){

  //This constructore uses an event list to reduce the number
  //of entries that are read from the tree during loops by
  //cutting on min and max z vertex locations

  //Set Pointers to NULL
  file = NULL;
  tree = NULL;
  eventList = NULL;
  eventInfo = NULL;
  vertexArray = NULL;
  trackArray = NULL;
  
  ConstructTree(fileNameWithFullPath,minVz,maxVz);

  tree->Draw(">>eventList",Form("VertexInfo.zVertex>=%g&&VertexInfo.zVertex<=%g",minVz,maxVz));
  eventList = (TEventList *)gDirectory->Get("eventList");
  tree->SetEventList(eventList);
  totalEntries = eventList->GetN();
  currentEntry = eventList->GetEntry(0);
    
}

//____________________________________________________________
DavisDstReader::~DavisDstReader(){

  tree->Delete();

}

//____________________________________________________________
void DavisDstReader::ConstructTree(TString fileNameWithFullPath, Double_t minVz, Double_t maxVz){

  //Tree and Chain are used synonomously
  tree = new TChain("DataTree");

  //If the fileNameWithFullPathVariable contains one of the extensions:
  //".list",".txt",or ".input" the the file is a filelist list so loop over
  //it and add the files to the chain
  if (fileNameWithFullPath.Contains(".list") ||
      fileNameWithFullPath.Contains(".txt")  ||
      fileNameWithFullPath.Contains(".input")){

    ifstream fileList(fileNameWithFullPath);
    if (!fileList.is_open()){
      fputs("ERROR:: DavisDstReader::ConstructTree() - File list not open! Exiting\n",stderr);
      exit (EXIT_FAILURE);
    }

    string fileName;
    while (fileList.good()){
      getline(fileList, fileName);
      if (fileName.length() != 0 && fileName.find(".root") != std::string::npos)
	cout <<"Adding: " <<fileName <<endl;
	tree->Add(fileName.c_str());
    };

    fileList.close();

  }
  else {
    tree->Add(fileNameWithFullPath.Data());
  }

  if (!tree){
    fputs("ERROR: DavisDstReader::ConstructTree() - Requested tree is not found!",stderr);
    exit (EXIT_FAILURE);
  }
  
  //Set the Branch Addresses
  tree->SetBranchAddress("EventInfo",&eventInfo);
  tree->SetBranchAddress("VertexInfo",&vertexArray);
  tree->SetBranchAddress("TrackInfo",&trackArray);

  //Tree Properties
  tree->SetEstimate(tree->GetEntries()+1);
  tree->SetCacheSize(100000000); //100 MBytes
  tree->AddBranchToCache("*",true);
  
}

//____________________________________________________________
void DavisDstReader::SetBranchStatus(TString branchName, Int_t active){

  tree->SetBranchStatus(branchName,active);

}

//____________________________________________________________
TTree *DavisDstReader::GetTree(){

  //Return a pointer to the tree
  //Returns NULL if pointer has not been set
  return (TTree *)tree;

}

//____________________________________________________________
EventInfo *DavisDstReader::GetEntry(Long64_t iEntry){

  if (iEntry < 0){
    fputs("ERROR: DavisDstReader::GetEntry() - entry number must be positive!\n",stderr);
    exit (EXIT_FAILURE);
  }

  if (iEntry >= totalEntries){
    fputs("ERROR: DavisDstReader::GetEntry() - entry number is greater than total entries!\n",stderr);
    exit (EXIT_FAILURE);
  }

  //Get the Entry from the Tree
  Int_t bytesRead(0);
  long entryToGet(0);
  if (!eventList)
    entryToGet = iEntry;
  else{
    entryToGet = eventList->GetEntry(iEntry);
  }

  bytesRead = tree->GetEntry(entryToGet);
  currentEntry = entryToGet;
  
  //SUCCESS
  //If the number of bytes read from the tree is greater than
  //zero then it was a successful read
  if (bytesRead > 0){
    eventInfo->SetVertexArrayPtr(vertexArray);
    eventInfo->SetTrackArrayPtr(trackArray);
    return eventInfo;
  }

  //FAILURE
  fputs("WARNING: DavisDstReader::GetEntry() - no bytes were read from the tree.\n",stderr);
  fputs("         Returning NULL pointer to event!.\n",stderr);
  return NULL;

}

//_______________________________________________________________
EventInfo *DavisDstReader::FindEvent(Long64_t runNumber, Long64_t eventNumber){

  //Return a pointer to the event matching runNumber and eventNumber 
  //NOTE: If, for whatever reason, there are multiple entries in the tree
  //      with the same combination of run and event numbers, this will
  //      return the first instance.
  
  tree->Draw("Entry$",Form("EventInfo.runNumber==%lld&&EventInfo.eventNumber==%lld",
			   runNumber,eventNumber),"goff");

  Int_t nSelected = tree->GetSelectedRows();
  if (nSelected <= 0){
    return NULL;
  }
  else if (nSelected > 1){
    cout <<"WARNING: DavisDstReader::FindEvent() - More than one entry matches these criteria!\n"
	 <<"         Only the first matching event will be returned. This is weird!\n";
  }
  
  Double_t *list = tree->GetV1();
  tree->GetEntry((long)list[0]);
  currentEntry = (long)list[0];
  eventInfo->SetVertexArrayPtr(vertexArray);
  eventInfo->SetTrackArrayPtr(trackArray);

  return eventInfo;
  
}

//__________________________________________________________________
std::vector<long> DavisDstReader::GetEventsWithRunNumber(Long64_t runNumber){

  //Return a vector of indices corresponding to events with matching
  //runNumbers. WARNING: No other cuts are applied to these events!

  tree->Draw("Entry$",Form("EventInfo.runNumber==%lld",runNumber),"goff");
  
  Long64_t nSelected = tree->GetSelectedRows();
  Double_t *list = tree->GetV1();

  //By defualt return retuns an array of Doubles and since we want a vec
  //of longs we can't just use the constructor
  std::vector<long> vec(nSelected,0);
  for (Long64_t i=0; i<nSelected; i++){
    vec.at(i)=list[i];
  }
  
  return vec;
}

//____________________________________________________________________
std::vector<long> DavisDstReader::GetEventsWithTriggerID(unsigned int triggerID){

  //Return a vector of indices corresponding to events whose trigger
  //vectors contain triggerID. WARNING: No other cuts are applied

  tree->Draw("Entry$",Form("EventInfo.triggerIDs==%u",triggerID),"goff");

  Int_t nSelected = tree->GetSelectedRows();
  Double_t *list = tree->GetV1();

  //By defualt return retuns an array of Doubles and since we want a vec
  //of longs we can't just use the constructor
  std::vector<long> vec(nSelected,0);
  for (int i=0; i<nSelected; i++){
    vec.at(i)=list[i];
  }

  return vec;

}

