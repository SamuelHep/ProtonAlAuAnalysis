
#include "TH1D.h"
#include "TGraph2D.h"
#include <vector>
#include <iostream>

#include "ReturnStructs.h"

using namespace std;

//___________________________________________
void Delete(NegBinomialSearchResults *result){

  //This Function is responsible for deleting
  //the object passed into it
  delete result->bestFitHisto;
  delete result->chi2Graph;

}

//____________________________________________
ClassImp(RefMultCorrInfo);

RefMultCorrInfo::RefMultCorrInfo(){

}

RefMultCorrInfo::~RefMultCorrInfo(){

}
