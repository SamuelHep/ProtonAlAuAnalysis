#include <iostream>
#include <vector>

#include "TH1D.h"

#include "inc/ReturnStructs.h"

using namespace std;

void FindCentralityBinCuts(NegBinomialSearchResults *results){

  //This function computes the particle multiplicity cuts for each of the
  //centrality bins defined below.
  TH1D *simHisto = results->bestFitHisto;
  
  //The Centrality cut vector (gets returned)
  std::vector <double> centralityCuts;

  //Centrality Bin Definitions
  const Int_t nCentBins(16);
  Double_t centBinDefs[nCentBins] =
    {.80, .75, .70, .65, .60, .55, .50, .45, .40,
     .35, .30, .25, .20, .15, .10, .05};
  //Top % Centrality (must bin in order of {high%,...,low%}
  
  //Centrality Bin Definitions
  results->nCentralityBins = nCentBins;
  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++)
    results->centralityBinDefinitions.push_back(centBinDefs[iCentBin]);

  //The Total Integral of the simulation histogram
  Double_t totalIntegral = simHisto->Integral();
  
  //For Each Centrality bin find the particle multiplicity
  //correspoing to the % centrality
  for (Int_t iCentBin=0; iCentBin<nCentBins; iCentBin++){
    Double_t tempIntegral(0);
    Int_t histoBin(1);
    while (tempIntegral / totalIntegral < 1.0-results->centralityBinDefinitions[iCentBin]
	   && histoBin >0){

      tempIntegral = simHisto->Integral(1,histoBin);
      histoBin++;

    }//End Loop Over Histo Bins
    
    results->centralityBinCuts.push_back(simHisto->GetBinLowEdge(histoBin));
    
  }//End Loop Over Centrality Bins

}
