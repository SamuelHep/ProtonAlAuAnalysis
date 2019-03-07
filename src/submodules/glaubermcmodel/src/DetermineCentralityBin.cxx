#include <iostream>
#include <vector>

#include <DetermineCentralityBin.h>

using namespace std;

//_________________________________________________________________________________
Int_t DetermineCentralityBin16(std::vector<double> centralityCuts, Double_t centralityVar){

  //This Function loops over the centrality bins and determines which bin this event
  //belongs to based on the value of its centralityVar

  //Check that centralityVar isn't less than the minimum centrality cut
  if (centralityVar < centralityCuts.at(0))
    return -1;

  //Loop Over the centrality cuts
  for (int iCentBin=(int)centralityCuts.size()-1; iCentBin>=0; iCentBin--){

    if (centralityVar > centralityCuts.at(iCentBin))
      return iCentBin;
    
  }
  
  //If for some reason we don't find a centrality bin then return something silly
  return -1;

}

//___________________________________________________________________________________
Int_t DetermineCentralityBin9(std::vector<double> centralityCuts, Double_t centralityVar){

  //This function simply converts the 16 bin centrality down to 9 bins
  
  Int_t centralityBin16 = DetermineCentralityBin16(centralityCuts, centralityVar);

  //Consider Special Cases

  //If this is the highest cent bin
  if (centralityBin16 == 15)
    return 8;
  else if (centralityBin16 == 14)
    return 7;
  else if (centralityBin16 < 0) //This is the case where a bin was not found
    return centralityBin16;

  //The General Case
  return (Int_t)(0.5 * centralityBin16);
  

}

//__________________________________________________________________________________
Int_t DetermineCentralityBin(std::vector<double> centralityCuts, Double_t centralityVar, Int_t nCentBins){

  //This Function simply allows the user to choose from the two above

  if (nCentBins == 16)
    return DetermineCentralityBin16(centralityCuts, centralityVar);

  //Otherwise
  return DetermineCentralityBin9(centralityCuts, centralityVar);

}
