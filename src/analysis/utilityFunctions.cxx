//Functions that are potentially useful in all scripts

#include <iostream>

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TList.h>

#include "globalDefinitions.h"
#include "utilityFunctions.h"

//___________________________________________________________
Double_t GetRapidityRangeLow(Int_t rapidityIndex){

  return rapidityMin + (rapidityIndex * rapidityBinWidth);
}

//___________________________________________________________
Double_t GetRapidityRangeHigh(Int_t rapidityIndex){

  return GetRapidityRangeLow(rapidityIndex)+rapidityBinWidth;
}

//___________________________________________________________
Double_t GetRapidityRangeCenter(Int_t rapidityIndex){

  return (GetRapidityRangeLow(rapidityIndex) + GetRapidityRangeHigh(rapidityIndex)) / 2.0;
}

//___________________________________________________________
Int_t GetRapidityIndex(Double_t rapidity){

  //Make sure the value of rapidity is in the [min,max] range
  if (rapidity < rapidityMin || rapidity > rapidityMax)
    return -1;

  //Loop through the rapidity bins and return the index
  //of the bin that this value of rapidity falls into
  for (Int_t index=0; index<nRapidityBins; index++){
    
    if (rapidity >= GetRapidityRangeLow(index) &&
	rapidity < GetRapidityRangeHigh(index))
      return index;
  }

  return -1;
}

//___________________________________________________________
Double_t GetmTm0RangeLow(Int_t mTm0Index){

  return mTm0Min + (mTm0Index * mTm0BinWidth);
}

//___________________________________________________________
Double_t GetmTm0RangeHigh(Int_t mTm0Index){

  return GetmTm0RangeLow(mTm0Index) + mTm0BinWidth;

}

//___________________________________________________________
Double_t GetmTm0RangeCenter(Int_t mTm0Index){

  return (GetmTm0RangeLow(mTm0Index) + GetmTm0RangeHigh(mTm0Index)) / 2.0;
}

//___________________________________________________________
Int_t GetmTm0Index(Double_t mTm0){

  //Make sure the value of mTm0 is in the [min,max] range
  if (mTm0 < mTm0Min || mTm0 > mTm0Max)
    return -1;

  //Loop through the mTm0 bins and return the index
  //of the bin that this value of mTm0 falls into
  for (Int_t index=0; index<nmTm0Bins; index++){

    if (mTm0 >= GetmTm0RangeLow(index) &&
	mTm0 < GetmTm0RangeHigh(index))
      return index;
  }

  return -1;
}

//___________________________________________________________
Int_t GetCentralityRangeLow(int centBin){

  if (centBin < 0 || centBin > 8){
    cout <<"WARNING - UtilityFunctions::GetCentralityRangeLow() - Centrality bin not defined.\n";
    return -1;
  }
  
  int centPercents[] = {80,70,60,50,40,30,20,10,5};

  return centPercents[centBin];

}

//___________________________________________________________
Int_t GetCentralityRangeHigh(int centBin){

  if (centBin < 0 || centBin > 8){
    cout <<"WARNING - UtilityFunctions::GetCentralityRangeHigh() - Centrality bin not defined.\n";
    return -1;
  }
    
  int centPercents[] = {70,60,50,40,30,20,10,5,0};

  return centPercents[centBin];
}


//___________________________________________________________
Double_t ComputemTm0(Double_t pT, Double_t mass){

  return sqrt(pow(pT,2) + pow(mass,2)) - mass;
}

//___________________________________________________________
Double_t ConvertmTm0ToPt(Double_t mTm0, Double_t mass){

  return sqrt( pow(mTm0,2) + 2.0*mass*mTm0 );
}

//___________________________________________________________
Double_t ComputepTotal(Double_t pT, Double_t pZ){

  return sqrt(pow(pT,2) + pow(pZ,2));
}

//___________________________________________________________
Double_t ComputepZ(Double_t mT, Double_t rapidity){

  return mT*sinh(rapidity);
}

//___________________________________________________________
void TGraphChop(TGraphErrors *graph, Double_t threshold, Bool_t below){

  //Remove points from graph that are below (or above) the
  //the x value of threshold

  for (Int_t iPoint=graph->GetN()-1; iPoint>=0; iPoint--){

    //If the user wants to remove points above the threshold
    if (!below && graph->GetX()[iPoint] > threshold){
      graph->RemovePoint(iPoint);
    }

    else if (below && graph->GetX()[iPoint] < threshold) {
      graph->RemovePoint(iPoint);
    }

  }

}

//___________________________________________________________
TGraphErrors *TGraphScale(TGraphErrors *graph, Double_t scaleFactor, Bool_t returnNew){

  //Returns a new TGraph scaled by scaleFactor
  //By default returnNew==true and this function returns
  //a pointer to a new TGraphErrors object

  TGraphErrors *scaledGraph = NULL;
  if (returnNew){
    scaledGraph = new TGraphErrors();
    scaledGraph->SetName(graph->GetName());
    scaledGraph->SetTitle(graph->GetTitle());
    scaledGraph->SetMarkerStyle(graph->GetMarkerStyle());
    scaledGraph->SetMarkerColor(graph->GetMarkerColor());
  }
  else
    scaledGraph = graph;
  
  Int_t nPoints = graph->GetN();
  for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
    
    scaledGraph->SetPoint(iPoint,
			  graph->GetX()[iPoint],
			  graph->GetY()[iPoint]*scaleFactor);
    scaledGraph->SetPointError(iPoint,
			       graph->GetEX()[iPoint],
			       graph->GetEY()[iPoint]*scaleFactor);
    
  }
  
  return scaledGraph;
}

//_____________________________________________________________________
Int_t TGraphFindPoint(TGraph *graph, Double_t xVal){

  //Return the index of the point in graph whose x value is
  //nearest to xVal

  //Make sure the graph has points
  if (graph->GetN() < 1){
    cout <<"WARNING: UtilityFunctions::TGraphFindPoint - Graph has no points!" <<endl;
    return -1;
  }
  
  std::vector<double> diff(graph->GetN(),99999);
  for (int i=0; i<graph->GetN(); i++){
    diff.at(i) = fabs(graph->GetX()[i] - xVal);
  }

  return TMath::LocMin(diff.size(),&diff.at(0));

}

//_____________________________________________________________________
int FindElementWithMinDistance(int n, double *arr, double val){

  //Return the index of the element of the array which is
  //nearest the value of val
  std::vector<double> diffs(n);
  for (int i=0; i<n; i++){
    diffs.at(i) = fabs(arr[i] - val);    
  }

  return TMath::LocMin(n,&diffs.at(0));
  
}


//_____________________________________________________________________
int MultiGraphGetN(TMultiGraph *graph){

  //Return the number of points in the multi-graph
  //This should really be included in root
  int nPoints(0);
  int nGraphs = graph->GetListOfGraphs()->GetSize();
  TGraph *g;

  for (int iGraph=0; iGraph<nGraphs; iGraph++){
    g = (TGraph *)((TList *)graph->GetListOfGraphs())->At(iGraph);
    nPoints += g->GetN();
  }
 
  return nPoints;
}


//_____________________________________________________________________
TGraphErrors *GetConfidenceIntervalOfFit(TF1 *fitFunc, Double_t sigma){

  //Get the Confidence interval resulting from a fit
  //involving fitFunc. NOTE: This must be called IMMEDIDIATELY
  //after the fit.
  
  const int npx = fitFunc->GetNpx();
  Double_t fitMin(0), fitMax(0);
  fitFunc->GetRange(fitMin,fitMax);
  Double_t stepSize =  fabs(fitMax-fitMin)/(double)npx;

  TGraphErrors *confGraph = new TGraphErrors(npx);
  confGraph->SetName(Form("%s_Conf",fitFunc->GetName()));
  confGraph->SetFillColor(fitFunc->GetLineColor()-3);
  confGraph->SetFillStyle(3002);

  for (int i=0; i<npx; i++){
    confGraph->SetPoint(i,fitMin+i*stepSize,0);
  }

  ((TVirtualFitter::GetFitter())->GetConfidenceIntervals(confGraph,sigma));
  
  return confGraph;
  
}


//_____________________________________________________________________
Double_t dEdxBetaGammaFit(Double_t *x,  Double_t *par){

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

//___________________________________________________________________
TGraphErrors *MovingAverageGraph(TGraphErrors *g,TString direction,Double_t range){

  //Return a new TGraphErrors containing the moving average of graph g. The moving
  //average can be computed either direction (forward,backward). In addition, a range
  //can be specified over which the running average will be computed with respect
  //to each point.

  TGraphErrors *avgGraph = new TGraphErrors();

  Bool_t forward = direction.Contains("forward");
  int min = 0;
  int max = g->GetN() - 1;

  for (int i=(forward ? min:max); (forward ? (i<=max):(i>=min)); (forward ? i++:i--)){

    if (avgGraph->GetN() == 0){
      avgGraph->SetPoint(avgGraph->GetN(),g->GetX()[i],g->GetY()[i]);
      continue;
    }


    Double_t xVal = g->GetX()[i];
    Int_t innerMin = min;
    Int_t innerMax = i;
    if (range != -999){

      Double_t tempMin = xVal -range;
      Double_t tempMax = xVal + range;

      if (tempMin < g->GetX()[min]){
	innerMin = min;
      }
      else
	innerMin = TMath::BinarySearch(g->GetN(),g->GetX(),tempMin);

      if (tempMax > g->GetX()[max])
	innerMax = max;
      else
	innerMax = TMath::BinarySearch(g->GetN(),g->GetX(),tempMax);

    }

    std::vector<double> vals;
    if (forward){
      for (int j=innerMin; j<=innerMax; j++){
	vals.push_back(g->GetY()[j]);
      }
    }
    else {
      for (int j=innerMax; j>=innerMin; j--){
	vals.push_back(g->GetY()[j]);
      }
    }

    Double_t avg = TMath::Mean(vals.size(),&vals.at(0));
    Double_t rms = TMath::RMS(vals.size(),&vals.at(0));

    avgGraph->SetPoint(avgGraph->GetN(),g->GetX()[i],avg);
    avgGraph->SetPointError(avgGraph->GetN()-1,g->GetEX()[i],rms/sqrt(vals.size()));

  }

  avgGraph->SetLineWidth(4);
  if (forward)
    avgGraph->SetLineColor(kRed);
  else
    avgGraph->SetLineColor(kBlack);

  return avgGraph;
}

//__________________________________________________________________
TGraphErrors *MovingAverageGraphFB(TGraphErrors *g,Double_t range){

  //Return a new TGraphErrors object which will contain the average
  //of the forward and backward moving averages of the graphs created
  //in MovingAverageGraph() above. The errors are the quad sum of the two graph's errors

  TGraphErrors *gForward = MovingAverageGraph(g,"forward",range);
  TGraphErrors *gBackward = MovingAverageGraph(g,"backward",range);

  TGraphErrors *gCombined = new TGraphErrors();
  for (int i=0; i<gForward->GetN(); i++){

    int iB = gBackward->GetN()-1-i;

    Double_t avg = (gForward->GetY()[i] + gBackward->GetY()[iB]) / 2.0;

    Double_t err = sqrt(pow(gForward->GetEY()[i],2) + pow(gBackward->GetEY()[iB],2));

    gCombined->SetPoint(gCombined->GetN(),gForward->GetX()[i],avg);
    gCombined->SetPointError(gCombined->GetN()-1,gForward->GetEX()[i],err);

  }

  gCombined->SetLineColor(g->GetMarkerColor());
  gCombined->SetFillColor(g->GetMarkerColor());
  gCombined->SetFillStyle(3002);

  delete gForward;
  delete gBackward;

  return gCombined;

}

//__________________________________________________________________
Double_t GetErrorAtValue(TGraphErrors *g, Double_t xVal, Double_t epsilon){

  //Return the Error at the point corresponding to xVal.
  //If no point exists between [xVal-epsilon,xVal+epsilon],
  //then use a linear interpolation for the error

  if (!g){
    cout <<"WARNING! - UtilityFunctions::GetErrorAtValue() - graph does not exist! Returning 0!\n";
    return 0.0;
  }
  
    
  Int_t iPoint = TMath::BinarySearch(g->GetN(),g->GetX(),xVal);

  if (fabs(g->GetX()[iPoint] - xVal) < epsilon){
    return g->GetEY()[iPoint];
  }

  Int_t nextPoint(0);
  if (g->GetX()[iPoint] < xVal)
    nextPoint = iPoint +1;
  else
    nextPoint = iPoint -1;

  if (nextPoint < 0 || nextPoint > g->GetN()-1)
    return g->GetEY()[iPoint];

  Double_t avgErr = (g->GetEY()[iPoint] + g->GetEY()[nextPoint]) / 2.0;

  return avgErr;
}

//_________________________________________________________________
TGraphErrors *GetMovingAverageErrorBand(TGraphErrors *mvAvgGraph){

  TGraphErrors *errorBandGraph = new TGraphErrors();

  if (!mvAvgGraph)
    return errorBandGraph;
  
  if (mvAvgGraph->GetN() == 0)
    return errorBandGraph;
  
  int nSteps = mvAvgGraph->GetN()*4;
  double min = mvAvgGraph->GetX()[0]-mvAvgGraph->GetEX()[0];
  double max = mvAvgGraph->GetX()[mvAvgGraph->GetN()-1]+mvAvgGraph->GetEX()[mvAvgGraph->GetN()-1];
  double stepSize = fabs(max-min)/(double)nSteps;
  for (int i=0; i<nSteps; i++){
    double xVal = min + i*stepSize;
    errorBandGraph->SetPoint(errorBandGraph->GetN(),xVal,mvAvgGraph->Eval(xVal));
    errorBandGraph->SetPointError(errorBandGraph->GetN()-1,stepSize/2.0,GetErrorAtValue(mvAvgGraph,xVal));
  }

  errorBandGraph->SetFillStyle(3001);
  errorBandGraph->SetFillColor(mvAvgGraph->GetMarkerColor());
  
  return errorBandGraph;
}


//Eta-Phi Mappings - Maps a track's eta and phi to tpc sector index
double sectorStartBoundary[24];
double sectorStopBoundary[24];
bool sectorsAreDefined=false;
void DefineTPCSectorBoundaries(){

  //Each TPC Sector covers 0.16666*Pi or 30 degrees
  //Phi = 0 along the x Axis and points to the center
  //of sector 3 in if Eta>0 and to the center of sector 21 if Eta<0
  const int nSectors(24);
  const double sectorAngle((2.0*TMath::Pi())/((double)nSectors/2.0)); //Angles are in Radians
  const double halfSectorAngle = sectorAngle/2.0;

  //Sector Boundaries: Sector ID is 1+arrayIndex for PosEta and arrayIndex
  //double sectorStartBoundary[nSectors];
  //double sectorStopBoundary[nSectors];

  //**** Assign Sector Boundaries Positive Eta****
  //Sector 3 //Points Along the Positive X Axis
  sectorStartBoundary[2] =  0.0 - halfSectorAngle;                
  sectorStopBoundary[2]  = sectorStartBoundary[2]+sectorAngle;
  
  //Sector 2
  sectorStartBoundary[1] = sectorStopBoundary[2];
  sectorStopBoundary[1]  = sectorStartBoundary[1]+sectorAngle;
  
  //Sector 1
  sectorStartBoundary[0] = sectorStopBoundary[1];
  sectorStopBoundary[0]  = sectorStartBoundary[0]+sectorAngle;
  
  //Sector 12
  sectorStartBoundary[11] = sectorStopBoundary[0];
  sectorStopBoundary[11]  = sectorStartBoundary[11]+sectorAngle;
  
  //Sector 11
  sectorStartBoundary[10] = sectorStopBoundary[11];
  sectorStopBoundary[10]  = sectorStartBoundary[10]+sectorAngle;
  
  //Sector 10
  sectorStartBoundary[9] = sectorStopBoundary[10];
  sectorStopBoundary[9]  = sectorStartBoundary[9]+sectorAngle;
  
  //Sector 9 is a special exception

  //Sector 4
  sectorStartBoundary[3] = sectorStartBoundary[2];
  sectorStopBoundary[3]  = sectorStartBoundary[3]-sectorAngle;
  
  //Sector 5
  sectorStartBoundary[4] = sectorStopBoundary[3];
  sectorStopBoundary[4]  = sectorStartBoundary[4]-sectorAngle;
  
  //Sector 6
  sectorStartBoundary[5] = sectorStopBoundary[4];
  sectorStopBoundary[5]  = sectorStartBoundary[5]-sectorAngle;
  
  //Sector 7
  sectorStartBoundary[6] = sectorStopBoundary[5];
  sectorStopBoundary[6]  = sectorStartBoundary[6]-sectorAngle;
  
  //Sector 8
  sectorStartBoundary[7] = sectorStopBoundary[6];
  sectorStopBoundary[7]  = sectorStartBoundary[7]-sectorAngle;

  //**** Assign Sector Boundaries Negative Eta****
  //Sector 21 //Points Along the Positive X Axis
  sectorStartBoundary[20] =  0.0 - halfSectorAngle;
  sectorStopBoundary[20]  = sectorStartBoundary[2]+sectorAngle;

  //Sector 22
  sectorStartBoundary[21] = sectorStopBoundary[20];
  sectorStopBoundary[21]  = sectorStartBoundary[21]+sectorAngle;

  //Sector 23
  sectorStartBoundary[22] = sectorStopBoundary[21];
  sectorStopBoundary[22]  = sectorStartBoundary[22]+sectorAngle;

  //Sector 24
  sectorStartBoundary[23] = sectorStopBoundary[22];
  sectorStopBoundary[23]  = sectorStartBoundary[23]+sectorAngle;

  //Sector 13
  sectorStartBoundary[12] = sectorStopBoundary[23];
  sectorStopBoundary[12]  = sectorStartBoundary[12]+sectorAngle;

  //Sector 14
  sectorStartBoundary[13] = sectorStopBoundary[12];
  sectorStopBoundary[13]  = sectorStartBoundary[13]+sectorAngle;

  //Sector 15 is a special Exception

  //Sector 20
  sectorStartBoundary[19] = sectorStartBoundary[20];
  sectorStopBoundary[19]  = sectorStartBoundary[19]-sectorAngle;

  //Sector 19
  sectorStartBoundary[18] = sectorStopBoundary[19];
  sectorStopBoundary[18]  = sectorStartBoundary[18]-sectorAngle;

  //Sector 18
  sectorStartBoundary[17] = sectorStopBoundary[18];
  sectorStopBoundary[17]  = sectorStartBoundary[17]-sectorAngle;

  //Sector 17
  sectorStartBoundary[16] = sectorStopBoundary[17];
  sectorStopBoundary[16]  = sectorStartBoundary[16]-sectorAngle;

  //Sector 16
  sectorStartBoundary[15] = sectorStopBoundary[16];
  sectorStopBoundary[15]  = sectorStartBoundary[15]-sectorAngle;
  
  sectorsAreDefined = true;
  
}

//________________________________________________________________
Int_t GetSectorID(Double_t trackEta, Double_t trackPhi){

  //Make sure the TPC Sectors are defined
  if (sectorsAreDefined == false)
    DefineTPCSectorBoundaries();

  //________________________________________________________________
  //Return the sector index depending on positive or negative eta
  //________________________________________________________________
  
  //Positive Eta
  if (trackEta >= 0){
    
    //Sector 3 Exception
    if (fabs(trackPhi) < sectorStopBoundary[2])
      return 3;
    
    //Sector 9 Exception
    if (fabs(trackPhi) > sectorStopBoundary[9])
      return 9;
    
    for (int iSector=0; iSector<12; iSector++){
      
      //Skip Sectors 3 and 9
      if (iSector == 2 || iSector == 8)
	continue;
      
      if (trackPhi/fabs(trackPhi) > 0){
	if (trackPhi > sectorStartBoundary[iSector] && trackPhi <= sectorStopBoundary[iSector])
	  return iSector+1;
      }
      else {
	if (trackPhi <= sectorStartBoundary[iSector] && trackPhi > sectorStopBoundary[iSector])
	  return iSector+1;
      }
      
    }
    
  }//End Positive Eta

  //Negative Eta
  if (trackEta < 0){
    
    //Sector 21 Exception
    if (fabs(trackPhi) < sectorStopBoundary[20])
      return 21;
    
    //Sector 15 Exception
    if (fabs(trackPhi) > sectorStopBoundary[13])
      return 15;
    
    for (int iSector=12; iSector<24; iSector++){
      
      //Skip Sectors 21 and 15
      if (iSector == 14 || iSector == 20)
	continue;
      
      if (trackPhi/fabs(trackPhi) > 0){
	if (trackPhi > sectorStartBoundary[iSector] && trackPhi <= sectorStopBoundary[iSector])
	  return iSector+1;
      }
      else {
	if (trackPhi <= sectorStartBoundary[iSector] && trackPhi > sectorStopBoundary[iSector])
	  return iSector+1;
      }
      
    }
    
  }//End Negative Eta
  
  //Sector Not found
  return -1;
  
}

//__________________________________________________________
TGraphErrors *SymmeterizedNdyGraph(TGraphErrors *gr){

  //Create a new TGraphErrors containing points
  //that have been averaged across symValue. The
  //points are symmeterized by taking the average. 

  TGraphErrors *symGraph = new TGraphErrors(gr->GetN());

  //Loop Over the points of the input graph and do the symmeterization
  //skipping the symPoint
  for (int iPoint=0; iPoint<((gr->GetN()-1)/2.0); iPoint++){

    //Get the Points to be Symmeterized
    Double_t x1 = gr->GetX()[iPoint];
    Double_t y1 = gr->GetY()[iPoint];
    Double_t ex1 = gr->GetEX()[iPoint];
    Double_t ey1 = gr->GetEY()[iPoint];
    
    Double_t x2 = gr->GetX()[gr->GetN()-1-iPoint];
    Double_t y2 = gr->GetY()[gr->GetN()-1-iPoint];
    Double_t ex2 = gr->GetEX()[gr->GetN()-1-iPoint];
    Double_t ey2 = gr->GetEY()[gr->GetN()-1-iPoint];
    
    //Take the Average
    Double_t avg    = (y1+y2)/2.0;
    Double_t avgErr = sqrt(pow(ey1,2)+pow(ey2,2))/2.0; //(fabs(y1-y2)/2.0)/sqrt(2.0);
    
    symGraph->SetPoint(iPoint,x1,avg);
    symGraph->SetPointError(iPoint,ex1,avgErr);

    symGraph->SetPoint(gr->GetN()-1-iPoint,x2,avg);
    symGraph->SetPointError(gr->GetN()-1-iPoint,ex2,avgErr);
    
  }

  //Now do the mid-rapidity point (the error on the mid rapidity point is
  //defined to be the same as the error on the preceeding point
  int midYPoint = ((gr->GetN()-1)/2.0);
  symGraph->SetPoint(midYPoint,gr->GetX()[midYPoint],gr->GetY()[midYPoint]);
  symGraph->SetPointError(midYPoint,gr->GetEX()[midYPoint],symGraph->GetEY()[midYPoint-1]);
			  
  return symGraph;

}

//_____________________________________________________________
TGraphErrors *SystematicdNdyGraph(TGraphErrors *gr){

  //Computes the Systematic Error Associated with
  //the left/right asymmetry of the dNdy distributions

  TGraphErrors *sysGraph = new TGraphErrors(gr->GetN());

  //Loop Over the points of the input graph and do the symmeterization
  //skipping the symPoint
  for (int iPoint=0; iPoint<((gr->GetN()-1)/2.0); iPoint++){

    //Get the Points to be Symmeterized
    Double_t x1 = gr->GetX()[iPoint];
    Double_t y1 = gr->GetY()[iPoint];
    Double_t ex1 = gr->GetEX()[iPoint];
    Double_t ey1 = gr->GetEY()[iPoint];

    Double_t x2 = gr->GetX()[gr->GetN()-1-iPoint];
    Double_t y2 = gr->GetY()[gr->GetN()-1-iPoint];
    Double_t ex2 = gr->GetEX()[gr->GetN()-1-iPoint];
    Double_t ey2 = gr->GetEY()[gr->GetN()-1-iPoint];

    //Take the Average
    Double_t vals[2] = {y1,y2};
    Double_t avg    = TMath::Mean(2,vals);
    Double_t sysErr = TMath::RMS(2,vals)/2.0;

    sysGraph->SetPoint(iPoint,x1,avg);
    sysGraph->SetPointError(iPoint,rapidityBinWidth/4.0,sysErr);

    sysGraph->SetPoint(gr->GetN()-1-iPoint,x2,avg);
    sysGraph->SetPointError(gr->GetN()-1-iPoint,rapidityBinWidth/4.0,sysErr);

  }

  //Now do the mid-rapidity point (the error on the mid rapidity point is
  //defined to be the same as the error on the preceeding point
  int midYPoint = ((gr->GetN()-1)/2.0);
  sysGraph->SetPoint(midYPoint,gr->GetX()[midYPoint],gr->GetY()[midYPoint]);
  sysGraph->SetPointError(midYPoint,rapidityBinWidth/4.0,sysGraph->GetEY()[midYPoint-1]);

  return sysGraph;
  


}


//_____________________________________________________________
TString GetParticleSymbolUsingGeantID(int geantID){

  //Return the particle symbol belonging to geantID
  //NOTE: Geant only defines particles with indices
  //      from [1,50]. We reserve the index 0 for all others.
  //      See here: http://www.star.bnl.gov/public/comp/simu/gstar/Manual/particle_id.html

  //Check to make sure the value of geantID is in the range
  if (geantID < 0 || geantID > 50)
    return TString("No Match");

  //Defind the Symbol Array
  TString pidSymbol[51] = {
    "NA", //Not Applicable
    "#gamma", //gamma
    "e^{+}", //positron
    "e^{-}", //electron
    "#nu", //neutrino
    "#mu^{+}", //muon+
    "#mu^{-}", //muon-
    "#pi^{0}", //pion 0
    "#pi^{+}", //pion +
    "#pi^{-}", //pion -
    "K^{0}_{L}", //Kaon 0 Long
    "K^{+}", //Kaon +
    "K^{-}", //Kaon -
    "n", //Neutron
    "p", //Proton
    "#bar{p}", //AntiProton
    "K^{0}_{S}", //Kaon 0 Short
    "#eta", //Eta
    "#Lambda", //Lambda
    "#Sigma^{+}", //Sigma +
    "#Sigma^{0}", //Sigma 0
    "#Sigma^{-}", //Sigma -
    "#Xi^{0}", //Xi 0
    "#Xi^{-}", //Xi -
    "#Omega", //Omega -
    "#bar{n}", //Antineutron
    "#bar{#Lambda}", //Antilambda
    "#bar{#Sigma^{-}}", //Antisigma -
    "#bar{#Sigma^{0}}", //Antisigma 0
    "#bar{#Sigma^{0}}", //Antisigma +
    "#bar{#Xi^{0}}", //Anti Xi 0
    "#bar{#Xi^{+}}", //Anti Xi +
    "#bar{#Omega^{+}}", //Anti Omega +
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "d", //Deuteron
    "t", //Triton
    "#alpha", //Alpha
    "Geantino", //Geantino
    "{}^{3}He", //He3
    "Cherenkov", //Cherenkov    
  };

  return pidSymbol[geantID];
  
}
