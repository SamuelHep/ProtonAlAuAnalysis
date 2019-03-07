//Functions are that potentially useful in all sripts

#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TVirtualFitter.h>

//Kinematics
Double_t GetRapidityRangeLow(Int_t rapidityIndex);
Double_t GetRapidityRangeHigh(Int_t rapidityIndex); 
Double_t GetRapidityRangeCenter(Int_t rapidityIndex);
Int_t GetRapidityIndex(Double_t rapidity);
Double_t GetmTm0RangeLow(Int_t mTm0Index);
Double_t GetmTm0RangeHigh(Int_t mTm0Index);
Double_t GetmTm0RangeCenter(Int_t mTm0Index);
Int_t GetmTm0Index(Double_t mTm0);
Double_t ComputemTm0(Double_t pT, Double_t mass);
Double_t ConvertmTm0ToPt(Double_t mTm0, Double_t mass);
Double_t ComputepTotal(Double_t pT, Double_t pZ);
Double_t ComputepZ(Double_t mT, Double_t rapidity);

Int_t GetCentralityRangeLow(int centBin);
Int_t GetCentralityRangeHigh(int centBin);

//General
int FindElementWithMinDistance(int n, double *arr, double val);

//TGraph Operations
void TGraphChop(TGraphErrors *graph, Double_t threshold, Bool_t below=true);
TGraphErrors *TGraphScale(TGraphErrors *graph, Double_t scaleFactor, Bool_t returnNew=true);
Int_t TGraphFindPoint(TGraph *graph, Double_t xVal);

//TMultiGraph Operations
int MultiGraphGetN(TMultiGraph *graph);

//Confidence Intervals
TGraphErrors *GetConfidenceIntervalOfFit(TF1 *fitFunc,Double_t sigma=.683);

//dEdx Parameterizations
Double_t dEdxBetaGammaFit(Double_t *x, Double_t *par);

//Moving Average Graphs
TGraphErrors *MovingAverageGraph(TGraphErrors *g,TString direction="forward",Double_t range=-999);
TGraphErrors *MovingAverageGraphFB(TGraphErrors *g,Double_t range=-999);
Double_t GetErrorAtValue(TGraphErrors *g, Double_t xVal, Double_t epsilon=.025/2.0);
TGraphErrors *GetMovingAverageErrorBand(TGraphErrors *mvAvgGraph);

//Eta-Phi Mappings - Maps a track's eta and phi to tpc sector index
void DefineTPCSectorBoundaries();
Int_t GetSectorID(Double_t trackEta, Double_t trackPhi);

//TGraph Symmetrization (usefull for dNdy plots)
TGraphErrors *SymmeterizedNdyGraph(TGraphErrors *gr);
TGraphErrors *SystematicdNdyGraph(TGraphErrors *gr);

//Geant ID Symbol
TString GetParticleSymbolUsingGeantID(int geantID);

TString GetSpectraAxisTitle(){
  return TString("#frac{1}{N_{Evt}}#times#frac{1}{2#pim_{T}}#times#frac{d^{2}N}{dm_{T}dy} (GeV/c^{2})^{-2}");
}

#endif
