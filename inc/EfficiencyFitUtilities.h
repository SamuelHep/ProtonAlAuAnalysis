//Some functions useful for fitting efficiency curves
Double_t GetOneSigmaEfficiencyDeviateAt(TGraphAsymmErrors *gr, Double_t xVal);
Int_t CountNAbove(Int_t n, Double_t *vec, Double_t threshold);
Int_t CountNBelow(Int_t n, Double_t *vec, Double_t threshold);
Double_t FindFractionalCut(Int_t n, Double_t *vec, Double_t desiredFrac,
			   bool forward=true, const int nTrials=1000);
TMatrixD GetRandomDeviateVector(TF1 *func, TFitResultPtr fitResult);
TGraphAsymmErrors *ConstructErrorFuncConfidenceBand(TF1 *fitFunc, TFitResultPtr *fitResult, Double_t percenttile);
