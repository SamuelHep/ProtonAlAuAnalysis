//Header for the files contained in src/spectraUtilities/SpectraFitUtilities.cxx

//Compute the Total dNdy
std::pair<double,double> ComputeTotaldNdy(SpectraClass *spectrum, TF1 *lowmTFit, TF1 *highmTFit);

//Vary the fit lengths_____________________________________________
Double_t VaryFitLengthUpper(SpectraClass *spectrum, TF1 *fitFunc, TGraphErrors *chi2Graph=NULL);
Double_t VaryFitLengthLower(SpectraClass *spectrum, TF1 *fitFunc, TGraphErrors *chi2Graph=NULL);

//Errors
Double_t GetPercentErrorFromFixedSlope(SpectraClass *spectrum, TF1 *fitFunc,
				       TGraphErrors *slopeConfInterval, Double_t rapidity,
				       Double_t minFitRange, Double_t maxFitRange);
