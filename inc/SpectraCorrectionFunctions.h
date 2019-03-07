//These functions apply the various corrections to the spectra.
//They are used in src/analysis/correctSpectra.cxx

class CorrectionCurves;

TGraphErrors *ApplyCorrectionChain(TGraphErrors *rawSpectrum, CorrectionCurves *corrCurve,
				   Double_t minTPCEfficiency, Bool_t doSystematics=false);
TGraphErrors *ApplyInterBinTransverseMassLocation(TGraphErrors *spectrum, TF1 *interBinAverageCurve,
						  TGraphErrors *interBinAverageCurveSys=NULL);

TGraphErrors *ApplyFiniteBinWidthCorrection(TGraphErrors *spectrum, TF1 *energyLossCurve, Double_t pidMass,
					    TGraphErrors *energyLossCurveSys=NULL);
TGraphErrors *ApplyEnergyLossCorrection(TGraphErrors *spectrum, TF1 *energyLossCurve, Double_t pidMass,
					TGraphErrors *energyLossCurveSys=NULL);
TGraphErrors *ApplyEfficiencyCorrection(TGraphErrors *spectrum, TF1 *efficiencyCurve, Double_t minTPCEfficiency,
					TGraphErrors *efficiencyCurveConf=NULL, TF1 *efficiencyCurveSys=NULL,
					TGraphErrors *efficiencyCurveSysConf=NULL);
TGraphErrors *ApplyKnockoutProtonCorrection(TGraphErrors *spectrum, TF1 *backgroundCurve,
					    Bool_t doSystematics=false);
TGraphErrors *ApplyMuonBackgroundCorrection(TGraphErrors *spectrum, TF1 *muonBackgroundCurveDefault,
					    TF1 *muonBackgroundCurveAlt=NULL,
					    TGraphErrors *muonBackgroundCurveDefaultSys=NULL,
					    TGraphErrors *muonBackgroundCurveAltSys=NULL);
TGraphErrors *ApplyFeedDownBackgroundCorrection(TGraphErrors *spectrum, TF1 *feedDownBackgroundCurveDefault,
						TF1 *feedDownBackgroundCurveAlt=NULL,
						TGraphErrors *feedDownBackgroundCurveDefaultSys=NULL,
						TGraphErrors *feedDownBackgroundCurveAltSys=NULL);
//TGraphErrors *MakeSystematicErrorSpectrum(TGraphErrors *nominalSpectrum, TGraphErrors *widthSpectrum,
//					  const int nSpectra, std::vector<TGraphErrors *> *sysSpectra);
TGraphErrors *MakeSystematicErrorSpectrum(TGraphErrors *nominalSpectrum,
					  const int nSpectra, std::vector<TGraphErrors *> *sysSpectra);
TGraphErrors *MakeTotalErrorSpectrum(TGraphErrors *nominalSpectrum, TGraphErrors *systematicSpectrum);

TF1 *GetOneSigmaEfficiencyDeviateFunction(TF1 *func);
std::vector<double> GetCorrectedBinWidthVector(TGraphErrors *spectrum, TF1 *interBinMeanCurve, TF1 *energyLossCurve,
					       Double_t pidMass, TF1 *efficiencyCurve, Double_t minTPCEfficiency, TF1 *feedDownCurve);

std::vector<int> CleanSpectrum(TGraphErrors *gr, int pid);
