//This is a helper class to load, access, and
//keep track of all the correction functions
//that are applied to the spectra.

class ParticleInfo;
class TFile;
class TF1;
class TGraphErrors;

class CorrectionCurves {

private:
  int pid;
  int charge;
  int iCentBin;
  int yIndex;

  ParticleInfo *particleInfo = NULL;

  TFile *corrFile = NULL;
  
  TF1 *interBinMeanCurve = NULL;
  TF1 *energyLossCurve = NULL;
  TF1 *efficiencyCurve = NULL;
  TF1 *efficiencyCurveSys = NULL;
  TF1 *knockoutProtonCurve = NULL;
  TF1 *muonBackgroundCurvePowerLaw = NULL;
  TF1 *muonBackgroundCurveExpo = NULL;
  TF1 *feedDownBackgroundCurveExpo = NULL;
  TF1 *feedDownBackgroundCurvePowerLaw = NULL;
  
  TGraphErrors *interBinMeanCurveConf = NULL;
  TGraphErrors *energyLossCurveConf = NULL;
  TGraphErrors *efficiencyCurveConf = NULL;
  TGraphErrors *efficiencyCurveSysConf = NULL;
  TGraphErrors *knockoutProtonCurveConf = NULL;
  TGraphErrors *muonBackgroundCurvePowerLawConf = NULL;
  TGraphErrors *muonBackgroundCurveExpoConf = NULL;
  TGraphErrors *feedDownBackgroundCurveExpoConf = NULL;
  TGraphErrors *feedDownBackgroundCurvePowerLawConf;
  
  void LoadCurves();

public:
  CorrectionCurves(TFile *CORRECTIONFILE, int PID, int CHARGE, int CENTINDEX, int YINDEX);
  virtual ~CorrectionCurves();

  int GetPid(){return pid;}
  int GetCharge(){return charge;}
  int GetCentIndex(){return iCentBin;}
  int GetRapidityIndex(){return yIndex;}
  
  TF1 *GetInterBinMeanCurve(){return interBinMeanCurve;}
  TF1 *GetEnergyLossCurve(){return energyLossCurve;}
  TF1 *GetEfficiencyCurve(){return efficiencyCurve;}
  TF1 *GetEfficiencyCurveSys(){return efficiencyCurveSys;}
  TF1 *GetKnockoutProtonCurve(){return knockoutProtonCurve;}
  TF1 *GetMuonBackgroundCurvePowerLaw(){return muonBackgroundCurvePowerLaw;}
  TF1 *GetMuonBackgroundCurveExpo(){return muonBackgroundCurveExpo;}
  TF1 *GetFeedDownBackgroundCurvePowerLaw(){return feedDownBackgroundCurvePowerLaw;}
  TF1 *GetFeedDownBackgroundCurveExpo(){return feedDownBackgroundCurveExpo;}
  
  TGraphErrors *GetInterBinMeanCurveConf(){return interBinMeanCurveConf;}
  TGraphErrors *GetEnergyLossCurveConf(){return energyLossCurveConf;}
  TGraphErrors *GetEfficiencyCurveConf(){return efficiencyCurveConf;}
  TGraphErrors *GetEfficiencyCurveSysConf(){return efficiencyCurveSysConf;}
  TGraphErrors *GetKnockoutProtonCurveConf(){return knockoutProtonCurveConf;}
  TGraphErrors *GetMuonBackgroundCurvePowerLawConf(){return muonBackgroundCurvePowerLawConf;}
  TGraphErrors *GetMuonBackgroundCurveExpoConf(){return muonBackgroundCurveExpoConf;}
  TGraphErrors *GetFeedDownBackgroundCurvePowerLawConf(){return feedDownBackgroundCurvePowerLawConf;}
  TGraphErrors *GetFeedDownBackgroundCurveExpoConf(){return feedDownBackgroundCurveExpoConf;}


};
