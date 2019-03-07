//Various colors and marker styles for energies, particles, centralities, ect.

Double_t GetEnergy(int iEnergy){

  if (iEnergy < 0 || iEnergy >=7){
    cout <<"StyleSettings::GetEnergy() - Energy Index not defined. Returning 0\n";
    return 0.0;
  }
  
  Double_t energies[7] = {7.7,11.5,14.5,19.6,27.0,39.0,62.4};

  return energies[iEnergy];

}

TString GetEventConfig(int iConfig){

  if (iConfig < 0 || iConfig >= 3){
    cout <<"StyleSettings::GetEventConfig() - Config index not defined. Returning empty string.\n";
    return TString("");
  }
  
  TString configs[3] = {"ColliderCenter","ColliderPosY","ColliderNegY"};

  return configs[iConfig];
  
}

Style_t GetCentralityMarker(){

  

}

Color_t GetCentralityColor(int centBin){

  if (centBin==8)
    return kRed;
  else if (centBin==7)
    return kOrange+7;
  else if (centBin==6)
    return kOrange-3;
  else if (centBin==5)
    return kGreen+2;
  else if (centBin==4)
    return kAzure+7;
  else if (centBin==3)
    return kBlue+2;
  else if (centBin==2)
    return kViolet-3;
  else if (centBin==1)
    return kMagenta+2;
  else if (centBin==0)
    return kMagenta-2;

  
}

Color_t GetEnergyColor(Double_t energy){

  if (fabs(energy) - 7.7 < 2.)
    return kRed;
  else if (fabs(energy - 11.5) < 2.)
    return kOrange+7;
  else if (fabs(energy - 14.5) < 2.)
    return kOrange-3;
  else if (fabs(energy - 19.6) < 2.)
    return kGreen+2;
  else if (fabs(energy - 27.) < 2.)
    return kAzure+7;
  else if (fabs(energy - 39.) < 2.)
    return kBlue+2;
  else if (fabs(energy - 62.4) < 2.)
    return kViolet-3;
  
}
