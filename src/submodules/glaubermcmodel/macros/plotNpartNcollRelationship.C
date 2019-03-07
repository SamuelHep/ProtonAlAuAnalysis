void plotNpartNcollRelationship(TString filename){

  TFile *file = new TFile(filename);
  TTree *tree = (TTree *)file->Get("GlauberTree");

  TCanvas *canvas = new TCanvas("canvas","canvas",20,20,800,600);
  canvas->SetLogz();
  
  TH2D *nPartnCollHisto = new TH2D("nPartnCollHisto","nPartnCollHisto",200,0,1200,225,0,450);
  tree->Draw("nParticipants:nBinaryCollisions>>nPartnCollHisto","","COLZ");

  TGraphErrors *nPartnCollGraph = new TGraphErrors();
  TH1D *htemp;
  for (Int_t i=1; i<nPartnCollHisto->GetNbinsX(); i++){

    htemp = nPartnCollHisto->ProjectionY("",i,i);

    if (htemp->GetEntries() < 2)
      continue;
    
    nPartnCollGraph->SetPoint(nPartnCollGraph->GetN(),
			      nPartnCollHisto->GetXaxis()->GetBinCenter(i),
			      htemp->GetMean());
    nPartnCollGraph->SetPointError(nPartnCollGraph->GetN()-1,
				   0.5*nPartnCollHisto->GetXaxis()->GetBinWidth(i),
				   htemp->GetRMS());

  }

  nPartnCollGraph->SetMarkerStyle(7);
  nPartnCollGraph->SetFillColor(kRed);
  nPartnCollGraph->SetFillStyle(3001);
  nPartnCollGraph->Draw("3");
  
}
