void PlotGlauberResultsSingleSpecies(){

  gSystem->Load("../bin/GlauberClass_cxx");

  //Glauber File to read
  TFile *file = new TFile("../data/Glauber_197_197_42mb_WoodsSaxon.root");

  TTree *glauberTree = (TTree *)file->Get("GlauberTree");
  GlauberClass *glauberEvent = 0;
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  
  gStyle->SetOptStat(0);

  //Canvases
  TCanvas *canvas1 = new TCanvas("canvas1","canvas1",20,20,800,600);
  TCanvas *canvas2 = new TCanvas("canvas2","canvas2",20,20,800,600);
  TCanvas *canvas3 = new TCanvas("canvas3","canvas3",20,20,800,600);
  TCanvas *canvas4 = new TCanvas("canvas4","canvas4",20,20,800,600);

  canvas1->SetTopMargin(.05);
  canvas1->SetRightMargin(.05);
  canvas2->SetTopMargin(.05);
  canvas2->SetRightMargin(.05);
  canvas3->SetTopMargin(.05);
  canvas3->SetRightMargin(.05);
  canvas4->SetTopMargin(.05);

  canvas1->SetLogy();
  canvas2->SetLogy();
  canvas4->SetLogz();

  canvas1->SetTicks();
  canvas2->SetTicks();
  canvas3->SetTicks();
  canvas4->SetTicks();
  
  //Histograms
  TH1D *nPartCentralHisto = new TH1D("nPartCentralHisto","",225,0,450);
  TH1D *nPartMinBiasHisto = new TH1D("nPartMinBiasHisto","",225,0,450);
  TH1D *nCollCentralHisto = new TH1D("nCollCentralHisto","",300,0,1500);
  TH1D *nCollMinBiasHisto = new TH1D("nCollMinBiasHisto","",300,0,1500);
  TH1D *impactParamCentralHisto = new TH1D("impactParamCentralHisto","",60,0,20);
  TH1D *impactParamMinBiasHisto = new TH1D("impactParamMinBiasHisto","",60,0,20);
  TH2F *mScalingFactorHisto = new TH2F("mScalingFactorHisto","",225,0,450,400,0,400);

  nPartMinBiasHisto->SetTitle(";N_{part};dN/dN_{part}");
  nCollMinBiasHisto->SetTitle(";N_{coll};dN/dN_{coll}");
  impactParamMinBiasHisto->SetTitle(";Impact Prameter, b (fm);dN/db");
  impactParamMinBiasHisto->GetYaxis()->SetTitleOffset(1.35);
  mScalingFactorHisto->SetTitle(";N_{part};m");
  mScalingFactorHisto->GetYaxis()->SetTitleOffset(1.35);
  
  nPartMinBiasHisto->SetTitleFont(62,"X");
  nPartMinBiasHisto->SetTitleFont(62,"Y");
  nCollMinBiasHisto->SetTitleFont(62,"X");
  nCollMinBiasHisto->SetTitleFont(62,"Y");
  impactParamMinBiasHisto->SetTitleFont(62,"X");
  impactParamMinBiasHisto->SetTitleFont(62,"Y");
  mScalingFactorHisto->SetTitleFont(62,"X");
  mScalingFactorHisto->SetTitleFont(62,"Y");

  nPartMinBiasHisto->SetLabelFont(62,"X");
  nPartMinBiasHisto->SetLabelFont(62,"Y");
  nCollMinBiasHisto->SetLabelFont(62,"X");
  nCollMinBiasHisto->SetLabelFont(62,"Y");
  impactParamMinBiasHisto->SetLabelFont(62,"X");
  impactParamMinBiasHisto->SetLabelFont(62,"Y");
   
  nPartCentralHisto->SetLineColor(kBlack);
  nPartMinBiasHisto->SetLineColor(kBlack);
  nCollCentralHisto->SetLineColor(kBlack);
  nCollMinBiasHisto->SetLineColor(kBlack);
  impactParamCentralHisto->SetLineColor(kBlack);
  impactParamMinBiasHisto->SetLineColor(kBlack);

  nPartCentralHisto->SetFillColor(kRed);
  nCollCentralHisto->SetFillColor(kRed);
  impactParamCentralHisto->SetFillColor(kRed);

  nPartCentralHisto->SetFillStyle(3001);
  nCollCentralHisto->SetFillStyle(3001);
  impactParamCentralHisto->SetFillStyle(3001);
  
  //Loop Over the Entries in the tree
  for (Int_t iEvent=0; iEvent<glauberTree->GetEntries(); iEvent++){
    
    glauberTree->GetEntry(iEvent);

    nPartMinBiasHisto->Fill(glauberEvent->GetNParticipants());
    nCollMinBiasHisto->Fill(glauberEvent->GetNBinaryCollisions());
    impactParamMinBiasHisto->Fill(glauberEvent->GetImpactParameter());

    Double_t x(.12);
    Int_t m = TMath::Ceil(x*glauberEvent->GetNBinaryCollisions() +
			  (1.0 - x)*(0.5 * glauberEvent->GetNParticipants()));

    mScalingFactorHisto->Fill(glauberEvent->GetNParticipants(),m);
    
    if (glauberEvent->GetImpactParameter() <= 2){
      nPartCentralHisto->Fill(glauberEvent->GetNParticipants());
      nCollCentralHisto->Fill(glauberEvent->GetNBinaryCollisions());
      impactParamCentralHisto->Fill(glauberEvent->GetImpactParameter());
    }
    
  }//End Loop over Entries

  canvas1->cd();
  nPartMinBiasHisto->Draw();
  nPartCentralHisto->Draw("SAME");

  canvas2->cd();
  nCollMinBiasHisto->Draw();
  nCollCentralHisto->Draw("SAME");

  canvas3->cd();
  impactParamMinBiasHisto->Draw();
  impactParamCentralHisto->Draw("SAME");

  TLegend *leg = new TLegend(.14,.71,.45,.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(.038);
  leg->SetTextFont(62);
  leg->SetHeader("Glauber MC: ^{197}Au+^{197}Au");
  leg->AddEntry(impactParamMinBiasHisto,"Min-Bias","LF");
  leg->AddEntry(impactParamCentralHisto,"b #leq 2 fm","F");
  leg->Draw("SAME");


  canvas4->cd();
  TF1 *oneToOne = new TF1("oneToOne","x",0,400);
  oneToOne->SetNpx(100);
  oneToOne->SetLineColor(kBlack);
  
 
  mScalingFactorHisto->Draw("COLZ");
  oneToOne->Draw("SAME");  
}
