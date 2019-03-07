//Draw the Woods Saxon Density Profiles of different
//species of nuclei


void DrawDensityProfiles(){

  gSystem->Load("../bin/GlauberClass_cxx.so");


  GlauberClass *g1 = new GlauberClass(197,27,0);
  GlauberClass *g2 = new GlauberClass(208,238,0);

  TF1 *gold  = g1->GetWoodsSaxon1DGivenTheta(0,0);
  TF1 *alum  = g1->GetWoodsSaxon1DGivenTheta(0,1);
  TF1 *lead  = g2->GetWoodsSaxon1DGivenTheta(0,0);
  TF1 *uran1 = g2->GetWoodsSaxon1DGivenTheta(0,1);
  TF1 *uran2 = g2->GetWoodsSaxon1DGivenTheta(TMath::Pi()/2.0,1);
  
  gold->SetLineWidth(3);
  alum->SetLineWidth(3);
  lead->SetLineWidth(3);
  uran1->SetLineWidth(3);
  uran2->SetLineWidth(3);
  
  TCanvas *canvas = new TCanvas("canvas","canvas",20,20,800,600);
  canvas->SetTopMargin(.05);
  canvas->SetRightMargin(.05);
  canvas->SetLeftMargin(.12);
  canvas->SetTicks();

  TH1F *frame = canvas->DrawFrame(0.0,0.0,10.5,.195);
  frame->SetTitle(";r (fm); #rho(r,#theta) (fm^{-3})");
  frame->GetXaxis()->SetTitleSize(.045);
  frame->GetYaxis()->SetTitleSize(.045);
  frame->SetTitleFont(42,"X");
  frame->SetTitleFont(42,"Y");
  frame->GetYaxis()->SetTitleOffset(1.3);
  
  alum->SetLineColor(kGray);
  lead->SetLineColor(kBlack);
  gold->SetLineColor(kGray+1);
  uran1->SetLineColor(kGray+2);
  uran2->SetLineColor(kGray+3);

  alum->SetLineStyle(2);
  lead->SetLineStyle(5);
  gold->SetLineStyle(1);
  uran1->SetLineStyle(3);
  uran2->SetLineStyle(7);

  alum->Draw("SAME");
  lead->Draw("SAME");
  gold->Draw("SAME");
  uran1->Draw("SAME");
  uran2->Draw("SAME");

  TPaveText *title = new TPaveText(.47,.84,.92,.92,"brNDC");
  title->SetFillColor(kWhite);
  title->SetBorderSize(0);
  title->SetTextSize(.045);
  title->SetTextFont(62);
  title->AddText("Nuclear Density Distributions");
  title->Draw("SAME");
  
  TLegend *leg = new TLegend(.72,.57,.92,.84);
  leg->SetTextSize(.038);
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->AddEntry(alum,"^{27}Al, #theta=0","L");
  leg->AddEntry(gold,"^{197}Au, #theta=0","L");
  leg->AddEntry(lead,"^{208}Pb, #theta=0","L");
  leg->AddEntry(uran1,"^{238}U, #theta=0","L");
  leg->AddEntry(uran2,"^{238}U, #theta=#pi/2","L");
  leg->Draw("SAME");

  
}
