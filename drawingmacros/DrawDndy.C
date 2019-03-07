
void DrawDndy(){

  //will print to ps and gif directory

  bool print = true;
  
  //Load dndy from the 895, 917 and FXT rootfiles 

  //  TFile * file895     = new TFile("../userfiles/comparisonData/895.root" ,"UPDATE");
  //TFile * file917     = new TFile("../userfiles/comparisonData/917.root" ,"UPDATE");
  //TFile * file917_ref = new TFile("../userfiles/comparisonData/917_reflect.root" ,"UPDATE");

  TFile * fileFxt    = new TFile("../userfiles/AlAu_4_9GeV_2015/analysis_toddcuts/protondNdyFit.root" ,"UPDATE");
  //TFile * fileGaus    = new TFile("../userfiles/AuAu_4_5GeV_2015/analysis/gaus_dNdyFit.root" ,"UPDATE");

  //TGraphErrors * graph895     = (TGraphErrors*) file895->Get("895_dndy");
  //TGraphErrors * graph917     = (TGraphErrors*) file917->Get("917_dndy");
  //TGraphErrors * graph917_ref = (TGraphErrors*) file917_ref->Get("917_reflect_dndy");
  TGraphErrors * graphFxt     = (TGraphErrors*) fileFxt->Get("dNdyGraph_ProtonPlus/dNdyGraphSys_Cent00");
  //TGraphErrors * graphFxt_ref = (TGraphErrors*) fileFxt->Get("dNdyGraph_ProtonPlus/dNdyGraphSysRef_Cent00");
  //TF1 * doubleGaus = (TF1*) fileGaus->Get("doubleGaus");
  //TF1 * PosGaus = (TF1*) fileGaus->Get("positiveGaus");
  //TF1 * NegGaus = (TF1*) fileGaus->Get("negativeGaus");
  
  graphFxt->SetLineWidth(2);
  graphFxt->SetLineColor(kRed);
  graphFxt->SetMarkerColor(kRed);
  graphFxt->SetMarkerStyle(kFullStar);
  graphFxt->SetMarkerSize(2.2);
  
  //graphFxt_ref->SetLineWidth(2);
  //graphFxt_ref->SetLineColor(kRed);
  //graphFxt_ref->SetMarkerColor(kRed);
  //graphFxt_ref->SetMarkerStyle(kOpenStar);
  //graphFxt_ref->SetMarkerSize(2.1);

  //graph895->SetLineWidth(2);
  //graph895->SetMarkerStyle(21);
  //graph895->SetMarkerSize(1.5);

  //graph917_ref->SetLineWidth(2);
  //graph917_ref->SetLineColor(kBlack);
  //graph917_ref->SetMarkerColor(kBlack);
  //graph917_ref->SetMarkerStyle(kOpenCircle);
  //graph917_ref->SetMarkerSize(1.5);
  
  //graph917->SetLineWidth(2);
  //graph917->SetLineColor(kBlack);
  //graph917->SetMarkerColor(kBlack);
  //graph917->SetMarkerStyle(kFullCircle);
  //graph917->SetMarkerSize(1.5);

  TCanvas * canvas = new TCanvas("dNdy","dNdy",800,850);
  canvas->SetLeftMargin(.15);
  canvas->SetBottomMargin(.1);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameLineWidth(2);

  TH1F * frame;

  frame = canvas->DrawFrame(-1.5,0,0,110);
  frame->GetXaxis()->SetLabelSize(.040);
  frame->GetXaxis()->SetTitle("y - y_{CM}");
  frame->GetYaxis()->SetLabelSize(.040);
  frame->GetYaxis()->SetTitle("dN/dy");
  frame->GetYaxis()->SetLabelOffset(0.01);
  frame->GetYaxis()->SetTitleOffset(1.7);

  //graph895->Draw("P");
  //graph917->Draw("P");
  //graph917_ref->Draw("P");
  graphFxt->Draw("P");
  //graphFxt_ref->Draw("P");
  
  //doubleGaus->Draw("SAME");
  //NegGaus->Draw("SAME");
  //PosGaus->Draw("SAME");

  TPaveText * pt3 = new TPaveText(0.43,0.83,0.63,0.86,"BRNBNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(0);
  pt3->SetTextFont(63);
  pt3->SetTextSize(23);
  TText * text3 = pt3->AddText("Open symbols are reflected");
  pt3->Draw();

  //Build legend and add 895, 917, FXT entries

  TLegend *leg = new TLegend(0.36,0.70,0.56,0.82,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(.030);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TLegendEntry *entry1=leg->AddEntry("TMarker","STAR FXT #sqrt{s_{NN}} = 4.9 GeV","P");
  entry1->SetLineColor(1);
  entry1->SetLineStyle(1);
  entry1->SetLineWidth(1);
  entry1->SetMarkerColor(kRed);
  entry1->SetMarkerStyle(29);
  entry1->SetMarkerSize(1.3);
  entry1->SetTextFont(42);
  //  TLegendEntry * entry2=leg->AddEntry("TMarker","E895 #sqrt{s_{NN}} = 4.3 GeV","P");
  //entry2->SetLineColor(1);
  //entry2->SetLineStyle(1);
  //entry2->SetLineWidth(1);
  //entry2->SetMarkerColor(1);
  //entry2->SetMarkerStyle(21);
  //entry2->SetMarkerSize(1.3);
  //entry2->SetTextFont(42);
  //TLegendEntry * entry3=leg->AddEntry("TMarker","E917 #sqrt{s_{NN}} = 4.3 GeV","P");
  //entry3->SetLineColor(1);
  //entry3->SetLineStyle(1);
  //entry3->SetLineWidth(1);
  //entry3->SetMarkerColor(1);
  //entry3->SetMarkerStyle(20);
  //entry3->SetMarkerSize(1.3);
  //entry3->SetTextFont(42);
  leg->Draw();

  //Print

  if (print==true){
    canvas->Print("../userfiles/AlAu_4_9GeV_2015/images/Proton_AlAu_FXT_dndy.gif");
    canvas->Print("../userfiles/AlAu_4_9GeV_2015/images/Proton_AlAu_FXT_dndy.ps");
    canvas->Print("../userfiles/AlAu_4_9GeV_2015/images/Proton_AlAu_FXT_dndy.eps");
    canvas->Print("../userfiles/AlAu_4_9GeV_2015/images/Proton_AlAu_FXT_dndy.pdf");
  }


}
