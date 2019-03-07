//_______________________________________________________________________
TPaveText *MakeTitle(TGraphErrors *graph, Double_t energy,
		     Double_t x1=.46, Double_t y1=.71, Double_t x2=.91, Double_t y2=.91, double size=.06){

  double textSize(size);
  TPaveText *title = new TPaveText(x1,y1,x2,y2,"NDC");
  title->SetTextFont(62);
  title->SetTextSize(textSize);
  title->SetTextAlign(11);
  title->SetFillColor(kWhite);
  title->SetBorderSize(0);
  TString graphTitle(graph->GetTitle());
  TString graphSubTitle(graphTitle(0,graphTitle.First(";")));
  TObjArray *tokens = graphSubTitle.Tokenize("|");

  if (tokens->GetEntries() == 0){
    cout <<"WARNING - DrawingFunctions::MakeTitle() - No Tokens Found. You may need to add a ;\n";
    return title;
  }

  int line(0);
  while (line < tokens->GetEntries()){
    if (line == 1){
      title->AddText(Form("#sqrt{s_{NN}} = %.3g GeV |%s",
			  energy,
			  ((TObjString*)tokens->At(line))->GetString().Data()));
    }
    else {
      title->AddText(((TObjString*)tokens->At(line))->GetString().Data());
      
    }

    if (line != 0)
      title->GetLine(line)->SetTextSize(0.8 * textSize);

    
    line++;
  }


  delete tokens;
  return title;
  
  title->AddText(((TObjString*)tokens->At(0))->GetString().Data());
  title->AddText(Form("#sqrt{s_{NN}} = %.3g GeV |%s",
		      energy,
		      ((TObjString*)tokens->At(1))->GetString().Data()));
  title->GetLine(1)->SetTextSize(.045);
  
  TString nextLine = TString::Format("%s",TString(((TObjString*)tokens->At(2))->GetString())(1,100).Data());


  if (tokens->GetEntries() > 3){
    nextLine += TString::Format("| %s",TString(((TObjString*)tokens->At(3))->GetString())(1,100).Data());
  }
  title->AddText(Form("%s",nextLine.Data()));
  title->GetLine(2)->SetTextSize(.045);

  if (tokens->GetEntries() > 4){
    title->AddText(Form("%s",TString(((TObjString*)tokens->At(4))->GetString())(1,100).Data()));
    title->GetLine(3)->SetTextSize(.045);
  }

  delete tokens;

  return title;

}

//_______________________________________________________________________
TH1F *MakeFrame(TPad *pad, Double_t yMax=.025, Double_t xMax = 2.0){

  pad->SetTopMargin(.02);
  pad->SetRightMargin(.02);
  pad->SetBottomMargin(.1);
  pad->SetLeftMargin(.15);
  pad->SetTicks();

  double textSize(.045);
  
  TH1F *frame = pad->DrawFrame(0,0,xMax,yMax);  
  frame->SetTitleSize(textSize,"X");
  frame->SetTitleFont(62,"X");
  frame->GetXaxis()->SetLabelSize(textSize);
  
  frame->SetTitleSize(textSize,"Y");
  frame->SetTitleFont(62,"Y");
  frame->GetYaxis()->SetLabelSize(textSize);
  
  frame->SetTitleOffset(1.15,"Y");

  return frame;

}

//_______________________________________________________________________
TLegend *MakeLegend(float x0=.58, float y0=.52, float x1=.88, float y1=.72){

  TLegend *leg = new TLegend(x0,y0,x1,y1);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(.045);

  return leg;
}
