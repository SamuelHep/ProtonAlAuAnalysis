void GetGlauberStats(TFile *file, Double_t &nPartCentralMean, Double_t &nPartCentralRMS,
		     Double_t &nCollCentralMean, Double_t &nCollCentralRMS,
		     Double_t &nPartMinBiasMean, Double_t &nPartMinBiasRMS,
		     Double_t &nCollMinBiasMean, Double_t &nCollMinBiasRMS,
		     Double_t &nCentral, Double_t &nMinBias);

TGraphErrors *ratio(TGraphErrors *gNum, TGraphErrors *gDenom);

void PlotGlauberResultsVariousSpecies(){

  gSystem->Load("../bin/GlauberClass_cxx");

  TCanvas *canvas = new TCanvas("canvas","canvas",20,20,700,1000);
  
  //Upper Pad
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->SetTopMargin(.05);
  pad1->SetRightMargin(.15);
  pad1->SetLeftMargin(.15);
  pad1->SetTickx();
  pad1->Draw();
  pad1->cd();
  
  Double_t leftMax = 575;
  TH1F *frame = gPad->DrawFrame(0,0,275,leftMax);
  frame->SetTitle(";Mass Number (A); <N_{part}>");
  frame->GetYaxis()->SetTitleSize(.045);
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetTitleFont(62,"X");
  frame->SetTitleFont(62,"Y");
  frame->SetLabelFont(62,"X");
  frame->SetLabelFont(62,"Y");
  frame->GetXaxis()->SetLabelOffset(999);
  frame->GetXaxis()->SetLabelSize(0);

  //draw an axis on the right side
  Double_t rightMax = 1650;
  Double_t scale = leftMax/rightMax;
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),0,rightMax,510,"+L");
  axis->SetTitle("<N_{coll}>");
  axis->SetTitleSize(.045);
  axis->SetTitleColor(kRed);
  axis->SetTitleOffset(1.6);
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->Draw();

  //Lower Pad
  canvas->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,.00,1,0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.28);
  pad2->SetRightMargin(.15);
  pad2->SetLeftMargin(.15);
  pad2->SetTicks();
  pad2->Draw();
  pad2->cd();

  TH1F *frame1 = gPad->DrawFrame(0,0.75,275,1.15);
  frame1->SetTitle(";Mass Number (A); FD/Sphere");
  frame1->GetYaxis()->SetTitleSize(25);
  frame1->GetXaxis()->SetTitleSize(25);
  frame1->GetYaxis()->SetTitleOffset(2.4);
  frame1->GetXaxis()->SetTitleOffset(3.8);
  frame1->GetYaxis()->CenterTitle();
  frame1->SetTitleFont(63,"X");
  frame1->SetTitleFont(63,"Y");
  frame1->SetLabelFont(63,"X");
  frame1->SetLabelFont(63,"Y");
  frame1->GetYaxis()->SetLabelSize(25);
  frame1->GetXaxis()->SetLabelSize(25);
  frame1->GetYaxis()->SetNdivisions(5);
  
  //return;

  //Graphs for the Woods Saxon Distributions
  TGraphErrors *nPartGraphCentral = new TGraphErrors();
  TGraphErrors *nPartGraphMinBias = new TGraphErrors();
  nPartGraphCentral->SetMarkerStyle(kFullCircle);
  nPartGraphMinBias->SetMarkerStyle(kOpenCircle);
  nPartGraphCentral->SetMarkerColor(kBlack);
  nPartGraphMinBias->SetMarkerColor(kBlack);
  nPartGraphCentral->SetMarkerSize(1.2);
  nPartGraphMinBias->SetMarkerSize(1.2);
  
  TGraphErrors *nCollGraphCentral = new TGraphErrors();
  TGraphErrors *nCollGraphMinBias = new TGraphErrors();
  nCollGraphCentral->SetMarkerStyle(kFullSquare);
  nCollGraphMinBias->SetMarkerStyle(kOpenSquare);
  nCollGraphCentral->SetMarkerColor(kRed);
  nCollGraphMinBias->SetMarkerColor(kRed);
  nCollGraphCentral->SetMarkerSize(1.2);
  nCollGraphMinBias->SetMarkerSize(1.2);

  //Graphs for the Uniform Sphere Distributions
  TGraphErrors *nPartGraphCentral_Sphere = new TGraphErrors();
  TGraphErrors *nPartGraphMinBias_Sphere = new TGraphErrors();
  nPartGraphCentral_Sphere->SetMarkerStyle(kFullDiamond);
  nPartGraphMinBias_Sphere->SetMarkerStyle(kOpenDiamond);
  nPartGraphCentral_Sphere->SetMarkerColor(kBlack);
  nPartGraphMinBias_Sphere->SetMarkerColor(kBlack);
  nPartGraphCentral_Sphere->SetMarkerSize(1.2);
  nPartGraphMinBias_Sphere->SetMarkerSize(1.2);
  
  TGraphErrors *nCollGraphCentral_Sphere = new TGraphErrors();
  TGraphErrors *nCollGraphMinBias_Sphere = new TGraphErrors();
  nCollGraphCentral_Sphere->SetMarkerStyle(kFullCross);
  nCollGraphMinBias_Sphere->SetMarkerStyle(kOpenCross);
  nCollGraphCentral_Sphere->SetMarkerColor(kRed);
  nCollGraphMinBias_Sphere->SetMarkerColor(kRed);
  nCollGraphCentral_Sphere->SetMarkerSize(1.2);
  nCollGraphMinBias_Sphere->SetMarkerSize(1.2);

  canvas->cd();
  pad1->cd();
  
  TLegend *leg = new TLegend(.18,.71,.58,.92);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(.03);
  leg->SetTextFont(62);
  leg->SetHeader("Fermi-Distributions");
  leg->AddEntry(nCollGraphCentral,"<N_{coll}> (b #leq 2 fm)","P");
  leg->AddEntry(nPartGraphCentral,"<N_{part}> (b #leq 2 fm)","P");
  leg->AddEntry(nCollGraphMinBias,"<N_{coll}> (Min Bias)","P");
  leg->AddEntry(nPartGraphMinBias,"<N_{part}> (Min Bias)","P");
  leg->Draw("SAME");

  TLegend *leg1 = new TLegend(.18,.5,.58,.71);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(kWhite);
  leg1->SetTextSize(.03);
  leg1->SetTextFont(62);
  leg1->SetHeader("Uniform Hard-Sphere");
  leg1->AddEntry(nCollGraphCentral_Sphere,"<N_{coll}> (b #leq 2 fm)","P");
  leg1->AddEntry(nPartGraphCentral_Sphere,"<N_{part}> (b #leq 2 fm)","P");
  leg1->AddEntry(nCollGraphMinBias_Sphere,"<N_{coll}> (Min Bias)","P");
  leg1->AddEntry(nPartGraphMinBias_Sphere,"<N_{part}> (Min Bias)","P");
  leg1->Draw("SAME");
  
  const int nFiles = 14;
  Int_t nNucleons[nFiles] =
    {16, 27, 35, 40, 56, 63, 70, 93, 152, 184, 197, 208, 238, 244};

  //WOODS SAXON GRAPHS
  //Loop Over the Files 
  for (Int_t iFile=0; iFile<nFiles; iFile++){

    TFile file(Form("../data/Glauber_%d_%d_42mb_WoodsSaxon.root",nNucleons[iFile],
		    nNucleons[iFile]));

    Double_t nPartCentralMean(0), nPartCentralRMS(0);
    Double_t nCollCentralMean(0), nCollCentralRMS(0);
    Double_t nPartMinBiasMean(0), nPartMinBiasRMS(0);
    Double_t nCollMinBiasMean(0), nCollMinBiasRMS(0);
    Double_t nCentral(0), nMinBias(0);

    GetGlauberStats(&file,nPartCentralMean,nPartCentralRMS,
		    nCollCentralMean, nCollCentralRMS,
		    nPartMinBiasMean, nPartMinBiasRMS,
		    nCollMinBiasMean, nCollMinBiasRMS,
		    nCentral, nMinBias);

    file.Close();

    cout <<Form("%d\t%g\t%g\t%g\t%g\n",
		nNucleons[iFile],nPartCentralMean,nPartMinBiasMean,
		nCollCentralMean,nCollMinBiasMean);
    
    //nPart Graphs
    nPartGraphCentral->SetPoint(nPartGraphCentral->GetN(),nNucleons[iFile],
				nPartCentralMean);
    nPartGraphCentral->SetPointError(nPartGraphCentral->GetN()-1,0,
				     nPartCentralRMS/sqrt(nCentral));

    nPartGraphMinBias->SetPoint(nPartGraphMinBias->GetN(),nNucleons[iFile],
				nPartMinBiasMean);
    nPartGraphMinBias->SetPointError(nPartGraphMinBias->GetN()-1,0,
				     nPartMinBiasRMS/sqrt(nMinBias));

    //nColl Graphs
    nCollGraphCentral->SetPoint(nCollGraphCentral->GetN(),nNucleons[iFile],
				nCollCentralMean*scale);
    nCollGraphCentral->SetPointError(nCollGraphCentral->GetN()-1,0,
				     (nCollCentralRMS/sqrt(nCentral))*scale);

    nCollGraphMinBias->SetPoint(nCollGraphMinBias->GetN(),nNucleons[iFile],
				nCollMinBiasMean*scale);
    nCollGraphMinBias->SetPointError(nCollGraphMinBias->GetN()-1,0,
				     (nCollMinBiasRMS/sqrt(nMinBias))*scale);

    canvas->cd();
    pad1->cd();
    nPartGraphCentral->Draw("P");
    nPartGraphMinBias->Draw("P");
    nCollGraphCentral->Draw("P");
    nCollGraphMinBias->Draw("P");
    canvas->Update();    
    
  }//End Loop Over Files


  //UNIFORM SPHERE GRAPHS
  //Loop Over the Files 
  for (Int_t iFile=0; iFile<nFiles; iFile++){

    TFile file(Form("../data/Glauber_%d_%d_42mb_UniformSphere.root",nNucleons[iFile],
		    nNucleons[iFile]));

    Double_t nPartCentralMean(0), nPartCentralRMS(0);
    Double_t nCollCentralMean(0), nCollCentralRMS(0);
    Double_t nPartMinBiasMean(0), nPartMinBiasRMS(0);
    Double_t nCollMinBiasMean(0), nCollMinBiasRMS(0);
    Double_t nCentral(0), nMinBias(0);

    GetGlauberStats(&file,nPartCentralMean,nPartCentralRMS,
		    nCollCentralMean, nCollCentralRMS,
		    nPartMinBiasMean, nPartMinBiasRMS,
		    nCollMinBiasMean, nCollMinBiasRMS,
		    nCentral, nMinBias);

    file.Close();

    cout <<Form("%d\t%g\t%g\t%g\t%g\n",
		nNucleons[iFile],nPartCentralMean,nPartMinBiasMean,
		nCollCentralMean,nCollMinBiasMean);
    
    //nPart Graphs
    nPartGraphCentral_Sphere->SetPoint(nPartGraphCentral_Sphere->GetN(),nNucleons[iFile],
				nPartCentralMean);
    nPartGraphCentral_Sphere->SetPointError(nPartGraphCentral_Sphere->GetN()-1,0,
				     nPartCentralRMS/sqrt(nCentral));

    nPartGraphMinBias_Sphere->SetPoint(nPartGraphMinBias_Sphere->GetN(),nNucleons[iFile],
				nPartMinBiasMean);
    nPartGraphMinBias_Sphere->SetPointError(nPartGraphMinBias_Sphere->GetN()-1,0,
				     nPartMinBiasRMS/sqrt(nMinBias));

    //nColl Graphs
    nCollGraphCentral_Sphere->SetPoint(nCollGraphCentral_Sphere->GetN(),nNucleons[iFile],
				nCollCentralMean*scale);
    nCollGraphCentral_Sphere->SetPointError(nCollGraphCentral_Sphere->GetN()-1,0,
				     (nCollCentralRMS/sqrt(nCentral))*scale);

    nCollGraphMinBias_Sphere->SetPoint(nCollGraphMinBias_Sphere->GetN(),nNucleons[iFile],
				nCollMinBiasMean*scale);
    nCollGraphMinBias_Sphere->SetPointError(nCollGraphMinBias_Sphere->GetN()-1,0,
				     (nCollMinBiasRMS/sqrt(nMinBias))*scale);

    canvas->cd();
    pad1->cd();    
    nPartGraphCentral_Sphere->Draw("P");
    nPartGraphMinBias_Sphere->Draw("P");
    nCollGraphCentral_Sphere->Draw("P");
    nCollGraphMinBias_Sphere->Draw("P");
    canvas->Update();    
    
  }//End Loop Over Files

  //Plot the Ratios

  TGraphErrors *nPartCentralRatio = ratio(nPartGraphCentral,nPartGraphCentral_Sphere);
  TGraphErrors *nCollCentralRatio = ratio(nCollGraphCentral,nCollGraphCentral_Sphere);
  TGraphErrors *nPartMinBiasRatio = ratio(nPartGraphMinBias,nPartGraphMinBias_Sphere);
  TGraphErrors *nCollMinBiasRatio = ratio(nCollGraphMinBias,nCollGraphMinBias_Sphere);

  nPartCentralRatio->SetMarkerStyle(kFullCircle);
  nPartMinBiasRatio->SetMarkerStyle(kOpenCircle);
  nCollCentralRatio->SetMarkerStyle(kFullSquare);
  nCollMinBiasRatio->SetMarkerStyle(kOpenSquare);

  nPartCentralRatio->SetMarkerColor(kBlack);
  nPartMinBiasRatio->SetMarkerColor(kBlack);
  nCollCentralRatio->SetMarkerColor(kRed);
  nCollMinBiasRatio->SetMarkerColor(kRed);
 

  canvas->cd();
  pad2->cd();

  TLine *unity = new TLine(0,1,275,1);
  unity->SetLineWidth(2);
  unity->SetLineColor(kGray+2);
  unity->SetLineStyle(7);
  unity->Draw("SAME");
  
  nPartCentralRatio->Draw("P");
  nCollCentralRatio->Draw("P");
  nPartMinBiasRatio->Draw("P");
  nCollMinBiasRatio->Draw("P");
  
}



//_______________________________________________________________________________________
void GetGlauberStats(TFile *file, Double_t &nPartCentralMean, Double_t &nPartCentralRMS,
		     Double_t &nCollCentralMean, Double_t &nCollCentralRMS,
		     Double_t &nPartMinBiasMean, Double_t &nPartMinBiasRMS,
		     Double_t &nCollMinBiasMean, Double_t &nCollMinBiasRMS,
		     Double_t &nCentral, Double_t &nMinBias){


  TTree *glauberTree = (TTree *)file->Get("GlauberTree");
  GlauberClass *glauberEvent = 0;
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  
  std::vector<int> nPartCentral;
  std::vector<int> nPartMinBias;
  std::vector<int> nCollCentral;
  std::vector<int> nCollMinBias;
  
  //Loop Over the Entries in the tree
  for (Int_t iEvent=0; iEvent<glauberTree->GetEntries(); iEvent++){
    
    glauberTree->GetEntry(iEvent);
    
    nPartMinBias.push_back(glauberEvent->GetNParticipants());
    nCollMinBias.push_back(glauberEvent->GetNBinaryCollisions());
    
    if (glauberEvent->GetImpactParameter() <= 2){
      nPartCentral.push_back(glauberEvent->GetNParticipants());
      nCollCentral.push_back(glauberEvent->GetNBinaryCollisions());
    }
    
  }//End Loop over Entries
  
  
  nPartCentralMean = TMath::Mean(nPartCentral.size(),
				 &nPartCentral.at(0));
  nPartCentralRMS   = TMath::RMS(nPartCentral.size(),
				 &nPartCentral.at(0));
  
  nPartMinBiasMean = TMath::Mean(nPartMinBias.size(),
				 &nPartMinBias.at(0));
  nPartMinBiasRMS  = TMath::RMS(nPartMinBias.size(),
				&nPartMinBias.at(0));
  
  nCollCentralMean = TMath::Mean(nCollCentral.size(),
				 &nCollCentral.at(0));
  nCollCentralRMS   = TMath::RMS(nCollCentral.size(),
				 &nCollCentral.at(0));
  
  nCollMinBiasMean = TMath::Mean(nCollMinBias.size(),
				 &nCollMinBias.at(0));
  nCollMinBiasRMS  = TMath::RMS(nCollMinBias.size(),
				&nCollMinBias.at(0));

  nCentral = (double)nPartCentral.size();
  nMinBias = (double)nPartMinBias.size();

  delete glauberTree;
  delete glauberEvent;
  
  return;
    
}

//________________________________________________________________
TGraphErrors *ratio(TGraphErrors *gNum, TGraphErrors *gDenom){

  TGraphErrors *g = new TGraphErrors();
  
  for (Int_t i=0; i<gNum->GetN(); i++){

    Double_t x(0), ratio(0), ratioErr(0);

    x = gNum->GetX()[i];
    ratio = gNum->GetY()[i] / gDenom->GetY()[i];
    ratioErr = fabs(ratio) * sqrt(pow(gNum->GetEY()[i]/gNum->GetY()[i],2) +
				  pow(gDenom->GetEY()[i]/gDenom->GetY()[i],2));

    g->SetPoint(i,x,ratio);
    g->SetPointError(i,0,ratioErr);
    
  }

  return g;

}
