Double_t rapidityMin = -2.05;
Double_t rapidityMax = 2.05;
Double_t rapidityBinWidth =  0.1;

Double_t mTm0Min = 0.0;
Double_t mTm0Max = 2.0;
Double_t mTm0BinWidth = 0.025;

const int nRapidityBins = 41;
const int nmTm0Bins = 80;
const int nCentBins = 6;
const int funcBins = 700;


//__________________________________________________________________________

Double_t ThermalFitFunc(Double_t *x, Double_t *par){

  //NOTE: par2 corresponding to the particle mass should be fixed
  
  Double_t xx =x[0];

  Double_t dNdy = par[0];
  Double_t tSlope = par[1];
  Double_t pMass = par[2];  //This Should be Fixed

  return (1.0/(2.0*TMath::Pi())) * (dNdy/tSlope) *(xx+pMass) *(1.0/(2.0*tSlope*tSlope + 2.0*tSlope*pMass + pMass*pMass)) * exp(-(xx)/tSlope);

}

//__________________________________________________________________________


TGraphErrors *TGraphScale(TGraphErrors *graph, Double_t scaleFactor, Bool_t returnNew=true){

  //Returns a new TGraph scaled by scaleFactor
  //By default returnNew==true and this function returns
  
  TGraphErrors *scaledGraph = NULL;
  if (returnNew){
    scaledGraph = new TGraphErrors();
    scaledGraph->SetName(graph->GetName());
    scaledGraph->SetTitle(graph->GetTitle());
    scaledGraph->SetMarkerStyle(graph->GetMarkerStyle());
    scaledGraph->SetMarkerColor(graph->GetMarkerColor());
  }
  else
    scaledGraph = graph;

  Int_t nPoints = graph->GetN();
  for (Int_t iPoint=0; iPoint<nPoints; iPoint++){

    scaledGraph->SetPoint(iPoint,
                          graph->GetX()[iPoint],
                          graph->GetY()[iPoint]*scaleFactor);
    scaledGraph->SetPointError(iPoint,
                               graph->GetEX()[iPoint],
                               graph->GetEY()[iPoint]*scaleFactor);

  }

  return scaledGraph;
}

//___________________________________________________________________________
void RemovePointsWithLargeErrors(TGraphErrors *spectrum, Double_t maxRelativeError=.1){

  //Loop Over the Points of the spectrum. Remove any point which is found
  //to have a relativeError larger than maxRelativeError
  for (int iPoint=spectrum->GetN()-1; iPoint>=0; iPoint--){

    if (spectrum->GetEY()[iPoint] / spectrum->GetY()[iPoint] > maxRelativeError)
      spectrum->RemovePoint(iPoint);
  }

}


void TGraphChop(TGraphErrors *graph, Double_t threshold, Bool_t below){

  //Remove points from graph that are below (or above) the
  //the x value of threshold

  for (Int_t iPoint=graph->GetN()-1; iPoint>=0; iPoint--){

    //If the user wants to remove points above the threshold
    if (!below && graph->GetX()[iPoint] > threshold){
      graph->RemovePoint(iPoint);
    }

    else if (below && graph->GetX()[iPoint] < threshold) {
      graph->RemovePoint(iPoint);
    }

  }

}

//___________________________________________________________                                                                                                                                               
Double_t GetRapidityRangeLow(Int_t rapidityIndex){

  return rapidityMin + (rapidityIndex * rapidityBinWidth);
}

//___________________________________________________________                                                                                                                                               
Double_t GetRapidityRangeHigh(Int_t rapidityIndex){

  return GetRapidityRangeLow(rapidityIndex)+rapidityBinWidth;
}

//___________________________________________________________                                                                                                                                               
Double_t GetRapidityRangeCenter(Int_t rapidityIndex){

  return (GetRapidityRangeLow(rapidityIndex) + GetRapidityRangeHigh(rapidityIndex)) / 2.0;
}



void DrawAllCent(TString spectraFile, TString eventConfig, TString system, Double_t energy=4.5,Int_t centBin = 0,
		  Double_t startmTm0 = 0.0,
		  Bool_t corrected=true){

  cout << "DEBUG::0" << endl;

  Bool_t save = true;        //Should the canvas be saved?
  Double_t scaleFactor = 5.0; //The factor to scale the non- midrapidity spectra by
  Int_t yMin = 7;            //Minimum rapidity bin index to draw //Collider 11, PosY 18, NegY
  Int_t yMax = 18;            //Maximum rapidity bin index to draw //Collider 29, PosY 33, NegY
  Double_t maxmTm0=1.0;       //Maximum mT-m0 value to draw
  Double_t midY =-1.5;

  Double_t yScaleMin(0), yScaleMax(0);
  
  yScaleMin = .05;
  yScaleMax = 300000000000;

  //Load the root file
  TFile *file = new TFile(spectraFile,"READ");

  //Create the Spectra Name
  TString type = "corrected";
  TString Type = "Corrected";
  TString name = "Corrected";
 
  TString speciesName = "ProtonPlus";
  Int_t midRapidityIndex = 5;

  //Create Array of TGraphErrors
  TGraphErrors *spectra;
  TGraphErrors *scaledSpectra;
  TF1 *spectraFit;

  TPaveText *label;
  TPaveText *label2;

  //Loop Over the centrality Bins
  TCanvas * canvas = new TCanvas("SpectraPlots","SpectraPlots",1200,900);

  canvas->Divide(3,2,0.,0.);

  for (Int_t centIndex=0; centIndex<=5; centIndex++){

    canvas->cd(centIndex+1);
    
    //Create the Canvas and Frame
    //    TPad *pad = (TPad*) gROOT->FindObject(name);
    //    if (pad) delete pad;

    //        TPad = new TPad(Form("%s%02d_%s_%s_Cent%02d",
    //                            system.Data(),(int)energy,speciesName.Data(),
    //                            eventConfig.Data(),centBin),"title",0,400,300,800);
    

      
    //    gPad->SetTopMargin(.05);
    //    gPad->SetBottomMargin(.17);
    //    gPad->SetRightMargin(.17);
    //    gPad->SetLeftMargin(.17);
    //    gPad->SetLogy();
    //    gPad->SetTicks(1,1);
          
    switch (centIndex){
      
    case 0:
      gPad->SetTopMargin(.05);
      gPad->SetBottomMargin(.0);
      gPad->SetRightMargin(.0);
      gPad->SetLeftMargin(.17);
      gPad->SetLogy();
      gPad->SetTicks(1,1);
      break;
    case 1:
      gPad->SetTopMargin(.05);
      gPad->SetBottomMargin(.0);
      gPad->SetRightMargin(0);
      gPad->SetLeftMargin(0);
      gPad->SetLogy();
      gPad->SetTicks(1,1);
      break;
    case 2:
      gPad->SetTopMargin(.05);
      gPad->SetBottomMargin(0);
      gPad->SetRightMargin(0.02);
      gPad->SetLeftMargin(0);
      gPad->SetLogy();
      gPad->SetTicks(1,1);
      break;
    case 3:
      gPad->SetTopMargin(0);
      gPad->SetBottomMargin(.17);
      gPad->SetRightMargin(.0);
      gPad->SetLeftMargin(.17);
      gPad->SetLogy();
      gPad->SetTicks(1,1);
      break;
    case 4:
      gPad->SetTopMargin(0);
      gPad->SetBottomMargin(.17);
      gPad->SetRightMargin(0);
      gPad->SetLeftMargin(0);
      gPad->SetLogy();
      gPad->SetTicks(1,1);
      break;
    case 5:
      gPad->SetTopMargin(0);
      gPad->SetBottomMargin(.17);
      gPad->SetRightMargin(0.02);
      gPad->SetLeftMargin(0);
      gPad->SetLogy();
      gPad->SetTicks(1,1);
      break;
    }
      
    //        TH1F *frame = pad->DrawFrame(0,yScaleMin,maxmTm0+.4,yScaleMax);
    TH1F *frame = gPad->DrawFrame(0.01,yScaleMin,maxmTm0+.3,yScaleMax);
    frame->GetXaxis()->SetTitle("m_{T}-m_{p} (GeV/c^{2})");
    frame->GetYaxis()->SetTitle("#frac{1}{N_{Evt}}#times#frac{1}{2#pim_{T}}#times#frac{d^{2}N}{dm_{T}dy} (GeV/c^{2})^{-2}");
    frame->GetYaxis()->SetTitleOffset(2.9);
    frame->GetYaxis()->SetTitleFont(63);
    frame->GetYaxis()->SetTitleSize(21);
    frame->GetXaxis()->SetTitleFont(63);
    frame->GetXaxis()->SetTitleSize(21);
    frame->GetXaxis()->SetTitleOffset(2);

    TPaveText *title = new TPaveText(.56,.82,.91,.92,"NBNDCBR");
    title->SetFillColor(kWhite);
    title->SetBorderSize(0);
    title->SetTextSize(.05);
    title->SetTextAlign(32);
    title->AddText(Form("p Spectra %s",
                        system.Data()));
    title->AddText(Form("#sqrt{s_{NN}} = %.03g GeV",energy));
    
    title->GetLine(1)->SetTextSize(.05);

    if (centIndex == 0)    label2 = new TPaveText(.25,.04,.56,.22,"NBNDCBR");
    if (centIndex == 3)    label2 = new TPaveText(.25,.21,.56,.38,"NBNDCBR");
    if (centIndex == 1 || centIndex == 2)    label2 = new TPaveText(.21,.04,.4,.22,"NBNDCBR");
    if (centIndex == 4 || centIndex == 5)    label2 = new TPaveText(.21,.20,.4,.38,"NBNDCBR");
    label2->SetFillColor(kWhite);
    label2->SetBorderSize(0);
    label2->SetTextSize(.050);
    if (centIndex == 3)  label2->SetTextSize(.045);
    label2->SetTextAlign(32);
    label2->AddText(Form("%d-%d%% Central",
			 centIndex*5,
			 (centIndex+1)*5 )) ;

    label2->Draw("SAME");

    if (centIndex ==1  )title->Draw("SAME");

    TMarker *midYMarker = new TMarker(0,0,kFullCircle);
    midYMarker->SetMarkerColor(kRed);
    midYMarker->SetMarkerSize(1.2);

    TMarker *notMidYMarker = new TMarker(0,0,kFullCircle);
    notMidYMarker->SetMarkerColor(kBlack);
    notMidYMarker->SetMarkerSize(1.4);

    TLegend *leg = new TLegend(.65,.76,.92,.82);
    leg->SetLineColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(.035);
    

    leg->AddEntry(notMidYMarker,Form("x %g^{#pm n}",scaleFactor),"P");
    if (centIndex ==1) leg->Draw("SAME");





    TPaveText *centralityTitle = new TPaveText(.2,.13,.47,.19,"BRNBNDC");
    //        centralityTitle->SetFillColor(kWhite);
    //        centralityTitle->SetBorderSize(0);
    //        centralityTitle->SetTextFont(30);
    //        centralityTitle->SetTextSize(25);
    //        centralityTitle->AddText(Form("%d-%d%% Central",
    //                                centIndex*5,
    //                          (centIndex+1)*5));
    
    //if (centIndex == 0) centralityTitle->Draw("SAME");

    /*    TPaveText * starPrelim = new TPaveText(.55,.13,.90,.19,"BRNBNDC");
    starPrelim->SetFillColor(kWhite);
    starPrelim->SetBorderSize(0);
    starPrelim->SetTextFont(30);
    starPrelim->SetTextSize(18);
    starPrelim->AddText("STAR PRELIMINARY");
    if (centIndex ==0 ) starPrelim->Draw("SAME");

    */




    gPad->Update();

  }


  //Loop Over the rapidity bins
  for (Int_t centIndex=0; centIndex<=5; centIndex++){
    for (Int_t yIndex=yMin; yIndex<=yMax; yIndex++){
      cout <<"DEBUG::0.1" << endl;

      if (centIndex==1 && yIndex == 18) continue;

      canvas->cd(centIndex+1);

      spectra = NULL;
      scaledSpectra = NULL;
      spectraFit = NULL;


      TString spectraName(Form("CorrectedSpectra_ProtonPlus_Cent%02d_yIndex%02d_Total",
                               type.Data(),centIndex,yIndex));

      TString spectraFitName(Form("CorrectedSpectra_ProtonPlus_Cent%02d_yIndex%02d_Total_Fit0",
				  centIndex,yIndex));

      TString spectraFitName1(Form("CorrectedSpectra_ProtonPlus_Cent%02d_yIndex%02d_Total_Fit1",
                                  centIndex,yIndex));

      spectra = (TGraphErrors *)file->Get(Form("%sSpectra_ProtonPlus/%s",
					       Type.Data(),
					       spectraName.Data()));

      
      spectraFit = (TF1 *)file->Get(Form("ThermalFit_%s/%s",
					 speciesName.Data(),
					 spectraFitName.Data()));

      cout << Form("ThermalFit_%s/%s", speciesName.Data(), spectraFitName.Data()) << endl;

      cout << "DEBUG::1.4" << endl;
      if (spectra) cout << "Spectra exists!" << endl;
      else continue;
      if (spectraFit) cout << "Fit exists!" << endl;
      else continue;
      //Skip if no spectrum was found for this rapidity bin
      if (spectra == NULL)
        continue;

      cout << "DEBUG::2" << endl;
      //Make sure spectra have points
      if (spectra->GetN() == 0){
        delete spectra;
        continue;
      }

      //Remove high mTm0 Point
      TGraphChop(spectra,maxmTm0,false);

      //Remove points with large errors
      //      RemovePointsWithLargeErrors(spectra);

      //Scale the Spectra
      scaledSpectra = TGraphScale(spectra,pow(scaleFactor,yIndex-midRapidityIndex));

      if (centIndex == 1){

	if (yIndex == 9 ) {
	  //  scaledSpectra->RemovePoint(0);
	  //  scaledSpectra->RemovePoint(0);
	}
	//if (yIndex == 8 ) scaledSpectra->RemovePoint(0);
      }

      cout << "DEBUG::3" << endl;
      //Scale the Spectra Fit
      Double_t fitMin(0), fitMax(0);

      if (corrected && spectraFit){

	TGraphErrors * spectraGraph = scaledSpectra;
	cout << "DEBUG::3.1" << endl;

	//  Double_t *xarr = spectraGraph->GetX();
	int nPoints = spectraGraph->GetN();
	cout << "points " << nPoints << endl;
	if (nPoints > 0 ){
	    
	  Double_t xMin = TMath::MinElement(spectraGraph->GetN(),spectraGraph->GetX());
	  Double_t xMax = TMath::MaxElement(spectraGraph->GetN(),spectraGraph->GetX());
	  //Double_t xMin = TMath::MinElement(spectraGraph->GetN(),spectraGraph->GetX());
	  //Double_t xMax = TMath::MaxElement(spectraGraph->GetN(),spectraGraph->GetX());

	  cout << "DEBUG::3.2" << endl;
	  Int_t iFuncBin = centIndex*100 + yIndex;

	  //  funcToScale[centIndex][iFuncBin] = spectraFit;
	  spectraFit->SetRange(xMin,xMax);
	  spectraFit->GetRange(fitMin,fitMax);

	  cout << "DEBUG::3.3" << endl;
	  TF1 * scaledSpectraFit = new TF1("scaledFunc",ThermalFitFunc,xMin,xMax,3);
	  TF1 * scaledSpectraFitExp = new TF1("scaledExpFunc",ThermalFitFunc,0.01,xMax,3);

	  scaledSpectraFit->SetParameter(0,spectraFit->GetParameter(0)*pow(scaleFactor,yIndex-midRapidityIndex));
	  scaledSpectraFit->SetParameter(1,spectraFit->GetParameter(1));
	  scaledSpectraFit->SetParameter(2,spectraFit->GetParameter(2));
	  scaledSpectraFitExp->SetParameter(0,spectraFit->GetParameter(0)*pow(scaleFactor,yIndex-midRapidityIndex));
	  scaledSpectraFitExp->SetParameter(1,spectraFit->GetParameter(1));
	  scaledSpectraFitExp->SetParameter(2,spectraFit->GetParameter(2));


	  //  cout << "Spectra Fit: " << spectraFit->Print() <[yIndex]< endl; 
	  cout << "DEBUG::4" << endl;
	  
	  cout << "DEBUG::5" << endl;
	  //Create the Label
	  Double_t yLocation =  spectraGraph->GetY()[spectraGraph->GetN()-1];
	  Double_t xLocation = spectraGraph->GetX()[spectraGraph->GetN()-1]+.05;
	  cout << "xLoc= " << xLocation << endl;
	  if (yIndex < 10 ) xLocation = 1.025;
	  cout << "DEBUG::6" << endl;

	  label = new TPaveText(xLocation,yLocation,xLocation+.2,yLocation,"NB");
	  label->SetTextSize(0.030);
	  label->AddText(Form("y = %.1f",GetRapidityRangeCenter(yIndex)/-1.0));
	  cout << "DEBUG::7" << endl;
	  //Set the Mid Rapidity Spectrum to Red
	  if (yIndex == midRapidityIndex){
	    scaledSpectra->SetMarkerColor(kBlack);
	    label->SetTextColor(kRed);
	    if (corrected && scaledSpectraFit){
	      scaledSpectraFit->SetLineColor(kBlack);
	      scaledSpectraFit->SetLineWidth(2);
	      scaledSpectraFitExp->SetLineColor(kBlack);
	      scaledSpectraFitExp->SetLineWidth(2);
	      scaledSpectraFitExp->SetLineStyle(7);
	    }
	  }
	  else{
	    scaledSpectra->SetMarkerColor(kBlack);
	    if (corrected && scaledSpectraFit){
	      scaledSpectraFit->SetLineColor(kBlack);
	      scaledSpectraFit->SetLineWidth(2);
	      scaledSpectraFitExp->SetLineColor(kBlack);
	      scaledSpectraFitExp->SetLineWidth(2);
	      scaledSpectraFitExp->SetLineStyle(7);
	    }
	  }
	  cout << "DEBUG::8" << endl;
	  if (corrected) {
	    scaledSpectraFit->Draw("SAME");
	    scaledSpectraFitExp->Draw("SAME");
	  }
	  scaledSpectra->SetMarkerSize(1.3);
	  scaledSpectra->SetMarkerStyle(kFullCircle);
	  scaledSpectra->Draw("P");

	  label->Draw("SAME");

	}
      }//End Loop Over Spectra
    }

    cout << "DEBUG::9" << endl;
    gPad->Update();
    gPad->Modified();
  }


  cout << "DEBUG::7..1" << endl;
  //Create the Title

 
  TMarker *tpcMarker = new TMarker(0,0,kFullCircle);
  tpcMarker->SetMarkerColor(kGray+2);
  tpcMarker->SetMarkerSize(1.5);

  TMarker *tofMarker = new TMarker(0,0,kFullSquare);
  tofMarker->SetMarkerColor(kGray+2);
  tofMarker->SetMarkerSize(1.5);

  cout << "DEBUG::7..4" << endl;

  TLegend *leg1 = new TLegend(.67,.72,.89,.76);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(kWhite);
  leg1->SetTextSize(.035);
  leg1->SetNColumns(2);
  //leg1->AddEntry(tpcMarker,"TPC","P");
  //leg1->AddEntry(tofMarker,"TOF","P");
  //leg1->Draw("SAME");

  cout << "DEBUG::7..5" << endl;

  //    if (save){
  //      canvas->Print(Form("../images/gif/%s.gif",canvas->GetName()));
  //      canvas->Print(Form("../images/ps/%s.ps",canvas->GetName()));
  //  }


  for (int i =1; i <=6;i++){

    canvas->cd(i);
    //    TPaveText *title = new TPaveText(.56,.82,.91,.92,"NBNDCBR");
    TPaveText *centralityTitle = new TPaveText(.5,.82,.91,.92,"BRNBNDC");
    centralityTitle->SetFillColor(kWhite);
    centralityTitle->SetBorderSize(0);
    centralityTitle->SetTextFont(30);
    centralityTitle->SetTextSize(25);
    centralityTitle->AddText(Form("%d-%d%% Central",
				  i*5,
				  (i+1)*5));
 
    //  centralityTitle->Draw();
    /* 
 TPaveText * starPrelim = new TPaveText(.55,.13,.90,.19,"BRNBNDC");
 starPrelim->SetFillColor(kWhite);
 starPrelim->SetBorderSize(0);
 starPrelim->SetTextFont(30);
 starPrelim->SetTextSize(18);
 starPrelim->AddText("STAR PRELIMINARY");
 starPrelim->Draw("SAME");
    */

  }

  canvas->SaveAs("proton_AlAu_Spectra.pdf");
  canvas->SaveAs("proton_AlAu_Spectra.eps");


}
