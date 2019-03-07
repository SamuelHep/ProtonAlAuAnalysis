

Double_t rapidityMin = -2.05;
Double_t rapidityMax = 2.05;
Double_t rapidityBinWidth =  0.1;

Double_t mTm0Min = 0.0;
Double_t mTm0Max = 2.0;
Double_t mTm0BinWidth = 0.025;


const int nRapidityBins = 41;
const int nmTm0Bins = 80;

TF1 *funcToScale[nRapidityBins];

//__________________________________________________________________________

Double_t scaleFunction(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  Double_t scale = par[0];
  Int_t yIndex = (int)par[1];

  return funcToScale[yIndex]->Eval(xx) * scale;

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



void DrawProtonSpectraNoFit(TString spectraFile, TString imageDir, TString eventConfig, TString system, Double_t energy,
			    Int_t centIndex,
		       Bool_t corrected=true){

  Bool_t save = true;        //Should the canvas be saved?
  Double_t scaleFactor = 3.0; //The factor to scale the non- midrapidity spectra by
  Int_t yMin = 7;            //Minimum rapidity bin index to draw //Collider 11, PosY 18, NegY
  Int_t yMax = 19;            //Maximum rapidity bin index to draw //Collider 29, PosY 33, NegY
  Double_t maxmTm0=1.0;       //Maximum mT-m0 value to draw
  Double_t midY =-1.5;

  Double_t yScaleMin(0), yScaleMax(0);
  
  yScaleMin = .05;
  yScaleMax = 1000000000;

  

  //Load the root file
  TFile *file = new TFile(spectraFile,"READ");

  //Create the Spectra Name
  TString type = "corrected";
  TString Type = "Corrected";
  TString name = "Corrected";
 
  TString speciesName = "ProtonPlus";
  Int_t midRapidityIndex = 5;

  //Create Array of TGraphErrors
  TGraphErrors *spectra[nRapidityBins];
  TGraphErrors *scaledSpectra[nRapidityBins];

  TPaveText *label[nRapidityBins];

  //Loop Over the centrality Bins

    TCanvas *canvas = NULL;

    

    //Create the Canvas and Frame
    canvas = new TCanvas(Form("%s_%s_%s_Cent%02d",
                              system.Data(),speciesName.Data(),
                              eventConfig.Data(),centIndex),"title",0,20,600,800);

    canvas->SetTopMargin(.05);
    canvas->SetRightMargin(.05);
    canvas->SetLeftMargin(.17);
    canvas->SetLogy();
    canvas->SetTicks(1,1);

    
    //    TH1F *frame = canvas->DrawFrame(0,yScaleMin,maxmTm0+.4,yScaleMax);
    TH1F *frame = canvas->DrawFrame(0,yScaleMin,maxmTm0+.3,yScaleMax);
    frame->GetXaxis()->SetTitle("m_{T}-m_{p} (GeV/c^{2})");
    frame->GetYaxis()->SetTitle("#frac{1}{N_{Evt}}#times#frac{1}{2#pim_{T}}#times#frac{d^{2}N}{dm_{T}dy} (GeV/c^{2})^{-2}");
    frame->GetYaxis()->SetTitleOffset(2.2);
    frame->GetYaxis()->SetTitleFont(63);
    frame->GetYaxis()->SetTitleSize(25);
    frame->GetXaxis()->SetTitleFont(63);
    frame->GetXaxis()->SetTitleSize(25);
    frame->GetXaxis()->SetTitleOffset(1.3);

    canvas->Update();

    //Loop Over the rapidity bins
    for (Int_t yIndex=yMin; yIndex<=yMax; yIndex++){

      spectra[yIndex] = NULL;
      scaledSpectra[yIndex] = NULL;

      TString spectraName(Form("%sSpectra_ProtonPlus_Cent%02d_yIndex%02d",
                               type.Data(),centIndex,yIndex));


      spectra[yIndex] = (TGraphErrors *)file->Get(Form("%sSpectra_ProtonPlus/%s",
                                                       Type.Data(),
                                                       spectraName.Data()));



      if (spectra[yIndex]) cout << "Spectra exists!" << endl;

      //Skip if no spectrum was found for this rapidity bin
      if (spectra[yIndex] == NULL)
        continue;

      //Make sure spectra have points
      if (spectra[yIndex]->GetN() == 0){
	//        delete spectra[yIndex];
        continue;
      }

      //Remove high mTm0 Point
      TGraphChop(spectra[yIndex],maxmTm0,false);

      //Remove points with large errors
      RemovePointsWithLargeErrors(spectra[yIndex]);

      //Scale the Spectra
      scaledSpectra[yIndex] = TGraphScale(spectra[yIndex],pow(scaleFactor,yIndex-midRapidityIndex));

      //Scale the Spectra Fit
      Double_t fitMin(0), fitMax(0);
      if (yIndex < 20){
	if (corrected){

	  Double_t *xarr = scaledSpectra[yIndex]->GetX();

	  int nPoints = scaledSpectra[yIndex]->GetN();

	  if (nPoints < 1) continue;

	  Double_t xMin = TMath::MinElement(scaledSpectra[yIndex]->GetN(),scaledSpectra[yIndex]->GetX());

	  Double_t xMax = TMath::MaxElement(scaledSpectra[yIndex]->GetN(),scaledSpectra[yIndex]->GetX());



	  //	  scaledSpectraFit[yIndex] = new TF1(Form("%s_Scaled",spectraFit[yIndex]->GetName()),
	  //				     scaleFunction,fitMin,fitMax+.01,2);
	  //scaledSpectraFit[yIndex]->SetParameter(0,pow(scaleFactor,yIndex-midRapidityIndex));
	  //scaledSpectraFit[yIndex]->SetParameter(1,yIndex);

	  //scaledSpectraFitExp[yIndex] = new TF1(Form("%s_Scaled_Exp",spectraFit[yIndex]->GetName()),
	  //						scaleFunction,0,fitMax+.01,2);
	//scaledSpectraFitExp[yIndex]->SetParameter(0,pow(scaleFactor,yIndex-midRapidityIndex));
	//scaledSpectraFitExp[yIndex]->SetParameter(1,yIndex);

	}

	//Create the Label
	Double_t yLocation = scaledSpectra[yIndex]->GetY()[scaledSpectra[yIndex]->GetN()-1];
	Double_t xLocation = scaledSpectra[7]->GetX()[scaledSpectra[7]->GetN()-1]+.05;


	//	if (yIndex > 14) xLocation = scaledSpectra[yIndex]->GetX()[scaledSpectra[yIndex]->GetN()-1]+.1;


	label[yIndex] = new TPaveText(xLocation,yLocation,xLocation+.2,yLocation,"NB");
	label[yIndex]->SetTextSize(0.018);
	label[yIndex]->AddText(Form("y = %.1f",GetRapidityRangeCenter(yIndex)));

	//Set the Mid Rapidity Spectrum to Red
	//	if (yIndex == midRapidityIndex){
	  scaledSpectra[yIndex]->SetMarkerColor(kBlack);
	  // label[yIndex]->SetTextColor(kRed);
	  //	  if (corrected ){ //&& scaledSpectraFit[yIndex]){
	    //	    scaledSpectraFit[yIndex]->SetLineColor(kBlack);
	    //	    scaledSpectraFit[yIndex]->SetLineWidth(2);
	    //	    scaledSpectraFitExp[yIndex]->SetLineColor(kBlack);
	    //	    scaledSpectraFitExp[yIndex]->SetLineWidth(2);
	    //	    scaledSpectraFitExp[yIndex]->SetLineStyle(7);

       
	//	else{
	  scaledSpectra[yIndex]->SetMarkerColor(kBlack);
	  // if (corrected && scaledSpectraFit[yIndex]){
	  //scaledSpectraFit[yIndex]->SetLineColor(kBlack);
	  //	    scaledSpectraFit[yIndex]->SetLineWidth(2);
	  //	    scaledSpectraFitExp[yIndex]->SetLineColor(kBlack);
	  //	    scaledSpectraFitExp[yIndex]->SetLineWidth(2);
	  //	    scaledSpectraFitExp[yIndex]->SetLineStyle(7);
	  //	  }
	  //	}

	if (corrected) {
	  //	  scaledSpectraFit[yIndex]->Draw("SAME");
	  //	  scaledSpectraFitExp[yIndex]->Draw("SAME");
	}
	scaledSpectra[yIndex]->SetMarkerSize(1.3);
	scaledSpectra[yIndex]->SetMarkerStyle(kFullStar);
	scaledSpectra[yIndex]->Draw("P");

	label[yIndex]->Draw("SAME");


      }//End Loop Over Spectra
    }

    //Create the Title
    TPaveText *title = new TPaveText(.56,.82,.91,.92,"NBNDCBR");
    title->SetFillColor(kWhite);
    title->SetBorderSize(0);
    title->SetTextSize(.05);
    title->SetTextAlign(32);
    title->AddText(Form("p Spectra %s",
                        system.Data()));
    title->AddText(Form("#sqrt{s_{NN}} = %.03g GeV",energy));
    
    title->GetLine(1)->SetTextSize(.05);
    
    title->Draw("SAME");

    TMarker *midYMarker = new TMarker(0,0,kFullCircle);
    midYMarker->SetMarkerColor(kRed);
    midYMarker->SetMarkerSize(1.2);

    TMarker *notMidYMarker = new TMarker(0,0,kFullStar);
    notMidYMarker->SetMarkerColor(kBlack);
    notMidYMarker->SetMarkerSize(2);

    TLegend *leg = new TLegend(.65,.76,.92,.82);
    leg->SetLineColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(.035);
    
    leg->AddEntry(notMidYMarker,Form("x %g^{#pm n}",scaleFactor),"P");
    leg->Draw("SAME");

    TPaveText *centralityTitle = new TPaveText(.2,.13,.47,.19,"BRNBNDC");
    centralityTitle->SetFillColor(kWhite);
    centralityTitle->SetBorderSize(0);
    centralityTitle->SetTextFont(63);
    centralityTitle->SetTextSize(25);
    centralityTitle->AddText(Form("%d-%d%% Central",
                                  centIndex*5,
                                  (centIndex+1)*5));
    
    centralityTitle->Draw("SAME");

    TPaveText *starPrelim = new TPaveText(.55,.13,.90,.19,"BRNBNDC");
    starPrelim->SetFillColor(kWhite);
    starPrelim->SetBorderSize(0);
    starPrelim->SetTextFont(63);
    starPrelim->SetTextSize(18);
    starPrelim->AddText("STAR PRELIMINARY");
    starPrelim->Draw("SAME");

    TMarker *tpcMarker = new TMarker(0,0,kFullCircle);
    tpcMarker->SetMarkerColor(kGray+2);
    tpcMarker->SetMarkerSize(1.5);

    TMarker *tofMarker = new TMarker(0,0,kFullSquare);
    tofMarker->SetMarkerColor(kGray+2);
    tofMarker->SetMarkerSize(1.5);

    TLegend *leg1 = new TLegend(.67,.72,.89,.76);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(kWhite);
    leg1->SetTextSize(.035);
    leg1->SetNColumns(2);
    //leg1->AddEntry(tpcMarker,"TPC","P");
    //leg1->AddEntry(tofMarker,"TOF","P");
    //leg1->Draw("SAME");


    if (save){
      canvas->Print(Form("%s/%s.gif",imageDir.Data(),canvas->GetName()));
      canvas->Print(Form("%s/%s.ps",imageDir.Data(),canvas->GetName()));
      canvas->Print(Form("%s/%s.eps",imageDir.Data(),canvas->GetName()));
    }


}
