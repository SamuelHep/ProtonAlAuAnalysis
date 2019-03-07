#include <iostream>
#include <vector>
#include <utility>

#include <TCanvas.h>
#include <TString.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TMath.h>

#include "inc/ReturnStructs.h"

void MakeResultsTable(NegBinomialSearchResults *results){

  TCanvas *canvas = new TCanvas("resultsTableCanvas","resultsTableCanvas",20,20,700,900);

  canvas->SetTopMargin(.05);
  canvas->SetBottomMargin(.05);
  canvas->SetRightMargin(.05);
  canvas->SetLeftMargin(.05);

  TPaveText *title = new TPaveText(.1,.81,.9,.86,"brNDC");
  title->SetFillColor(kWhite);
  title->SetBorderSize(0);
  title->SetTextSize(.045);
  title->SetTextFont(62);
  title->AddText("Centrality Bin Results");
  title->Draw("SAME");
  
  TLine *horizLine = new TLine();
  horizLine->DrawLineNDC(.1,.8,.94,.8);

  TLine *horizLine1 = new TLine(*horizLine);
  horizLine1->DrawLineNDC(.1,.75,.94,.75);

  const int nCols(6);
  
  //Measures
  Double_t colStart[nCols] = {.1,.21,.37,.51,.66,.8};
  Double_t colStop[nCols]  = {.2,.36,.50,.65,.79,.94};

  Double_t rowHeight = .04;
		    
  //Column Title Bar
  TPaveText *colTitle[nCols];
  TLine *vertLines[nCols];
  TString colNames[nCols] = {"Index","% Central","N_{ch} Cut","<N_{part}>","<N_{coll}>","<b>"};
  for (Int_t iCol=0; iCol<nCols; iCol++){

    colTitle[iCol] = new TPaveText(colStart[iCol],.755,colStop[iCol],.785,"NDC");
    colTitle[iCol]->SetBorderSize(0);
    colTitle[iCol]->SetFillColor(kWhite);
    colTitle[iCol]->SetTextFont(42);
    colTitle[iCol]->SetTextSize(.03);
    colTitle[iCol]->AddText(colNames[iCol]);
    colTitle[iCol]->Draw("SAME");

    if (iCol != nCols-1){
      vertLines[iCol] = new TLine();
      Double_t locx = (colStop[iCol]+colStart[iCol+1])/2.0;
      vertLines[iCol]->DrawLineNDC(locx,.8,locx,.1);

    }

  }//End Loop to make title bar
  
  //Make all the Entries
  const int nRows = 16;//results->nCentralityBins;
  const int nEntries = nCols*nRows;
  TPaveText *entries[nEntries];
  int iEntry(0);
  for (Int_t iRow=0; iRow<nRows; iRow++){

    Double_t rowStart = .74 - (iRow+1)*rowHeight;
    
    for (Int_t iCol=0; iCol<nCols; iCol++){

      entries[iEntry] = new TPaveText(colStart[iCol],rowStart,colStop[iCol],rowStart+rowHeight,"NDC");
      entries[iEntry]->SetBorderSize(0);
      entries[iEntry]->SetTextSize(.03);
      entries[iEntry]->SetTextFont(62);
      entries[iEntry]->SetFillColor(kWhite);
          
      if (iCol == 0)
	entries[iEntry]->AddText(Form("%d",iRow));	
      else if (iCol == 1)
	entries[iEntry]->AddText(Form("%d-%d",iRow == nRows-1 ? 0 : TMath::Nint(results->centralityBinDefinitions.at(iRow+1)*100),
				      TMath::Nint(results->centralityBinDefinitions.at(iRow)*100)));
      else if (iCol == 2)
	entries[iEntry]->AddText(Form("%d",TMath::Nint(results->centralityBinCuts.at(iRow))));
      else if (iCol == 3)
	entries[iEntry]->AddText(Form("%d #pm %d",TMath::Nint(results->nPartMeans.at(iRow)),
				      TMath::Nint(sqrt(pow(results->nPartStatErrors.at(iRow),2)+
						       pow(results->nPartSysErrors.at(iRow),2)))));
      else if (iCol == 4)
	entries[iEntry]->AddText(Form("%d #pm %d",TMath::Nint(results->nCollMeans.at(iRow)),
				      TMath::Nint(sqrt(pow(results->nCollStatErrors.at(iRow),2)+
						       pow(results->nCollSysErrors.at(iRow),2)))));
      else if (iCol == 5)
	entries[iEntry]->AddText(Form("%.01F #pm %.01F",results->impactParamMeans.at(iRow),
				      sqrt(pow(results->impactParamStatErrors.at(iRow),2)+
						       pow(results->impactParamSysErrors.at(iRow),2))));
      
      entries[iEntry]->Draw("SAME");
      
      iEntry++;
    }//End Loop Over Column

  }//End Loop Over Row

  results->resultsTableCanvas = (TCanvas *)canvas->Clone(canvas->GetName());

  delete canvas;
}
