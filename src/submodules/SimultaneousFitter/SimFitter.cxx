#include <iostream>
#include <map>
#include <vector>

#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include <Foption.h>
//#include <Math/Minimizer.h>

#include <SimFitter.h>


//_____________________________________________________________________
ClassImp(SimFitter);
SimFitter::SimFitter(std::vector<TGraphErrors *> objs, std::vector<TF1 *> fits){

  //Check to make sure the vectors are the same size
  if (objs.size() != fits.size()){
    cout <<"ERROR - SimFitter::SimFitter() - The vector arguments must be the same size. Exiting!\n";
    exit (EXIT_FAILURE);
  }
  
  //Transfer Constructor Args to Internal Quantities
  objectsToFit = objs;
  nObjects = objectsToFit.size();
  isInitialized = false;
  dataSize = 0;
  printLevel = 0;
    
  //Create Necessary Internal Data Structures
  parameterMap.resize(nObjects);
  
  
  for (unsigned int i=0; i<nObjects; i++){

    fitFuncResults.push_back(new TF1 (*fits.at(i)));
    fitFuncResults.at(i)->SetName(Form("Fit%d",i));
    
    fitFuncs.push_back(new ROOT::Math::WrappedMultiTF1(*fitFuncResults.at(i),1));
    dataRange.push_back(new ROOT::Fit::DataRange());
    
    dataOption.push_back(new ROOT::Fit::DataOptions());
    dataOption.at(i)->fUseRange = true;

  }

}

//____________________________________________________________________
SimFitter::~SimFitter(){

  for (unsigned int iObj=0; iObj<nObjects; iObj++){
    delete fitFuncs.at(iObj);
    delete dataOption.at(iObj);
    delete dataRange.at(iObj);
    delete binData.at(iObj);
    delete chi2Func.at(iObj);
  }

  objectsToFit.clear();
  fitFuncs.clear();
  dataOption.clear();
  dataRange.clear();
  binData.clear();
  chi2Func.clear();

  delete globalChi2;
}

//_____________________________________________________________________
void SimFitter::Init(std::vector<double> parVals, std::vector<std::pair<float,float> > *range){

  if (range && range->size() != dataRange.size()){
    std::cout <<"ERROR - SimFitter::SetFitRange() - The size of the vector range must be the same as the number of objects to fit. EXITING!\n";
    exit (EXIT_FAILURE);
  }

  totalFitPars = parVals;
  globalChi2 = new GlobalChi2(this);

  for (unsigned int i=0; i<nObjects; i++){

    if (range)
      dataRange.at(i)->SetRange(range->at(i).first,range->at(i).second);
    else
      dataRange.at(i)->SetRange(TMath::MinElement(objectsToFit.at(i)->GetN(),
						  objectsToFit.at(i)->GetX()),
				TMath::MaxElement(objectsToFit.at(i)->GetN(),
						  objectsToFit.at(i)->GetX()));

    

    binData.push_back(new ROOT::Fit::BinData(*dataOption.at(i),*dataRange.at(i)));
    ROOT::Fit::FillData(*binData.at(i), objectsToFit.at(i));
    dataSize += binData.back()->Size();

    chi2Func.push_back(new ROOT::Fit::Chi2Function(*binData.at(i), *fitFuncs.at(i)));
        
  }

  fitter.Config().SetParamsSettings((int)totalFitPars.size(),&totalFitPars.at(0));
  fitter.Config().MinimizerOptions().SetPrintLevel(printLevel);
  fitter.Config().SetMinimizer("Minuit2","MigradImproved");
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(10000);

  //fitter.Config().MinosErrors();

  isInitialized = true;
  
}

//_____________________________________________________________________
void SimFitter::SetParameterRelationships(std::vector<std::vector<int> > parRelations){

  parameterRelationships = parRelations;
  
  for (unsigned int iFunc=0; iFunc<nObjects; iFunc++){
    for (unsigned int iPar=0; iPar<parRelations.at(iFunc).size(); iPar++){
      parameterMap.at(iFunc)[iPar] = parRelations.at(iFunc).at(iPar);
    }
  }
  
  
}

//_____________________________________________________________________
Double_t SimFitter::ComputeGlobalChi2(const Double_t *pars){
  
  for (unsigned int iPar=0; iPar<totalFitPars.size(); iPar++){
    totalFitPars.at(iPar) = pars[iPar];
  }

  double chi2(0);
  for (unsigned int i=0; i<nObjects; i++){

    const int nPars = fitFuncResults.at(i)->GetNpar();
    double vals[nPars];
    for (int iPar=0; iPar<nPars; iPar++){
      vals[iPar] = totalFitPars.at(parameterMap.at(i)[iPar]);
    }

    chi2 += (*chi2Func.at(i))(vals);
    
  }

  return chi2;
}

//_____________________________________________________________________
ROOT::Fit::FitResult SimFitter::Fit(){

  if (!isInitialized){
    cout <<"ERROR - SimFitter::Fit() - You must call SimFitter::Init before calling Fit. EXITING.\n";
    exit (EXIT_FAILURE);
  }

  fitter.FitFCN(totalFitPars.size(),*globalChi2,0,dataSize,true);

  fitResult = fitter.Result();
  fitResult.Print(std::cout);

  for (unsigned int iFunc=0; iFunc<nObjects; iFunc++){
    for (int iPar=0; iPar<fitFuncResults.at(iFunc)->GetNpar(); iPar++){

      fitFuncResults.at(iFunc)->SetParameter(iPar,fitResult.Parameter(parameterRelationships.at(iFunc).at(iPar)));
      fitFuncResults.at(iFunc)->SetParError(iPar,fitResult.ParError(parameterRelationships.at(iFunc).at(iPar)));
      
    }
  }

  ConstructCovarianceMatrix();
  
  return fitResult;

}

//______________________________________________________________________
void SimFitter::ConstructCovarianceMatrix(){

  for (unsigned int iFit=0; iFit<nObjects; iFit++){

    int nPars = fitFuncResults.at(iFit)->GetNpar();
    covarianceMatrix.push_back(new TMatrixDSym(nPars));

    double matrixArray[nPars*nPars];
    for (int i=0; i<nPars; i++){
      for (int j=0; j<nPars; j++){

	matrixArray[i*nPars+j] = fitResult.CovMatrix(parameterMap.at(iFit)[i],
						      parameterMap.at(iFit)[j]);
	
      }//End Loop Over j
    }//End Loop Over i

    covarianceMatrix.back()->SetMatrixArray(matrixArray);
    
  }//End Loop Over Fits

}


//______________________________________________________________________
TMatrixD SimFitter::GetRandomDeviateVector(int iFit){

  int dim = fitFuncResults.at(iFit)->GetNpar();

  TMatrixDSym covMatrix = *covarianceMatrix.at(iFit);
  TDecompChol decompMatrix(covMatrix);
  decompMatrix.Decompose();

  TMatrixD upper = decompMatrix.GetU();

  TMatrixD A(dim,dim);
  A.Transpose(upper);

  TMatrixD meanVector(dim,1);
  for (int iDim=0; iDim<dim; iDim++)
    meanVector[iDim][0] = fitFuncResults.at(iFit)->GetParameter(iDim);

  TMatrixD indepDeviates(dim,1);
  TRandom3 rand(0);
  for (int iDim=0; iDim<dim; iDim++)
    indepDeviates[iDim][0] = rand.Gaus(0,1);

  TMatrixD deviateVector(dim,1);
  deviateVector = meanVector + A*indepDeviates;

  return deviateVector;
}

//______________________________________________________________________
TGraphErrors *SimFitter::GetConfidenceInterval(int iFit, float cl){

  if (iFit < 0 || iFit >= (int)fitFuncResults.size()){
    cout <<"ERROR - SimFitter::GetConfidenceInterval() - Argument iFit is outside range. EXITING!\n";
    exit (EXIT_FAILURE);
  }
  
  TGraphErrors *gr = new TGraphErrors();
  gr->SetName(Form("%s_Conf",fitFuncResults.at(iFit)->GetName()));
  gr->SetFillColor(fitFuncResults.at(iFit)->GetLineColor()-3);
  gr->SetFillStyle(3002);

  const int nPar = fitFuncResults.at(iFit)->GetNpar();
  const int npx = fitFuncResults.at(iFit)->GetNpx();
  Double_t fitMin(0), fitMax(0);
  fitFuncResults.at(iFit)->GetRange(fitMin,fitMax);
  Double_t stepSize = fabs(fitMax-fitMin)/(double)npx;

  const int nIters(100);
  double vals[npx][nIters];
  for (int iIter=0; iIter<nIters; iIter++){

    TMatrixD randomDeviates = GetRandomDeviateVector(iFit);

    TF1 tempFunc(*fitFuncResults.at(iFit));
    for (int iPar=0; iPar<nPar; iPar++){
      tempFunc.FixParameter(iPar,randomDeviates[iPar][0]);
    }
    
    for (int iPoint=0; iPoint<npx; iPoint++){
      vals[iPoint][iIter] = tempFunc.Eval(iPoint * stepSize);
    }
    
  }//End Loop Over Iters

  for (int i=0; i<npx; i++){
    gr->SetPoint(gr->GetN(),i*stepSize,TMath::Mean(npx,vals[i]));
    gr->SetPointError(gr->GetN()-1,cl*TMath::RMS(npx,vals[i]));
  }

  
  /*
  for (int i=0; i<npx; i++){
    gr->SetPoint(i,fitMin+i*stepSize,fitFuncResults.at(iFit)->Eval(fitMin+i*stepSize));
  }
  */
  
  //fitResult.SetModelFunction(fitFuncs.at(iFit));
  //fitResult.GetConfidenceIntervals(gr->GetN(),1,1,gr->GetX(),gr->GetEY(),cl);

  return gr;
}

//_________________________________________________________________________
std::vector<TGraphErrors *> SimFitter::GetConfidenceIntervals(float cl){

  std::vector<TGraphErrors *> confIntervals;

  for (unsigned int iFit=0; iFit<nObjects; iFit++){
    confIntervals.push_back(GetConfidenceInterval(iFit,cl));
  }

  return confIntervals;
  
}
