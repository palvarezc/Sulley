#include<iostream>
#include<TMath.h>
#include<TFile.h>
#include<TTree.h>
#include<RooRealVar.h>
#include<RooArgSet.h>
#include<RooDataSet.h>
#include<RooKeysPdf.h>
#include<RooPlot.h>
#include<TCanvas.h>
#include<RooGlobalFunc.h>
#include<RooHistPdf.h>
#include<RooDataHist.h>
#include<RooAddPdf.h>
#include<RooArgList.h>
#include<TH1F.h>
#include<TH2F.h>
#include<RooNDKeysPdf.h>
#include<RooWorkspace.h>
#include<string>
#include<vector>
#include<RooExponential.h>
#include<TRandom.h>
#include<RooGaussian.h>
#include<RooMinuit.h>
#include<fstream>
#include<iomanip>
#include<RooFitResult.h>
#include<RooChebychev.h>
#include<RooProdPdf.h>
#include "RooMcorMvisTsallis.h"
#include "usefulFunctions.h"


using namespace std;
using namespace RooFit;

void toystudy1DMakeFit(string workspacename,  bool wantplot, 
      RooHistPdf& histPdfSignalZeroGamma, RooHistPdf& histPdfSignalOneGamma, RooHistPdf& histPdfSignalTwoGamma, 
                       RooHistPdf& histPdfPartReco, RooMcorMvisTsallis& McorMvis, RooMcorMvisTsallis& McorMvis_fit, 
      RooAbsPdf::GenSpec& GenSpecSigZeroGamma, RooAbsPdf::GenSpec& GenSpecSigOneGamma, RooAbsPdf::GenSpec& GenSpecSigTwoGamma, 
      RooAbsPdf::GenSpec& GenSpecPartReco, RooAbsPdf::GenSpec& GenSpecComb, int nGenSignalZeroGamma, int nGenSignalOneGamma,
      int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, double expoConstGen, RooRealVar& expoConst, ofstream& out, TTree* t, bool update);

void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartRec, int nGenComb, double expoConstGen, RooRealVar& nSignal,
      RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst);


