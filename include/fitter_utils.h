#ifndef FITTER_UTILS_H
#define FITTER_UTILS_H

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
#include<RooAbsPdf.h>
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


void initiateParams(RooArgSet& parset);

void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartRec, int nGenComb, double expoConstGen, RooRealVar& nSignal,
      RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst);

void prepare_PDFs(string workspacename, string trigStr, string BDTcut, bool fit2D,
                  string signalfile, string partrecofile, string combinatorialfile,
                  string signaltree = "DecayTree", string partrecotree = "DecayTree", string combinatorialtree = "DecayTree");

void plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataSet dataGenTot);



void generate_and_fit(string workspacename,  bool fit2D, bool wantplot, bool constPartReco,
              int nGenSignal,  int nGenPartReco,  int nGenComb,
              double nGenFracZeroGamma,  double nGenFracOneGamma,
              ofstream& out, TTree* t, bool update, string plotsfile);



#endif
