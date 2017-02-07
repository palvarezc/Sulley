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
#include "RooNumIntConfig.h"

using namespace std;
using namespace RooFit;

class FitterUtils
{

   public:

   FitterUtils(int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_, double nGenFracOneGamma_, bool fit2D_, string workspacename_);

   void prepare_PDFs(string trigStr, string BDTVar, double BDTcut,
         string signalfile, string partrecofile, string combinatorialfile, string JpsiLeakfile,
         double minBMass = 4880, double maxBMass = 5700,
         string signaltree = "DecayTree", string partrecotree = "DecayTree", string combinatorialtree = "DecayTree", string JpsiLeaktree = "DecayTree");


   void generate();
   void fit(bool wantplot, bool constPartReco,
         double fracPartReco_const, ofstream& out, TTree* t, bool update, string plotsfile);


   protected:

   void initiateParams(RooArgSet* parset);

   void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, RooRealVar const& expoConstGen,
         RooRealVar& nSignal, RooRealVar& nPartReco, 
         RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma);

   void plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataSet dataGenTot);

   void PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& B_plus_M_corr);
   void PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& B_plus_M_corr);
   void PlotShape1D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M);


   int nGenSignal;
   int nGenPartReco;
   int nGenComb;
   int nGenJpsiLeak;
   double nGenFracZeroGamma;  
   double nGenFracOneGamma;
   bool fit2D;
   string workspacename;

};


#endif
