#ifndef FITTER_UTILS_EXPOFPOLY_H
#define FITTER_UTILS_EXPOFPOLY_H

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
#include "RooPTMVis.h"
#include "usefulFunctions.h"
#include "RooNumIntConfig.h"

using namespace std;
using namespace RooFit;

class FitterUtilsExpOfPolyTimesX
{

   public:

   FitterUtilsExpOfPolyTimesX(int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_, double nGenFracOneGamma_, string workspacename_);

  void prepare_PDFs(string trigStr, string weightStr, string BDTVar, double BDTcut,
         string signalfile, string partrecofile, string combinatorialfile, string JpsiLeakfile,
         double minBMass = 4880, double maxBMass = 5700,
         string signaltree = "DecayTree", string partrecotree = "DecayTree", string combinatorialtree = "DecayTree", string JpsiLeaktree = "DecayTree");


   void generate(bool wantPlots = false, string plotsfile = "toErase.root");
   void fit(bool wantplot, bool constPartReco,
         double fracPartReco_const, ofstream& out, TTree* t, bool update, string plotsfile);


   protected:

   void initiateParams(RooArgSet* parset);

   void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, RooRealVar& nSignal, RooRealVar& nPartReco,
         RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma,
         RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee,
                       RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen , bool constFracs, bool constComb);


   void plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataSet dataGenTot);
   void plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataHist dataGenTot);

   void PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT);
   void PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT);
   void PlotShape1D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M);


   int nGenSignal;
   int nGenPartReco;
   int nGenComb;
   int nGenJpsiLeak;
   double nGenFracZeroGamma;  
   double nGenFracOneGamma;
   string workspacename;


};


#endif
