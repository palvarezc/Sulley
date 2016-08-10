#ifndef FITTER_UTILS_SIMULTANEOUS_EXPOFPOLYTIMESX_H
#define FITTER_UTILS_SIMULTANEOUS_EXPOFPOLYTIMESX_H

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
#include<RooSimultaneous.h>
#include<RooCategory.h>
#include "RooMcorMvisTsallis.h"
#include "RooPTMVis.h"
#include "usefulFunctions.h"
#include "RooNumIntConfig.h"
#include "RooExtendPdf.h"

#include "fitter_utils_ExpOfPolyTimesX.h"

using namespace std;
using namespace RooFit;

class FitterUtilsSimultaneousExpOfPolyTimesX: public FitterUtilsExpOfPolyTimesX
{

   public:

      FitterUtilsSimultaneousExpOfPolyTimesX(int nGenKemu_, int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_,
            double nGenFracZeroGamma_, double nGenFracOneGamma_, string workspacename_);


      void generate(bool wantplot = false, string plotsfile = "toErase.root");
      void fit(bool wantplot, bool constPartReco,
            double fracPartReco_const, ofstream& out, TTree* t, bool update, string plotsfile);


   protected:

      void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma,
            RooRealVar& nKemu, RooRealVar& nSignal, RooRealVar& nPartReco,
            RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne,
            RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma,
            RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee,
            RooRealVar& l1Kemu, RooRealVar& l2Kemu, RooRealVar& l3Kemu, RooRealVar& l4Kemu, RooRealVar& l5Kemu,
            RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen);

         void plot_kemu_fit_result(string plotsfile, RooAbsPdf &totKemuPdf, RooDataSet const& dataGenKemu);

      int nGenKemu;
};


#endif
