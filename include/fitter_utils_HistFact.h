#ifndef FITTER_UTILS_HISTFACT_H
#define FITTER_UTILS_HISTFACT_H

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
#include "TH2D.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"

using namespace std;
using namespace RooFit;
using namespace RooStats ;
using namespace HistFactory ;

class FitterUtilsHistFact
{

   public:

  FitterUtilsHistFact(int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_, double nGenFracOneGamma_, string workspacename_);

  void prepare_PDFs(string trigStr, string weightStr, string BDTVar, double BDTcut,
         string signalfile, string partrecofile, string combinatorialfile, string JpsiLeakfile,
         double minBMass = 4880, double maxBMass = 5700,
         string signaltree = "DecayTree", string partrecotree = "DecayTree", string combinatorialtree = "DecayTree", string JpsiLeaktree = "DecayTree");


   void generate(bool wantPlots = false, string plotsfile = "toErase.root");
   void fit(bool wantplot, bool constPartReco,
         double fracPartReco_const, ofstream& out, TTree* t, bool update, string plotsfile);


  TH2D* make_base_histogram(string componentName, string componentFile, string componentTree,
                      // string cut, string weight, RooArgSet varset, RooRealVar &x, RooRealVar &y, string binningName)
                      string cut, vector<string> varset, RooRealVar &x, RooRealVar &y, 
                      string binningName);
  
  void fit_combinatorial(string combfile, string combtree, string combcuts, vector<string> argsetcomb, 
                         RooRealVar &B_plus_M, RooRealVar &misPT, RooWorkspace *workspace);
  

   protected:

   void initiateParams(RooArgSet* parset);

   void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, RooRealVar& nSignal, RooRealVar& nPartReco,
         RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma,
         RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee,
                       RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen , bool constFracs, bool constComb);


  void plot_fit_result(string plotsfile, RooSimultaneous &totPdf, RooDataHist *dataGenTot,
                       double nsig, double npartreco, double ncomb, double njpsi);

   void PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT);
   void PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT);
   void PlotShape1D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M);



  double binWidthMass;
  double binWidthMassBroader;
  double binWidthMassExtended;
  int nMassBins;
  int nmisPTBins;
  int nGenSignal;
  int nGenPartReco;
  int nGenComb;
  int nGenJpsiLeak;
  double nGenFracZeroGamma;  
  double nGenFracOneGamma;
  string workspacename;
  


};


#endif
