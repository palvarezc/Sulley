#ifndef FITTER_UTILS_ANA_H
#define FITTER_UTILS_ANA_H

#include "Utils.h"
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
#include "RooBinning.h"
#include "Combinatorial.h"
#include "RooExtendPdf.h"
#include "RooSimultaneous.h"
#include "RooCBShape.h"
#include "RooAbsPdf.h"

using namespace std;
using namespace RooFit;
using namespace Sulley;

class FitterUtilsAna
{

public:

  FitterUtilsAna( string name, const Options &opts_);
  

  void prepare_PDFs(int fitmode);  
  void prepare_PDFs_ee();  
  void prepare_PDFs_mumu();  

  void generate(int fitmode);
  void generate_ee();
  void generate_mumu();

  void run_toy(double fracPartReco_const, 
	       ofstream& out, 
	       TTree* t, 
	       bool update,
	       int fitmode);

  void fit_ee(RooDataSet* dataset_ee, 
	      RooDataSet* dataset_Kemu, 
	      RooWorkspace* workspace,
	      double fracPartReco_const, 
	      ofstream& out, 
	      TTree* t, 
	      bool update);

  void fit_mumu(RooDataSet* datasetTot, 
		RooWorkspace* workspace,
		double fracPartReco_const, 
		ofstream& out, 
		TTree* t, 
		bool update);

  void fit_RK(RooDataSet* dataset_ee, 
	      RooDataSet* dataset_mumu, 
	      RooDataSet* dataset_Kemu, 
	      RooWorkspace* workspace_ee,
	      RooWorkspace* workspace_mumu,
	      double fracPartReco_const, 
	      ofstream& out, 
	      TTree* t, 
	      bool update);

  void enablePlotting() 
  {
    m_wantplot = true;
  }

  void disablePlotting() 
  {
    m_wantplot = false;
  }
  



  void initiateParams(RooArgSet parset, RooArgSet constraints = RooArgSet(), RooArgSet parset_const = RooArgSet());

   void initiateParams(RooRealVar const& expoConstGen,
                       RooRealVar& nSignal, RooRealVar& nPartReco, 
                       RooRealVar& nComb, RooRealVar& fracZero, 
                       RooRealVar& fracOne, RooRealVar& expoConst, 
                       RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma);



  double prepare_component_PDF_ee(string filename, string treename,
				  string compname, string cuts, 
				  RooArgSet varset, RooWorkspace *workspace,
				  string binningname);

  double prepare_component_PDF_mumu(string filename, string treename,
				  string compname, string cuts, 
				  RooArgSet varset, RooWorkspace *workspace,
				  string binningname);
    


   void plot_fit_result(RooAbsPdf &totPdf, RooDataSet dataGenTot);
   void plot_mumu_fit_result(RooAbsPdf &totPdf, RooDataSet dataGenTot);
  void plot_kemu_fit_result(RooAbsPdf &totKemuPdf, RooDataSet datasetKemu);
  
   void PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, 
                  string canvName, RooArgSet varset);
   void PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, 
                    string canvName, RooRealVar& B_plus_M, RooRealVar& misPT);
   void PlotShape1D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, 
                    string canvName, RooRealVar& B_plus_M);


  // void fit_combinatorial(string combcuts, RooArgSet *varset, RooWorkspace *workspace);
  
  void generate_component(string compname, RooArgSet varset, int nGenComp,
                          RooWorkspace* workspace, RooWorkspace* workspaceGen);

  void generate_component_ana(string compname, RooArgSet varset, int nGenComp,
                              RooWorkspace* workspaceMCFit, RooWorkspace* workspaceGen,
                              int trigCat, int PhotCat);
  
  RooAbsPdf* build_component_ana(string compname, RooArgSet varset, 
                                 RooWorkspace* workspaceMCFit, 
                                 int trigCat, int PhotCat,
                                 RooRealVar* meanShift, RooRealVar* sigmaScaleFactor);
  

  string name;
  Combinatorial combinatorial;
  Combinatorial Kemu;
  string defaultbinningname;
  Options opts;
  bool m_wantplot;

};


#endif
