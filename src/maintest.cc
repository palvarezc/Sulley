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
#include<TH1D.h>
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
#include "fitter_utils_1Dsimultaneous.h"


using namespace RooFit;
using namespace std;

int main()
{

   string fSignal("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/B2Kee_Strip21_MC_trigged.root");
   string fComb("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/B2Kemu_Strip21_trigged.root");
   string fJpsiLeak("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/total_signal_leakage_trigged.root");
   string fPartReco("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/BJpsiX_Strip21_MC2012_ctrl_noDTFMAllmass_trigged_rarebkgs.root");


   FitterUtils1DSimultaneous fuS(180,93, 500, 8, 600, 0.3, 0.4, 600, "workspacecaca.root");

   // fuS.prepare_PDFs("L0ETOSOnly_d", "BDTNewu4bR", 0.14,
   //       fSignal, fPartReco, fComb, fJpsiLeak, 4880, 6200);

   // ADDED WEIGHTSTR BUT NOT USED!
   fuS.prepare_PDFs("L0ETOSOnly_d", "L0ETOSOnly_d", "BDTNewu4bR", 0.14,
         fSignal, fPartReco, fComb, fJpsiLeak, 4880, 6200);

   fuS.display();


   bool wantplot(true);
   bool update(false);

   TFile f("resultscaca.root","update");
   TTree tKemu("tKemu", "tKemu");
   TTree tKee("tKee", "tKee");
   ofstream out("outcaca.dat");
   for(int i(0); i<10; ++i)
   {
      fuS.generate();
      fuS.fit(wantplot, false, 0.2, out, &tKee, &tKemu, update, "plotscaca.root");

      wantplot = false;
      update = true;
   }

   f.cd();
   tKemu.Write();
   tKee.Write();
   f.Close();
   return 0;
}

