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


using namespace RooFit;
using namespace std;

int main()
{



   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV/c^{2}");
   B_plus_DTFM_M_zero.setBins(100);
   RooRealVar BDTKeeBig3("BDTKeeBig3", "BDTKeeBig3", -1,1);
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", 4500, 6200, "MeV/c^{2}");
   RooRealVar B_plus_M_corr("B_plus_M_corr", "M_{cor}", 4600, 7000, "MeV/c^{2}");
   B_plus_M_corr.setBins(15);

   RooRealVar T("T", "T", 62.6567, 0, 200);
   RooRealVar n("n", "n", 4.36833, 1., 5.5);
   RooRealVar expoConst("expoConst", "expoConst", 0.000187787, -1, 1); 
//   RooRealVar T("T", "T", 97, 0, 200);
//   RooRealVar n("n", "n", 3.5, 1., 5.5);
//   RooRealVar expoConst("expoConst", "expoConst", -1e-3, -1, 1); 
   n.setConstant(true);

   //************2D Tsallis distributionn for combinatorial

   RooMcorMvisTsallis McorMvis("McorMvis", "McorMvis", B_plus_M_corr, B_plus_M, T, n, expoConst);

//   for(double cor(4600); cor<7000; ++cor)
//   {
//      for(double vis(4500); vis<6200; ++vis)
//      {
//         B_plus_M_corr.setVal(cor);
//         B_plus_M.setVal(vis);
//         cout<<McorMvis.evaluate()<<endl;
//      }
//   }




   RooArgSet argset2(B_plus_M, B_plus_M_corr);

   RooPlot* plotM = B_plus_M.frame();
   McorMvis.plotOn(plotM);
   
   TCanvas canv("canv", "canv", 800, 800);
   plotM->Draw();
   canv.Print("cacavis.pdf");

   RooPlot* plotMCorr = B_plus_M_corr.frame();
   McorMvis.plotOn(plotMCorr);
   
   TCanvas canv2("canv2", "canv2", 800, 800);
   plotMCorr->Draw();
   canv2.Print("cacacor.pdf");
   
   TCanvas canv3("canv3", "canv3", 800, 800);
   TH2F* th2fKey = (TH2F*)McorMvis.createHistogram("th2Shape", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
   th2fKey->Draw("surf");
   canv3.Print("caca2D.root");

   RooAbsPdf::GenSpec* GenSpecComb = McorMvis.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(20000));

   RooDataSet* dataGenComb = McorMvis.generate(*GenSpecComb);
 
   TH2F* th2fGenComb = (TH2F*)dataGenComb->createHistogram("th2fGenComb", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
   TH1D* th1dGenComb = th2fGenComb->ProjectionX();

   TCanvas cCaca("cCaca", "cCaca", 600, 600);
   th1dGenComb->Draw();
   cCaca.Print("caca1.pdf");

   RooPlot* plotM2 = B_plus_M.frame();
   dataGenComb->plotOn(plotM2);
   McorMvis.plotOn(plotM2);

   TCanvas cCaca2("cCaca2", "cCaca2", 600, 600);
   plotM2->Draw();
   cCaca2.Print("caca2.pdf");
   


   
   


//   expoConst.setVal(+1e-3);
//
//   dataGenComb = McorMvis.generate(*GenSpecComb);
//   
//   th2fGenComb = (TH2F*)dataGenComb->createHistogram("th2fGenComb", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
//   th1dGenComb = th2fGenComb->ProjectionX();
//
//   cCaca.cd();
//   th1dGenComb->Draw();
//   cCaca.Print("caca2.pdf");
//
//
//
//   //********************
//   expoConst.setVal(+1e-3);
//
//   RooAbsPdf::GenSpec* GenSpecComb2 = McorMvis.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(20000));
//   dataGenComb = McorMvis.generate(*GenSpecComb2);
//   
//   th2fGenComb = (TH2F*)dataGenComb->createHistogram("th2fGenComb", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
//   th1dGenComb = th2fGenComb->ProjectionX();
//
//   cCaca.cd();
//   th1dGenComb->Draw();
//   cCaca.Print("caca3.pdf");


   return 0;
}

