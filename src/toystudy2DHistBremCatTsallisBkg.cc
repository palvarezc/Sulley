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



int main(int argc, char* argv[])
{

   //***********Number of events to generate
  
  
  int nGenSignal(203);
  int nGenPartReco(79);
  int nGenComb(31);

   //*********** Get arguments and set stuff


   string trigStr("L0ETOSOnly_d");
   bool wantOldDataSet(false);
   string BDTcut("0.257");
   int ntoys(1000);

   if(argc == 8)
   {
      if(*argv[1] == '1') wantOldDataSet = true;

      nGenSignal = atoi(argv[2]);
      nGenPartReco = atoi(argv[3]);
      nGenComb = atoi(argv[4]);

      int trigInt(atoi(argv[5]));
      if(trigInt == 1) trigStr = "L0HTOSOnly_d";
      if(trigInt == 2) trigStr = "L0TISOnly_d";
      BDTcut = argv[7];
      ntoys = atoi(argv[6]);
   }

   int nGenSignalZeroGamma(floor(0.368*nGenSignal));
   int nGenSignalOneGamma(floor(0.484*nGenSignal));
   int nGenSignalTwoGamma(floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));


   //***********Get the datasets
   TFile* fSignal = new TFile("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldtrees/B2Kee_Strip21_BDT_ctrl_trigged.root");
   TTree* tSignal = (TTree*)fSignal->Get("DecayTree");
   TFile* fPartReco = new TFile("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldtrees/B2Kee_Strip21_BDT_prc_trigged.root");
   TTree* tPartReco = (TTree*)fPartReco->Get("DecayTree");

   TFile* fComb = new TFile("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldtrees/B2Kee_Strip21_piee_trigged.root");
   TTree* tComb = (TTree*)fComb->Get("DecayTree"); 


   RooRealVar B_plus_M_corr("B_plus_M_corr", "M_{cor}", 5000, 10000, "MeV/c^{2}");
   B_plus_M_corr.setBins(15);
   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV/c^{2}"); 
   B_plus_DTFM_M_zero.setBins(100); 
   RooRealVar BDTKeeBig("BDTKeeBig3", "BDTKeeBig3", -1,1);
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", 4880, 5700, "MeV/c^{2}");
   // RooRealVar L0ETOSOnly_d("L0ETOSOnly_d", "L0ETOSOnly_d", -10, 10);
   RooRealVar trigVar(trigStr.c_str(), trigStr.c_str(), -10, 10);
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);
   B_plus_M.setBins(20);


   // RooArgSet argset(BDTKeeBig, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, L0ETOSOnly_d, e_plus_BremMultiplicity, e_minus_BremMultiplicity);
   RooArgSet argset(BDTKeeBig, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity);

   cout<<"getting the datasets:"<<endl;

   RooDataSet* dataSetSignalZeroGamma;
   RooDataSet* dataSetSignalOneGamma;
   RooDataSet* dataSetSignalTwoGamma;
   RooDataSet* dataSetPartReco;
   RooDataSet* dataSetComb;

   TFile* fw;

   if(!wantOldDataSet)
   {
      dataSetSignalZeroGamma = new RooDataSet("dataSetSignalZeroGamma", "dataSetSignalZeroGamma", tSignal, argset,( " ("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5)").c_str());
      dataSetSignalOneGamma = new RooDataSet("dataSetSignalOneGamma", "dataSetSignalOneGamma", tSignal, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5)").c_str());
      dataSetSignalTwoGamma = new RooDataSet("dataSetSignalTwoGamma", "dataSetSignalTwoGamma", tSignal, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5)").c_str());
      dataSetPartReco = new RooDataSet("dataSetPartReco", "dataSetPartReco", tPartReco, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+")").c_str());
      dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("("+trigStr+"  > 0.9)").c_str());

   }

   if(wantOldDataSet)
   {
      fw = new TFile("toystudy2DHistBremCatTsallisBkg_L0ETOS.root");
      RooWorkspace* workspace = (RooWorkspace*)fw->Get("workspace");
      dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
      dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
      dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
      dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
      dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
   }


   cout<<"binning the datasets:"<<endl;

   RooArgSet argset2(B_plus_M, B_plus_M_corr);
   RooDataHist dataHistSignalZeroGamma("dataHistSignalZeroGamma", "dataHistSignalZeroGamma", argset2, *dataSetSignalZeroGamma); 
   RooDataHist dataHistSignalOneGamma("dataHistSignalOneGamma", "dataHistSignalOneGamma", argset2, *dataSetSignalOneGamma); 
   RooDataHist dataHistSignalTwoGamma("dataHistSignalTwoGamma", "dataHistSignalTwoGamma", argset2, *dataSetSignalTwoGamma); 
   RooDataHist dataHistComb("dataHistComb", "dataHistComb", argset2, *dataSetComb); 
   RooDataHist dataHistPartReco("dataHistPartReco", "dataHistPartReco", argset2, *dataSetPartReco); 

   //***************Create 2D histogram estimates from data


   cout<<"Preparing the 3 1D histPdf: 1"<<endl;
   //   RooArgSet argset2(B_plus_M);
   RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco,2); cout<<" 3"<<endl;

   RooRealVar T("T", "T", 97, 0, 200);
   RooRealVar n("n", "n", 3.5, 1., 5.5);
   RooRealVar T_fit("T_fit", "T_fit", 97, 0, 200);
   RooRealVar n_fit("n_fit", "n_fit", 3.5, 1., 5.5);
   RooRealVar expoConst("expoConst", "expoConst", -1e-3, -1, 1);
   //   n.setConstant(true);



   //************2D Tsallis distributionn for combinatorial

   RooMcorMvisTsallis McorMvis("McorMvis", "McorMvis", B_plus_M_corr, B_plus_M, T, n, expoConst);
   McorMvis.fitTo(*dataSetComb); // 

   T.setConstant(true);
   n.setConstant(true);

   std::cout<<"T generated is: "<<T.getVal()<<std::endl;


   RooMcorMvisTsallis McorMvis_fit("McorMvis_fit", "McorMvis_fit", B_plus_M_corr, B_plus_M, T_fit, n_fit, expoConst);

   // T_fit.setVal(50);
   T_fit.setVal(T.getVal());
   n_fit.setVal(n.getVal());
   T_fit.setConstant(true);
   n_fit.setConstant(true);

   std::cout<<"T fitted is: "<<T_fit.getVal()<<std::endl;



   double trueExpConst(expoConst.getValV());

   cout<<"Preparing the generation of events 1";

   RooAbsPdf::GenSpec* GenSpecSignalZeroGamma = histPdfSignalZeroGamma.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalZeroGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalOneGamma = histPdfSignalOneGamma.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalOneGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalTwoGamma = histPdfSignalTwoGamma.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalTwoGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecPartReco =  histPdfPartReco.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenPartReco)); cout<<" 3 "<<endl;
   RooAbsPdf::GenSpec* GenSpecComb = McorMvis.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenComb));


   //***************Save everything on a workspace

   if(!wantOldDataSet)
   {
      RooWorkspace workspace("workspace", "workspace");
      workspace.import(B_plus_DTFM_M_zero);
      //   workspace.import(B_plus_M);
      workspace.import(*dataSetSignalZeroGamma);
      workspace.import(*dataSetSignalOneGamma);
      workspace.import(*dataSetSignalTwoGamma);
      workspace.import(*dataSetPartReco);
      workspace.import(*dataSetComb);
      workspace.import(expoConst);
      workspace.writeToFile("toystudy2DHistBremCatTsallisBkg_L0ETOS.root");
   }

   //***************Prepare the stuff to generate events


   TFile f("toystudy2DHistBremCatTsallisBkg_results.root", "RECREATE");
   TTree t("paramsFloatingExpFloatingFracPartReco_L0ETOS", "paramsFloatingExpFloatingFracPartReco_L0ETOS");
   TH1F hSignal("hErrSignal", "hErrSignal", 100, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   TH1F hPartReco("hErrPartReco", "hErrPartReco", 100, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   TH1F hComb("hErrComb", "hErrComb", 100, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));

   TH1F hExpConst("hExpConst", "hExpConst", 100, -10, 10);
   TH1F hExpConstPull("hPullExpConst", "hPullExpConst", 100, -10, 10);


   TH1F hMinuitErr("hMinuitErr", "hMinuitErr", 100, 0, 40);
   TH1F hPull("hPull", "hPull", 100, -5, 5);

   bool wantPlots(true);
   bool update(false);

   ofstream out("fitResult2D.dat");

   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Generation and fit number "<<i<<endl;
      toystudy1DMakeFit("toystudy2DHistBremCatTsallisBkg_L0ETOS.root", wantPlots, histPdfSignalZeroGamma, histPdfSignalOneGamma, histPdfSignalTwoGamma, histPdfPartReco, McorMvis, McorMvis_fit, *GenSpecSignalZeroGamma,*GenSpecSignalOneGamma, *GenSpecSignalTwoGamma,  *GenSpecPartReco, *GenSpecComb, nGenSignalZeroGamma,nGenSignalOneGamma, nGenSignalTwoGamma, nGenPartReco, nGenComb, trueExpConst, expoConst,  out, &t,update);

      wantPlots = false;
      update = true;
   }

   out.close();

   f.cd();
   t.Write();

   ofstream out2("resultToy.dat", std::ios::app);
   out2<<"2D Mvis X MCorr with Tsallis bkg:"<<endl;
   makeTableResults(&t, nGenSignal, nGenPartReco, nGenComb, out2 );
   out2.close();


   fComb->Close();
   fSignal->Close();
   fPartReco->Close();
   f.Close();
   fw->Close();

   delete fComb;
   delete fSignal;
   delete fPartReco;
   delete fw;

   return 0;
}


void toystudy1DMakeFit(string workspacename,  bool wantplot, 
      RooHistPdf& histPdfSignalZeroGamma, RooHistPdf& histPdfSignalOneGamma, RooHistPdf& histPdfSignalTwoGamma, 
                       RooHistPdf& histPdfPartReco, RooMcorMvisTsallis& McorMvis, RooMcorMvisTsallis& McorMvis_fit, 
      RooAbsPdf::GenSpec& GenSpecSignalZeroGamma, RooAbsPdf::GenSpec& GenSpecSignalOneGamma, RooAbsPdf::GenSpec& GenSpecSignalTwoGamma, 
      RooAbsPdf::GenSpec& GenSpecPartReco, RooAbsPdf::GenSpec& GenSpecComb, int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, double expoConstGen, RooRealVar& expoConst, ofstream& out,  TTree* t, bool update)
{
   //***************Get the stuff from the workspace

  TFile f(workspacename.c_str());
  RooWorkspace* workspace = (RooWorkspace*)f.Get("workspace"); 

  RooRealVar B_plus_M(*workspace->var("B_plus_M")); 
  RooRealVar B_plus_M_corr(*workspace->var("B_plus_M_corr")); 
  RooRealVar B_plus_DTFM_M_zero(*workspace->var("B_plus_DTFM_M_zero")); 


   cout<<"Variable loaded:"<<endl;
   B_plus_M.Print(); B_plus_M_corr.Print(); B_plus_DTFM_M_zero.Print(); expoConst.Print();


   //***************Generate some datasets

   RooArgSet argset(B_plus_M);
   cout<<"Generating signal Zero Photon"<<endl;
   RooDataSet* dataGenSignalZeroGamma = histPdfSignalZeroGamma.generate(GenSpecSignalZeroGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating signal One Photon"<<endl;
   RooDataSet* dataGenSignalOneGamma = histPdfSignalOneGamma.generate(GenSpecSignalOneGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating signal two Photons"<<endl;
   RooDataSet* dataGenSignalTwoGamma = histPdfSignalTwoGamma.generate(GenSpecSignalTwoGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating combinatorial"<<endl;
   RooDataSet* dataGenComb = McorMvis.generate(GenSpecComb);//(argset, 100, false, true, "", false, true);
   cout<<"Generating PartReco"<<endl;
   RooDataSet* dataGenPartReco = histPdfPartReco.generate(GenSpecPartReco);//argset, 160, false, true, "", false, true);

 
   if(wantplot)
   {
      //**************Must get the datasets

      RooDataSet* dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
      RooDataSet* dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
      RooDataSet* dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
      RooDataSet* dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
      RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");

      //**************Prepare TFile to save the plots

      TFile f2("plots2DHistBremCatTsallisBkgfit.root", "UPDATE");

      //**************Plot Signal Zero Gamma

      TH2F* th2fKeySignalZeroGamma = (TH2F*)histPdfSignalZeroGamma.createHistogram("th2fKeySignalZeroGamma", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
      TH2F* th2fGenSignalZeroGamma = (TH2F*)dataGenSignalZeroGamma->createHistogram("th2fGenSignalZeroGamma", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));

      RooPlot* plotMSignalZeroGamma = B_plus_M.frame();
      dataSetSignalZeroGamma->plotOn(plotMSignalZeroGamma);
      histPdfSignalZeroGamma.plotOn(plotMSignalZeroGamma);

      RooPlot* plotMCorrSignalZeroGamma = B_plus_M_corr.frame();
      dataSetSignalZeroGamma->plotOn(plotMCorrSignalZeroGamma);
      histPdfSignalZeroGamma.plotOn(plotMCorrSignalZeroGamma);

      TCanvas* cSignalZeroGamma = new TCanvas("cSignalZeroGamma", "cSignalZeroGamma", 800, 800);
      cSignalZeroGamma->Divide(2,2);
      cSignalZeroGamma->cd(1); th2fGenSignalZeroGamma->Draw("lego");
      cSignalZeroGamma->cd(2); th2fKeySignalZeroGamma->Draw("surf");
      cSignalZeroGamma->cd(3); plotMSignalZeroGamma->Draw();
      cSignalZeroGamma->cd(4); plotMCorrSignalZeroGamma->Draw();

      cSignalZeroGamma->Write();

      //**************Plot Signal one Gamma

      TH2F* th2fKeySignalOneGamma = (TH2F*)histPdfSignalOneGamma.createHistogram("th2fKeySignalOneGamma", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
      TH2F* th2fGenSignalOneGamma = (TH2F*)dataGenSignalOneGamma->createHistogram("th2fGenSignalOneGamma", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));

      RooPlot* plotMSignalOneGamma = B_plus_M.frame();
      dataSetSignalOneGamma->plotOn(plotMSignalOneGamma);
      histPdfSignalOneGamma.plotOn(plotMSignalOneGamma);

      RooPlot* plotMCorrSignalOneGamma = B_plus_M_corr.frame();
      dataSetSignalOneGamma->plotOn(plotMCorrSignalOneGamma);
      histPdfSignalOneGamma.plotOn(plotMCorrSignalOneGamma);

      TCanvas* cSignalOneGamma = new TCanvas("cSignalOneGamma", "cSignalOneGamma", 800, 800);
      cSignalOneGamma->Divide(2,2);
      cSignalOneGamma->cd(1); th2fGenSignalOneGamma->Draw("lego");
      cSignalOneGamma->cd(2); th2fKeySignalOneGamma->Draw("surf");
      cSignalOneGamma->cd(3); plotMSignalOneGamma->Draw();
      cSignalOneGamma->cd(4); plotMCorrSignalOneGamma->Draw();

      cSignalOneGamma->Write();

      //**************Plot Signal two Gamma

      TH2F* th2fKeySignalTwoGamma = (TH2F*)histPdfSignalTwoGamma.createHistogram("th2fKeySignalTwoGamma", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
      TH2F* th2fGenSignalTwoGamma = (TH2F*)dataGenSignalTwoGamma->createHistogram("th2fGenSignalTwoGamma", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));

      RooPlot* plotMSignalTwoGamma = B_plus_M.frame();
      dataSetSignalTwoGamma->plotOn(plotMSignalTwoGamma);
      histPdfSignalTwoGamma.plotOn(plotMSignalTwoGamma);

      RooPlot* plotMCorrSignalTwoGamma = B_plus_M_corr.frame();
      dataSetSignalTwoGamma->plotOn(plotMCorrSignalTwoGamma);
      histPdfSignalTwoGamma.plotOn(plotMCorrSignalTwoGamma);

      TCanvas* cSignalTwoGamma = new TCanvas("cSignalTwoGamma", "cSignalTwoGamma", 800, 800);
      cSignalTwoGamma->Divide(2,2);
      cSignalTwoGamma->cd(1); th2fGenSignalTwoGamma->Draw("lego");
      cSignalTwoGamma->cd(2); th2fKeySignalTwoGamma->Draw("surf");
      cSignalTwoGamma->cd(3); plotMSignalTwoGamma->Draw();
      cSignalTwoGamma->cd(4); plotMCorrSignalTwoGamma->Draw();

      cSignalTwoGamma->Write();

      //**************Plot Part Reco 

      TH2F* th2fKeyPartReco = (TH2F*)histPdfPartReco.createHistogram("th2fKeyPartReco", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
      TH2F* th2fGenPartReco = (TH2F*)dataGenPartReco->createHistogram("th2fGenPartReco", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));

      RooPlot* plotMPartReco = B_plus_M.frame();
      dataSetPartReco->plotOn(plotMPartReco);
      histPdfPartReco.plotOn(plotMPartReco);

      RooPlot* plotMCorrPartReco = B_plus_M_corr.frame();
      dataSetPartReco->plotOn(plotMCorrPartReco);
      histPdfPartReco.plotOn(plotMCorrPartReco);

      TCanvas* cPartReco = new TCanvas("cPartReco", "cPartReco", 800, 800);
      cPartReco->Divide(2,2);
      cPartReco->cd(1); th2fGenPartReco->Draw("lego");
      cPartReco->cd(2); th2fKeyPartReco->Draw("surf");
      cPartReco->cd(3); plotMPartReco->Draw();
      cPartReco->cd(4); plotMCorrPartReco->Draw();

      cPartReco->Write();

      //**************Plot Combinatorial 

      TH2F* th2fKeyComb = (TH2F*)McorMvis.createHistogram("th2fKeyComb", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
      TH2F* th2fGenComb = (TH2F*)dataGenComb->createHistogram("th2fGenComb", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));

      RooPlot* plotMComb = B_plus_M.frame();
      dataSetComb->plotOn(plotMComb);
      McorMvis.plotOn(plotMComb);

      RooPlot* plotMCorrComb = B_plus_M_corr.frame();
      dataSetComb->plotOn(plotMCorrComb);
      McorMvis.plotOn(plotMCorrComb);

      TCanvas* cComb = new TCanvas("cComb", "cComb", 800, 800);
      cComb->Divide(2,2);
      cComb->cd(1); th2fGenComb->Draw("lego");
      cComb->cd(2); th2fKeyComb->Draw("surf");
      cComb->cd(3); plotMComb->Draw();
      cComb->cd(4); plotMCorrComb->Draw();

      cComb->Write();


   }

   f.cd();

   //***************Merge datasets

   RooDataSet* dataGenTot(dataGenPartReco);
   dataGenTot->append(*dataGenSignalZeroGamma);
   dataGenTot->append(*dataGenSignalOneGamma);
   dataGenTot->append(*dataGenSignalTwoGamma);
   dataGenTot->append(*dataGenComb);

   RooWorkspace workspaceb("workspaceb", "workspaceb");
   dataGenTot->SetNameTitle("myexperiment","myexperiment");
   workspaceb.import(*dataGenTot);
   workspaceb.writeToFile("my_experiment.root");
   
   

   //**************Prepare fitting function

   double nGenSignal(nGenSignalZeroGamma+nGenSignalOneGamma+nGenSignalTwoGamma);

      RooRealVar fracPartReco("fracPartReco", "fracPartReco", nGenPartReco/(1.*nGenSignal), 0, 2);
   RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   //RooFormulaVar nPartReco("nPartReco", "fracPartReco*nSignal", RooArgList(fracPartReco, nSignal)); 
   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);


   RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));

   RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", RooArgList(histPdfSignalZeroGamma, histPdfSignalOneGamma, histPdfSignalTwoGamma), RooArgList(fracZero, fracOneRec));
   RooAddPdf totPdf("totPdf", "totPdf", RooArgList(histPdfSignal, histPdfPartReco, McorMvis_fit), RooArgList(nSignal, nPartReco, nComb));

   //**************** Constrain the fraction of zero and one photon

   RooRealVar fracZeroConstMean("fracZeroConstMean", "fracZeroConstMean", nGenSignalZeroGamma/nGenSignal);
   RooRealVar fracZeroConstSigma("fracZeroConstSigma", "fracZeroConstSigma", sqrt(nGenSignalZeroGamma)/nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroConstMean, fracZeroConstSigma); 

   RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", nGenSignalOneGamma/nGenSignal);
   RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", sqrt(nGenSignalOneGamma)/nGenSignal);
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneConstMean, fracOneConstSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", nGenPartReco/(1.*nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", (1./(1.*nGenSignal))*sqrt(nGenPartReco + (nGenPartReco*nGenPartReco/(1.*nGenSignal))) ) ;
    RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);


   //**************** fit   

   initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, nGenPartReco, nGenComb, expoConstGen, nSignal, nPartReco, nComb, fracZero, fracOne, expoConst);

   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(RooArgSet(fracZeroConst, fracOneConst, fracPartRecoConst)));
   RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(RooArgSet(fracZeroConst, fracOneConst)));
   RooMinuit minuit(*nll);
   minuit.setStrategy(2);

   int migradRes(1);
   int hesseRes(4);
   double edm(10);
   int nrefit(0);

   RooFitResult* fitRes(0);

   for(int i(0); (i<10) && ( (migradRes != 0) || (hesseRes !=0) || (edm > 1e-4)); ++i)
   {
      initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, nGenPartReco, nGenComb, expoConstGen, nSignal, nPartReco, nComb, fracZero, fracOne, expoConst);
      cout<<"FITTING: starting with nsignal = "<<nSignal.getValV()<<" refit nbr. "<<i<<endl;
      if(fitRes != NULL && fitRes != 0) delete fitRes;

      migradRes = minuit.migrad();
      hesseRes = minuit.hesse();

      fitRes = minuit.save();
      edm = fitRes->edm();

      ++nrefit;
   }

   fillTreeResult(t, fitRes,  update);
   delete fitRes;
   //totPdf.fitTo(*dataGenTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(wantplot)
   {
      //**************Prepare TFile to save the plots

      TFile f2("plots2DHistBremCatTsallisBkgfit.root", "UPDATE");
      //**************Plot the results of the fit

      RooPlot* frameVis = B_plus_M.frame();
      dataGenTot->plotOn(frameVis);
      totPdf.plotOn(frameVis, Components("histPdfPartReco"), LineColor(kBlue));
      totPdf.plotOn(frameVis, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
      totPdf.plotOn(frameVis, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
      totPdf.plotOn(frameVis, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
      totPdf.plotOn(frameVis, Components("McorMvis"), LineColor(kBlack));
      totPdf.plotOn(frameVis, LineColor(kRed));


      RooPlot* frameCorr = B_plus_M_corr.frame();
      dataGenTot->plotOn(frameCorr);
      totPdf.plotOn(frameCorr, Components("histPdfPartReco"), LineColor(kBlue));
      totPdf.plotOn(frameCorr, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
      totPdf.plotOn(frameCorr, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
      totPdf.plotOn(frameCorr, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
      totPdf.plotOn(frameCorr, Components("McorMvis_fit"), LineColor(kBlack));
      totPdf.plotOn(frameCorr, LineColor(kRed));


      TCanvas* cFit = new TCanvas("cFit", "cFit", 600, 800);
      cFit->Divide(1,2);
      cFit->cd(1); frameCorr->Draw();
      cFit->cd(2); frameVis->Draw();

      cFit->Write();
      f2.Close();
   }

   f.cd();
   //delete and return

   delete dataGenSignalZeroGamma; 
   delete dataGenSignalOneGamma; 
   delete dataGenSignalTwoGamma; 
   delete dataGenPartReco; 
   delete dataGenComb; 
   delete workspace;
   delete nll;
   f.Close();
}

void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, double expoConstGen, RooRealVar& nSignal,
      RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst)
{
   TRandom rand;
   rand.SetSeed();

   int nGenSignal = nGenSignalZeroGamma + nGenSignalOneGamma + nGenSignalTwoGamma;
   double nGenSignal2 = rand.Uniform(nGenSignal-5*sqrt(nGenSignal), nGenSignal+5*sqrt(nGenSignal));
   double nGenPartReco2 = rand.Uniform(nGenPartReco-5*sqrt(nGenPartReco), nGenPartReco+5*sqrt(nGenPartReco));
   double nGenComb2 = rand.Uniform(nGenComb-5*sqrt(nGenComb), nGenComb+5*sqrt(nGenComb));


   nSignal.setVal(nGenSignal2);
   nSignal.setRange(TMath::Max(0.,nGenSignal2-10.*sqrt(nGenSignal)) , nGenSignal2+10*sqrt(nGenSignal));

   //fracPartReco.setVal(nGenPartReco2/nGenSignal2);
   //fracPartReco.setRange(0,2);

   nPartReco.setVal(nGenPartReco2);
   nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(nGenPartReco)), nGenPartReco2+10*sqrt(nGenPartReco));

   nComb.setVal(nGenComb2);
   nComb.setRange(TMath::Max(0.,nGenComb2-10.*sqrt(nGenComb)), nGenComb2+10*sqrt(nGenComb));

   double fracGenZero(nGenSignalZeroGamma/(1.*nGenSignal));
   double fracGenOne(nGenSignalOneGamma/(1.*nGenSignal));


   fracZero.setVal(rand.Gaus(fracGenZero, sqrt(nGenSignalZeroGamma)/(1.*nGenSignal))) ;
   fracZero.setRange(0., 1.);
   fracOne.setVal(rand.Gaus(fracGenOne, sqrt(nGenSignalOneGamma)/(1.*nGenSignal))) ;
   fracOne.setRange(0., 1.);

   expoConst.setVal(rand.Uniform(-TMath::Abs(5*expoConstGen), TMath::Abs(5*expoConstGen)));
   expoConst.setRange(-TMath::Abs(7*expoConstGen), TMath::Abs(7*expoConstGen));
}



