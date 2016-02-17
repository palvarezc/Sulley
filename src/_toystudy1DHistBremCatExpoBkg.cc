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

#include "usefulFunctions.h"


using namespace std;
using namespace RooFit;

//void toystudy2DCreateKey(string workspaceName, TFile* fileKey, RooNDKeysPdf& histPdfSignal, RooNDKeysPdf& histPdfPartReco, RooNDKeysPdf& histPdfComb);
vector<double> toystudy1DMakeFit(int indx, string workspacename,  bool wantplot, 
      RooHistPdf& histPdfSignalZeroGamma, RooHistPdf& histPdfSignalOneGamma, RooHistPdf& histPdfSignalTwoGamma, 
      RooHistPdf& histPdfPartReco, RooExponential& histPdfComb, double trueExpConst, 
      RooAbsPdf::GenSpec& GenSpecSigZeroGamma, RooAbsPdf::GenSpec& GenSpecSigOneGamma, RooAbsPdf::GenSpec& GenSpecSigTwoGamma, 
      RooAbsPdf::GenSpec& GenSpecPartReco, RooAbsPdf::GenSpec& GenSpecComb, int nGenSignalZeroGamma, int nGenSignalOneGamma,
      int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, ofstream& out, TTree* t, bool update);
void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartRec, int nGenComb, double expoConstGen, RooRealVar& nSignal,
      RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst2);

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
   TFile* fSignal = new TFile("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_BDT_ctrl_trigged.root");
   TTree* tSignal = (TTree*)fSignal->Get("DecayTree");
   TFile* fPartReco = new TFile("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_BDT_prc_trigged.root");
   TTree* tPartReco = (TTree*)fPartReco->Get("DecayTree");

   TFile* fComb = new TFile("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_piee_trigged.root");
   TTree* tComb = (TTree*)fComb->Get("DecayTree");

   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV/c^{2}"); 
   B_plus_DTFM_M_zero.setBins(100); 
   RooRealVar BDTKeeBig3("BDTKeeBig3", "BDTKeeBig3", -1,1);
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", 4880, 5700, "MeV/c^{2}");
   RooRealVar L0ETOSOnly_d("L0ETOSOnly_d", "L0ETOSOnly_d", -10, 10);
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);
   B_plus_M.setBins(20);
   RooRealVar B_plus_M_comb("B_plus_M", "M_{visible}", 4880, 5700, "MeV/c^{2}");
   B_plus_M_comb.setBins(4);


   RooArgSet argset(BDTKeeBig3, B_plus_DTFM_M_zero, B_plus_M, L0ETOSOnly_d, e_plus_BremMultiplicity, e_minus_BremMultiplicity);

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
      fw = new TFile("output/toystudy1DHistBremCatExpoBkg_L0ETOS.root");
      RooWorkspace* workspace = (RooWorkspace*)fw->Get("workspace");
      dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
      dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
      dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
      dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
      dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
   }


   cout<<"binning the datasets:"<<endl;

   RooArgSet argset2(B_plus_M);
   RooArgSet argset2_comb(B_plus_M_comb);
   RooDataHist dataHistSignalZeroGamma("dataHistSignalZeroGamma", "dataHistSignalZeroGamma", argset2, *dataSetSignalZeroGamma); 
   RooDataHist dataHistSignalOneGamma("dataHistSignalOneGamma", "dataHistSignalOneGamma", argset2, *dataSetSignalOneGamma); 
   RooDataHist dataHistSignalTwoGamma("dataHistSignalTwoGamma", "dataHistSignalTwoGamma", argset2, *dataSetSignalTwoGamma); 
   RooDataHist dataHistComb("dataHistComb", "dataHistComb", argset2_comb, *dataSetComb); 
   RooDataHist dataHistPartReco("dataHistPartReco", "dataHistPartReco", argset2, *dataSetPartReco); 

   //***************Create 1D histogram estimates from data


   cout<<"Preparing the 3 1D histPdf: 1"<<endl;
   //   RooArgSet argset2(B_plus_M);
   RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco,2); cout<<" 3"<<endl;
   RooRealVar expoConst("expoConst", "expoConst", -1e-3, -1, 1);
   RooExponential histPdfComb("histPdfComb", "histPdfComb", B_plus_M, expoConst); 
   histPdfComb.fitTo(*dataSetComb);
   double trueExpConst(expoConst.getValV());

   cout<<"Preparing the generation of events 1";

   RooAbsPdf::GenSpec* GenSpecSignalZeroGamma = histPdfSignalZeroGamma.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalZeroGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalOneGamma = histPdfSignalOneGamma.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalOneGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalTwoGamma = histPdfSignalTwoGamma.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalTwoGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecPartReco =  histPdfPartReco.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenPartReco)); cout<<" 3 "<<endl;
   RooAbsPdf::GenSpec* GenSpecComb = histPdfComb.prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenComb));


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
      workspace.writeToFile("output/toystudy1DHistBremCatExpoBkg_L0ETOS.root");
   }

   //***************Prepare the stuff to generate events


   TFile f("output/toystudy1DHistBremCatExpoBkg_results.root", "UPDATE");
   TTree t("paramsFloatingExpConstrainedPartReco_L0ETOS", "paramsFloatingExpConstrainedPartReco_L0ETOS");

   bool wantPlots(true);
   bool update(false);

   ofstream out("fitResult1D.dat");

   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Generation and fit number "<<i<<endl;
      toystudy1DMakeFit(i,"output/toystudy1DHistBremCatExpoBkg_L0ETOS.root", wantPlots, histPdfSignalZeroGamma, histPdfSignalOneGamma, histPdfSignalTwoGamma, histPdfPartReco, histPdfComb, trueExpConst, *GenSpecSignalZeroGamma,*GenSpecSignalOneGamma, *GenSpecSignalTwoGamma,  *GenSpecPartReco, *GenSpecComb, nGenSignalZeroGamma,nGenSignalOneGamma, nGenSignalTwoGamma, nGenPartReco, nGenComb, out, &t, update);


      wantPlots = true;
      update = true;
   }

   out.close();

   f.cd();
   t.Write();

   ofstream out2("output/resultToy.dat", std::ios::app);
   out2<<"1D Mvis:"<<endl;
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


vector<double> toystudy1DMakeFit(int indx, string workspacename,  bool wantplot, 
      RooHistPdf& histPdfSignalZeroGamma, RooHistPdf& histPdfSignalOneGamma, RooHistPdf& histPdfSignalTwoGamma, 
      RooHistPdf& histPdfPartReco, RooExponential& histPdfComb, double trueExpConst, 
      RooAbsPdf::GenSpec& GenSpecSignalZeroGamma, RooAbsPdf::GenSpec& GenSpecSignalOneGamma, RooAbsPdf::GenSpec& GenSpecSignalTwoGamma, 
      RooAbsPdf::GenSpec& GenSpecPartReco, RooAbsPdf::GenSpec& GenSpecComb, int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, ofstream& out, TTree* t, bool update)
{
   //***************Get the stuff from the workspace

   TFile f(workspacename.c_str());
   RooWorkspace* workspace = (RooWorkspace*)f.Get("workspace"); 

   RooRealVar B_plus_M(*workspace->var("B_plus_M")); 
   RooRealVar B_plus_DTFM_M_zero(*workspace->var("B_plus_DTFM_M_zero")); 


   cout<<"Variable loaded:"<<endl;
   B_plus_M.Print(); B_plus_DTFM_M_zero.Print();


   //***************Generate some datasets

   RooArgSet argset(B_plus_M);
   cout<<"Generating signal Zero Photon"<<endl;
   RooDataSet* dataGenSignalZeroGamma = histPdfSignalZeroGamma.generate(GenSpecSignalZeroGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating signal One Photon"<<endl;
   RooDataSet* dataGenSignalOneGamma = histPdfSignalOneGamma.generate(GenSpecSignalOneGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating signal two Photons"<<endl;
   RooDataSet* dataGenSignalTwoGamma = histPdfSignalTwoGamma.generate(GenSpecSignalTwoGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating combinatorial"<<endl;
   RooDataSet* dataGenComb = histPdfComb.generate(GenSpecComb);//(argset, 100, false, true, "", false, true);
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

      TFile f2("output/plots1DHistBremCatExpoBkgfit.root", "UPDATE");

      //**************Plot Signal ZeroGamma

      RooPlot* plotGenSignalZeroGamma = B_plus_M.frame(Binning(20));
      dataGenSignalZeroGamma->plotOn(plotGenSignalZeroGamma);

      RooPlot* plotMSignalZeroGamma = B_plus_M.frame();
      dataSetSignalZeroGamma->plotOn(plotMSignalZeroGamma);
      histPdfSignalZeroGamma.plotOn(plotMSignalZeroGamma);

      TCanvas* cSignalZeroGamma = new TCanvas("cSignalZeroGamma", "cSignalZeroGamma", 800, 800);
      cSignalZeroGamma->Divide(1,2);
      cSignalZeroGamma->cd(1); plotGenSignalZeroGamma->Draw();
      cSignalZeroGamma->cd(2); plotMSignalZeroGamma->Draw();

      cSignalZeroGamma->Write();

      //**************Plot Signal OneGamma

      RooPlot* plotGenSignalOneGamma = B_plus_M.frame(Binning(20));
      dataGenSignalOneGamma->plotOn(plotGenSignalOneGamma);

      RooPlot* plotMSignalOneGamma = B_plus_M.frame();
      dataSetSignalOneGamma->plotOn(plotMSignalOneGamma);
      histPdfSignalOneGamma.plotOn(plotMSignalOneGamma);

      TCanvas* cSignalOneGamma = new TCanvas("cSignalOneGamma", "cSignalOneGamma", 800, 800);
      cSignalOneGamma->Divide(1,2);
      cSignalOneGamma->cd(1); plotGenSignalOneGamma->Draw();
      cSignalOneGamma->cd(2); plotMSignalOneGamma->Draw();

      cSignalOneGamma->Write();

      //**************Plot Signal TwoGamma

      RooPlot* plotGenSignalTwoGamma = B_plus_M.frame(Binning(20));
      dataGenSignalTwoGamma->plotOn(plotGenSignalTwoGamma);

      RooPlot* plotMSignalTwoGamma = B_plus_M.frame();
      dataSetSignalTwoGamma->plotOn(plotMSignalTwoGamma);
      histPdfSignalTwoGamma.plotOn(plotMSignalTwoGamma);

      TCanvas* cSignalTwoGamma = new TCanvas("cSignalTwoGamma", "cSignalTwoGamma", 800, 800);
      cSignalTwoGamma->Divide(1,2);
      cSignalTwoGamma->cd(1); plotGenSignalTwoGamma->Draw();
      cSignalTwoGamma->cd(2); plotMSignalTwoGamma->Draw();

      cSignalTwoGamma->Write();

      //**************Plot PartReco

      RooPlot* plotGenPartReco = B_plus_M.frame(Binning(20));
      dataGenPartReco->plotOn(plotGenPartReco);

      RooPlot* plotMPartReco = B_plus_M.frame();
      dataSetPartReco->plotOn(plotMPartReco);
      histPdfPartReco.plotOn(plotMPartReco);

      TCanvas* cPartReco = new TCanvas("cPartReco", "cPartReco", 800, 800);
      cPartReco->Divide(1,2);
      cPartReco->cd(1); plotGenPartReco->Draw();
      cPartReco->cd(2); plotMPartReco->Draw();

      cPartReco->Write();

      //**************Plot Comb

      RooPlot* plotGenComb = B_plus_M.frame(Binning(20));
      dataGenComb->plotOn(plotGenComb);

      RooPlot* plotMComb = B_plus_M.frame(Binning(12));
      dataSetComb->plotOn(plotMComb);
      histPdfComb.plotOn(plotMComb);

      TCanvas* cComb = new TCanvas("cComb", "cComb", 800, 800);
      cComb->Divide(1,2);
      cComb->cd(1); plotGenComb->Draw();
      cComb->cd(2); plotMComb->Draw();

      cComb->Write();

   }

   f.cd();

   //***************Merge datasets

   RooDataSet* dataGenTot(dataGenPartReco);
   dataGenTot->append(*dataGenSignalZeroGamma);
   dataGenTot->append(*dataGenSignalOneGamma);
   dataGenTot->append(*dataGenSignalTwoGamma);
   dataGenTot->append(*dataGenComb);

   //**************Prepare fitting function

   double nGenSignal(nGenSignalZeroGamma+nGenSignalOneGamma+nGenSignalTwoGamma);
   RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   RooRealVar fracPartReco("fracPartReco", "fracPartReco", nGenPartReco/(1.*nGenSignal), 0, 2);

   //RooFormulaVar nPartReco("nPartReco", "nPartReco", "fracPartReco*nSignal", RooArgList(fracPartReco, nSignal));
   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);

   RooRealVar expoConst2("expoConst2", "expoConst2",0,-1,1);
   RooRealVar expoConst3("expoConst3", "expoConst3", trueExpConst,-1,1);
   expoConst3.setConstant(true);
   RooExponential expPdf2("expPdf2", "expPdf2", B_plus_M, expoConst2);
   double expoConstGen(-1e-3);
   RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));

   RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", RooArgList(histPdfSignalZeroGamma, histPdfSignalOneGamma, histPdfSignalTwoGamma), RooArgList(fracZero, fracOneRec));
   RooAddPdf totPdf("totPdf", "totPdf", RooArgList(histPdfSignal, histPdfPartReco, expPdf2), RooArgList(nSignal, nPartReco, nComb));

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

   cout<<"CONSTRAIN PARTRECO PARAM: "<<fracPartRecoMean.getValV()<<" "<<fracPartRecoSigma.getValV()<<endl;

   //**************** fit   

   initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, nGenPartReco, nGenComb, expoConstGen, nSignal, nPartReco, nComb, fracZero, fracOne, expoConst2);

   RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(RooArgSet(fracZeroConst, fracOneConst)));
   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(RooArgSet(fracZeroConst, fracOneConst, fracPartRecoConst)));
   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(RooArgSet(fracZeroConst, fracOneConst, fracPartRecoConst)));
   RooMinuit minuit(*nll);
   minuit.setStrategy(2);

   int migradRes(1);
   int hesseRes(4);
   double edm(10);
   int nrefit(0);
   
   RooFitResult* fitRes(0);

   for(int i(0); (i<10) && ( (migradRes != 0) || (hesseRes !=0) || (edm > 1e-4)); ++i)
   {
      initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, nGenPartReco, nGenComb, expoConstGen, nSignal, nPartReco, nComb, fracZero, fracOne, expoConst2);
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

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;



   vector<double> ret(8);

   if(wantplot)
   {
      //**************Prepare TFile to save the plots

      TFile f2("output/plots1DHistBremCatExpoBkgfit.root", "UPDATE");
      //**************Plot the results of the fit

      RooPlot* frameVis = B_plus_M.frame();
      dataGenTot->plotOn(frameVis);
      totPdf.plotOn(frameVis, Components("histPdfPartReco"), LineColor(kBlue));
      totPdf.plotOn(frameVis, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
      totPdf.plotOn(frameVis, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
      totPdf.plotOn(frameVis, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
      totPdf.plotOn(frameVis, Components("expPdf2"), LineColor(kBlack));
      totPdf.plotOn(frameVis, LineColor(kRed));


      TCanvas* cFit = new TCanvas("cFit", "cFit", 600, 600);
      frameVis->Draw();

      //write and erase
      cFit->Write(makestring(indx).c_str());
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

   return ret;
}

void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, double expoConstGen, RooRealVar& nSignal,
      RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst2)
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
   nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(nGenSignal)) , nGenPartReco2+10*sqrt(nGenSignal));

   nComb.setVal(nGenComb2);
   nComb.setRange(TMath::Max(0.,nGenComb2-10.*sqrt(nGenComb)), nGenComb2+10*sqrt(nGenComb));

   double fracGenZero(nGenSignalZeroGamma/(1.*nGenSignal));
   double fracGenOne(nGenSignalOneGamma/(1.*nGenSignal));


   fracZero.setVal(rand.Gaus(fracGenZero, sqrt(nGenSignalZeroGamma)/(1.*nGenSignal))) ;
   fracZero.setRange(0., 1.);
   fracOne.setVal(rand.Gaus(fracGenOne, sqrt(nGenSignalOneGamma)/(1.*nGenSignal))) ;
   fracOne.setRange(0., 1.);

   expoConst2.setVal(rand.Uniform(-TMath::Abs(5*expoConstGen), TMath::Abs(5*expoConstGen)));
   expoConst2.setRange(-TMath::Abs(7*expoConstGen), TMath::Abs(7*expoConstGen));
}

