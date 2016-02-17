#include "fitter_utils.h"
#include "RooBinning.h"
#include "RooRandom.h"


// using namespace std;
// using namespace RooFit;


// void initiateParams(RooArgSet& parset);

// void prepare_PDFs(string workspacename, string trigStr, bool fit2D,
//                   string signalfile, string partrecofile, string combinatorialfile,
//                   string signaltree = "DecayTree", string partrecotree = "DecayTree", string combinatorialtree = "DecayTree");

// void toystudy1DMakeFit(string workspacename,  bool wantplot, 
//       RooHistPdf& histPdfSignalZeroGamma, RooHistPdf& histPdfSignalOneGamma, RooHistPdf& histPdfSignalTwoGamma, 
//                        RooHistPdf& histPdfPartReco, RooMcorMvisTsallis& McorMvis, RooMcorMvisTsallis& McorMvis_fit, 
//       RooAbsPdf::GenSpec& GenSpecSigZeroGamma, RooAbsPdf::GenSpec& GenSpecSigOneGamma, RooAbsPdf::GenSpec& GenSpecSigTwoGamma, 
//       RooAbsPdf::GenSpec& GenSpecPartReco, RooAbsPdf::GenSpec& GenSpecComb, int nGenSignalZeroGamma, int nGenSignalOneGamma,
//       int nGenSignalTwoGamma, int nGenPartReco, int nGenComb, double expoConstGen, RooRealVar& expoConst, ofstream& out, TTree* t, bool update);



void initiateParams(RooArgSet* parset)
{
  RooRealVar *var;
  TIterator *iter = parset->createIterator();

  while((var = (RooRealVar*) iter->Next()))
  {
    if ( !var->isConstant() ) var->randomize();
  }
}


void initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, int nGenPartReco, 
                    int nGenComb, int nGenJpsiLeak, double expoConstGen, RooRealVar& nSignal,
                    RooRealVar& nPartReco, 
                    RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar& expoConst, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma)
{
   TRandom rand;
   rand.SetSeed();



   int nGenSignal = nGenSignalZeroGamma + nGenSignalOneGamma + nGenSignalTwoGamma;
   double fracPartReco2( nGenPartReco/(1.*nGenSignal));
   double nGenSignal2 = rand.Uniform(nGenSignal-5*sqrt(nGenSignal), nGenSignal+5*sqrt(nGenSignal));
   double nGenPartReco2;
   if(!constPartReco) nGenPartReco2 = rand.Uniform(nGenPartReco-5*sqrt(nGenPartReco), nGenPartReco+5*sqrt(nGenPartReco));
   if(constPartReco) nGenPartReco2 = rand.Uniform( (fracPartReco2-5*fracPartRecoSigma.getVal())*nGenSignal2, (fracPartReco2+5*fracPartRecoSigma.getVal())*nGenSignal2);
   double nGenComb2 = rand.Uniform(nGenComb-5*sqrt(nGenComb), nGenComb+5*sqrt(nGenComb));
   double nGenJpsiLeak2 = rand.Uniform(nGenJpsiLeak-5*sqrt(nGenJpsiLeak), nGenJpsiLeak+5*sqrt(nGenJpsiLeak));


   if(!constPartReco)
   {
      nSignal.setVal(nGenSignal2);
      nSignal.setRange(TMath::Max(0.,nGenSignal2-10.*sqrt(nGenSignal)) , nGenSignal2+10*sqrt(nGenSignal));

      nPartReco.setVal(nGenPartReco2);
      nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(nGenPartReco)), nGenPartReco2+10*sqrt(nGenPartReco));
   }

   if(constPartReco)
   {
      
   }

   nComb.setVal(nGenComb2);
   nComb.setRange(TMath::Max(0.,nGenComb2-10.*sqrt(nGenComb)), nGenComb2+10*sqrt(nGenComb));

   nJpsiLeak.setVal(nGenJpsiLeak2);
   nJpsiLeak.setRange(TMath::Max(0., nGenJpsiLeak2-10*sqrt(nGenJpsiLeak)), nGenJpsiLeak2+10*sqrt(nGenJpsiLeak));

   double fracGenZero(nGenSignalZeroGamma/(1.*nGenSignal));
   double fracGenOne(nGenSignalOneGamma/(1.*nGenSignal));


   fracZero.setVal(rand.Gaus(fracGenZero, sqrt(nGenSignalZeroGamma)/(1.*nGenSignal))) ;
   fracZero.setRange(0., 1.);
   fracOne.setVal(rand.Gaus(fracGenOne, sqrt(nGenSignalOneGamma)/(1.*nGenSignal))) ;
   fracOne.setRange(0., 1.);

   expoConst.setVal(rand.Uniform(-TMath::Abs(5*expoConstGen), TMath::Abs(5*expoConstGen)));
   expoConst.setRange(-TMath::Abs(7*expoConstGen), TMath::Abs(7*expoConstGen));
}





void prepare_PDFs(string workspacefile, string trigStr, double BDTcut, bool fit2D,
                  string signalfile, string partrecofile, string combfile, string JpsiLeakfile,
                  double minBMass, double maxBMass,
                  string signaltree, string partrecotree, string combtree, string JpsiLeaktree)
{


  //***********Get the datasets
  TFile* fSignal = new TFile(signalfile.c_str());
  TTree* tSignal = (TTree*)fSignal->Get(signaltree.c_str());
  TFile* fPartReco = new TFile(partrecofile.c_str());
  TTree* tPartReco = (TTree*)fPartReco->Get(partrecotree.c_str());
  TFile* fComb = new TFile(combfile.c_str());
  TTree* tComb = (TTree*)fComb->Get(combtree.c_str()); 
  TFile* fJpsiLeak = new TFile(JpsiLeakfile.c_str());
  TTree* tJpsiLeak = (TTree*)fJpsiLeak->Get(JpsiLeaktree.c_str()); 


  //**********Define variables
  RooRealVar trigVar(trigStr.c_str(), trigStr.c_str(), -10, 10);
  RooRealVar BDTKeeBig2("BDTKeeBig2", "BDTKeeBig2", -1,1);
  RooRealVar B_plus_M("B_plus_M", "M_{visible}", minBMass, maxBMass, "MeV/c^{2}");
  RooRealVar B_plus_M_corr("B_plus_M_corr", "M_{cor}", minBMass+100, 7000, "MeV/c^{2}");
  RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV/c^{2}"); 
  RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
  RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);

  RooRealVar weightPartReco("weightPartReco", "weightPartReco", 0, 10);
  RooRealVar weightLeakage("weightLeakage", "weightLeakage", 0, 10);

  //***********Set Binning
  
 
  RooBinning defaultMBins(floor((maxBMass-minBMass)/(40.)), B_plus_M.getMin(), B_plus_M.getMax() ); 
  RooBinning defaultMCorrBins(floor((7000-minBMass-100)/120.), B_plus_M_corr.getMin(), B_plus_M_corr.getMax()); 
  RooBinning broaderMBins(floor((maxBMass-minBMass)/(80.)), B_plus_M.getMin(), B_plus_M.getMax()); 
  RooBinning broaderMCorrBins(floor((7000-minBMass-100)/240.), B_plus_M_corr.getMin(), B_plus_M_corr.getMax()); 

  B_plus_M.setBinning( defaultMBins);
  B_plus_M_corr.setBinning( defaultMCorrBins );
  B_plus_M.setBinning( broaderMBins, "broaderBins");
  B_plus_M_corr.setBinning( broaderMCorrBins, "broaderBins" );



  B_plus_DTFM_M_zero.setBins(100);

  RooArgSet argset(BDTKeeBig2, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity);
  RooArgSet argsetPartReco(BDTKeeBig2, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightPartReco);
  RooArgSet argsetLeakage(BDTKeeBig2, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightLeakage);

  cout<<"getting the datasets:"<<endl;
  
  RooDataSet* dataSetSignalZeroGamma;
  RooDataSet* dataSetSignalOneGamma;
  RooDataSet* dataSetSignalTwoGamma;
  RooDataSet* dataSetPartReco;
  RooDataSet* dataSetJpsiLeak;
  RooDataSet* dataSetComb;

  TFile* fw;

  dataSetSignalZeroGamma = new RooDataSet("dataSetSignalZeroGamma", "dataSetSignalZeroGamma", tSignal, argset,( " ("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
  dataSetSignalOneGamma = new RooDataSet("dataSetSignalOneGamma", "dataSetSignalOneGamma", tSignal, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
  dataSetSignalTwoGamma = new RooDataSet("dataSetSignalTwoGamma", "dataSetSignalTwoGamma", tSignal, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
  dataSetPartReco = new RooDataSet("dataSetPartReco", "dataSetPartReco",  argsetPartReco, Import(*tPartReco),Cut(("("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("weightPartReco"));

  dataSetJpsiLeak = new RooDataSet("dataSetJpsiLeak", "dataSetJpsiLeak",  argsetLeakage, Import(*tJpsiLeak),Cut(("B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("weightLeakage"));

//  dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig2 > "+d2s(BDTcut-0.15)+")  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
  dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig2 > "+d2s(BDTcut-2)+")  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());

  cout<<"Number of zero: "<< dataSetSignalZeroGamma->sumEntries()<<endl;
  cout<<"Number of one: "<< dataSetSignalOneGamma->sumEntries()<<endl;
  cout<<"Number of two: "<< dataSetSignalTwoGamma->sumEntries()<<endl;
  cout<<"Number of PartReco: "<< dataSetPartReco->sumEntries()<<endl;
  cout<<"Number of Jpsi leaking:"<< dataSetJpsiLeak->sumEntries()<<endl;
  cout<<"Number of combinatorial events:"<< dataSetComb->sumEntries()<<endl;
  

   cout<<"binning the datasets:"<<endl;

   RooArgSet argset2(B_plus_M);
   if (fit2D) argset2.add(B_plus_M_corr);

   RooDataHist dataHistSignalZeroGamma("dataHistSignalZeroGamma", "dataHistSignalZeroGamma", argset2, *dataSetSignalZeroGamma); 
   RooDataHist dataHistSignalOneGamma("dataHistSignalOneGamma", "dataHistSignalOneGamma", argset2, *dataSetSignalOneGamma); 
   RooDataHist dataHistSignalTwoGamma("dataHistSignalTwoGamma", "dataHistSignalTwoGamma", argset2, *dataSetSignalTwoGamma); 
   RooDataHist dataHistComb("dataHistComb", "dataHistComb", argset2, *dataSetComb); 
   RooDataHist dataHistPartReco("dataHistPartReco", "dataHistPartReco", argset2, *dataSetPartReco); 
   RooDataHist dataHistJpsiLeak("dataHistJpsiLeak", "dataHistJpsiLeak", argset2, "broaderBins");
   dataHistJpsiLeak.add(*dataSetJpsiLeak); 

   //*************** Compute Error on J/psi leak

   double ErrorJpsi(0);
   if(dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()) > 0) ErrorJpsi = 1./sqrt(dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()));
   RooRealVar fractionalErrorJpsiLeak("fractionalErrorJpsiLeak", "fractionalErrorJpsiLeak", ErrorJpsi);
   cout<<"JPSI LEAK: "<<dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && (BDTKeeBig2> "+d2s(BDTcut)+") && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   cout<<"JPSI LEAK fractional Error: "<<ErrorJpsi<<endl;


   //***************Create 2D histogram estimates from data


   cout<<"Preparing the 3 2D histPdf: 1";
   //   RooArgSet argset2(B_plus_M);
   RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma,2); cout<<" 2";
   RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma,2); cout<<" 3";
   RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma,2); cout<<" 4";
   RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco,2); cout<<" 5";
   RooHistPdf histPdfJpsiLeak("histPdfJpsiLeak", "histPdfJpsiLeak", argset2, dataHistJpsiLeak,2); cout<<" 6";


   //***************Create combinatorial from fit to data

   RooRealVar expoConst("expoConst", "expoConst", -1e-3, -1, 1);
   RooRealVar T("T", "T", 97, 0, 200);
   RooRealVar n("n", "n", 3.5, 1., 5.5);
   RooAbsPdf *combPDF;

   if (fit2D)
   {  
     combPDF =  new RooMcorMvisTsallis("McorMvis", "McorMvis", B_plus_M_corr, B_plus_M, T, n, expoConst);
   }
   else
   {
      combPDF =  new RooExponential("histPdfComb", "histPdfComb", B_plus_M, expoConst);
   }



   combPDF->fitTo(*dataSetComb); // 

   if (fit2D)
   {
     T.setConstant(true);
     n.setConstant(true);
     std::cout<<"T generated is: "<<T.getVal()<<std::endl;
   }

   RooRealVar trueExp("trueExp","trueExp", expoConst.getVal());


   //***************Save everything on a workspace
   RooWorkspace workspace("workspace", "workspace");
   workspace.import(B_plus_DTFM_M_zero);
   workspace.import(B_plus_M);
   workspace.import(B_plus_M_corr);
   workspace.import(expoConst);
   workspace.import(trueExp);
   workspace.import(T);
   workspace.import(n);
   workspace.import(*dataSetSignalZeroGamma);
   workspace.import(*dataSetSignalOneGamma);
   workspace.import(*dataSetSignalTwoGamma);
   workspace.import(*dataSetPartReco);
   workspace.import(*dataSetComb);
   workspace.import(*dataSetJpsiLeak);
   workspace.import(histPdfSignalZeroGamma);
   workspace.import(histPdfSignalOneGamma);
   workspace.import(histPdfSignalTwoGamma);
   workspace.import(histPdfPartReco);
   workspace.import(histPdfJpsiLeak);
   workspace.import(fractionalErrorJpsiLeak);
   // workspace.import(*combPDF);
   // workspace.importClassCode();
   workspace.writeToFile(workspacefile.c_str());

   delete fComb;
   delete fSignal;
   delete fPartReco;
   delete fw;
   delete combPDF;
   


   delete dataSetSignalZeroGamma; 
   delete dataSetSignalOneGamma;
   delete dataSetSignalTwoGamma;
   delete dataSetPartReco;
   delete dataSetJpsiLeak;
   delete dataSetComb;
}


void generate_and_fit(string workspacename,  bool fit2D, bool wantplot, bool constPartReco,
      int nGenSignal,  int nGenPartReco,  int nGenComb, int nGenJpsiLeak,
      double nGenFracZeroGamma,  double nGenFracOneGamma, double fracPartReco_const,
      ofstream& out, TTree* t, bool update, string plotsfile)
{


   //***************Get the PDFs from the workspace

   TFile *fw = new TFile(workspacename.c_str());   
   RooWorkspace* workspace = (RooWorkspace*)fw->Get("workspace");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *B_plus_M_corr = workspace->var("B_plus_M_corr");
   RooRealVar *T = workspace->var("T");
   RooRealVar *n = workspace->var("n");
   RooRealVar *expoConst = workspace->var("expoConst");
   RooRealVar *trueExp = workspace->var("trueExp");
   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");

   RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
   RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
   RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
   RooHistPdf *histPdfJpsiLeak(0);
   if(nGenJpsiLeak>1) histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");


   RooAbsPdf *combPDF;

   if (fit2D)
   {  
      combPDF =  new RooMcorMvisTsallis("combPDF", "combPDF", *B_plus_M_corr, *B_plus_M, *T, *n, *expoConst);
   }
   else
   {
      combPDF =  new RooExponential("combPDF", "combPDF", *B_plus_M, *expoConst);
   }

   double trueExpConst(trueExp->getValV());
   expoConst->setVal(trueExpConst);


   //***************Prepare generation

   int nGenSignalZeroGamma(floor(nGenFracZeroGamma*nGenSignal));
   int nGenSignalOneGamma(floor(nGenFracOneGamma*nGenSignal));
   int nGenSignalTwoGamma(floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));


   RooArgSet argset2(*B_plus_M);
   if (fit2D) argset2.add(*B_plus_M_corr);

   cout<<"Preparing the generation of events 1";

   RooRandom::randomGenerator()->SetSeed();
   RooAbsPdf::GenSpec* GenSpecSignalZeroGamma = histPdfSignalZeroGamma->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalZeroGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalOneGamma = histPdfSignalOneGamma->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalOneGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalTwoGamma = histPdfSignalTwoGamma->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalTwoGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecPartReco =  histPdfPartReco->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenPartReco)); cout<<" 3 "<<endl;
   RooAbsPdf::GenSpec* GenSpecComb = combPDF->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenComb));
   RooAbsPdf::GenSpec* GenSpecJpsiLeak(0);
   if(nGenJpsiLeak>1) GenSpecJpsiLeak = histPdfJpsiLeak->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenJpsiLeak));


   cout<<"Variable loaded:"<<endl;
   B_plus_M->Print(); expoConst->Print(); //B_plus_DTFM_M_zero->Print();
   if (fit2D) B_plus_M_corr->Print(); 


   //***************Generate some datasets

   cout<<"Generating signal Zero Photon"<<endl;
   RooDataSet* dataGenSignalZeroGamma = histPdfSignalZeroGamma->generate(*GenSpecSignalZeroGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating signal One Photon"<<endl;
   RooDataSet* dataGenSignalOneGamma = histPdfSignalOneGamma->generate(*GenSpecSignalOneGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating signal two Photons"<<endl;
   RooDataSet* dataGenSignalTwoGamma = histPdfSignalTwoGamma->generate(*GenSpecSignalTwoGamma);//(argset, 250, false, true, "", false, true);
   cout<<"Generating combinatorial"<<endl;
   RooDataSet* dataGenComb = combPDF->generate(*GenSpecComb);//(argset, 100, false, true, "", false, true);
   cout<<"Generating PartReco"<<endl;
   RooDataSet* dataGenPartReco = histPdfPartReco->generate(*GenSpecPartReco);//argset, 160, false, true, "", false, true);
   RooDataSet* dataGenJpsiLeak(0);
   if(nGenJpsiLeak>1)
   {
      cout<<"Generating Leaking JPsi"<<endl;
      dataGenJpsiLeak = histPdfJpsiLeak->generate(*GenSpecJpsiLeak);//argset, 160, false, true, "", false, true);
   }

 
   if(wantplot)
   {
      //**************Must get the datasets

      RooDataSet* dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
      RooDataSet* dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
      RooDataSet* dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
      RooDataSet* dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
      RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
      RooDataSet* dataSetJpsiLeak = (RooDataSet*)workspace->data("dataSetJpsiLeak");

      //**************Plot all the different components

      PlotShape(*dataSetSignalZeroGamma, *dataGenSignalZeroGamma, *histPdfSignalZeroGamma, plotsfile, "cSignalZeroGamma", *B_plus_M, *B_plus_M_corr,fit2D);
      PlotShape(*dataSetSignalOneGamma, *dataGenSignalOneGamma, *histPdfSignalOneGamma, plotsfile, "cSignalOneGamma", *B_plus_M, *B_plus_M_corr,fit2D);
      PlotShape(*dataSetSignalTwoGamma, *dataGenSignalTwoGamma, *histPdfSignalTwoGamma, plotsfile, "cSignalTwoGamma", *B_plus_M, *B_plus_M_corr,fit2D);
      PlotShape(*dataSetPartReco, *dataGenPartReco, *histPdfPartReco, plotsfile, "cPartReco", *B_plus_M, *B_plus_M_corr,fit2D);
      PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cComb", *B_plus_M, *B_plus_M_corr,fit2D);
      if(nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeak, *dataGenJpsiLeak, *histPdfJpsiLeak, plotsfile, "cJpsiLeak", *B_plus_M, *B_plus_M_corr,fit2D);
   }

   //***************Merge datasets

   RooDataSet* dataGenTot(dataGenPartReco);
   dataGenTot->append(*dataGenSignalZeroGamma);
   dataGenTot->append(*dataGenSignalOneGamma);
   dataGenTot->append(*dataGenSignalTwoGamma);
   dataGenTot->append(*dataGenComb);
   if(nGenJpsiLeak>1) dataGenTot->append(*dataGenJpsiLeak);

   // RooWorkspace workspaceb("workspaceb", "workspaceb");
   // dataGenTot->SetNameTitle("myexperiment","myexperiment");
   // workspaceb.import(*dataGenTot);
   // workspaceb.writeToFile("my_experiment.root");
   
   

   //**************Prepare fitting function



   RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
   // RooFormulaVar nPartReco("nPartReco", "fracPartReco*nSignal", RooArgList(fracPartReco, nSignal)); // 
   RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);
   // RooRealVar fracPartReco("fracPartReco", "fracPartReco", nGenPartReco/(1.*nGenSignal), 0, 2);
   RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(nPartReco,nSignal));
   RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));


   // RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", RooArgList(*histPdfSignalZeroGamma, *histPdfSignalOneGamma, *histPdfSignalTwoGamma), RooArgList(fracZero, fracOneRec));
   RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", RooArgList(*histPdfSignalZeroGamma, *histPdfSignalOneGamma, *histPdfSignalTwoGamma), RooArgList(fracZero, fracOneRec));

   RooArgList pdfList(histPdfSignal, *histPdfPartReco, *combPDF);
   RooArgList yieldList(nSignal, nPartReco, nComb);
    
   if(nGenJpsiLeak>1)
   {
      pdfList.add(*histPdfJpsiLeak);
      yieldList.add(nJpsiLeak); 
   }
   RooAddPdf totPdf("totPdf", "totPdf", pdfList, yieldList);

   //**************** Constrain the fraction of zero and one photon

   RooRealVar fracZeroConstMean("fracZeroConstMean", "fracZeroConstMean", nGenSignalZeroGamma*1./nGenSignal);
   RooRealVar fracZeroConstSigma("fracZeroConstSigma", "fracZeroConstSigma", sqrt(nGenSignalZeroGamma)/nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroConstMean, fracZeroConstSigma); 

   // RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", nGenSignalOneGamma*1./nGenSignal);
   // RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", sqrt(nGenSignalOneGamma)/nGenSignal);
   RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", nGenSignalOneGamma*1./nGenSignal/(1-fracZeroConstMean.getVal()));
   RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", sqrt(nGenSignalOneGamma)/nGenSignal/(1-fracZeroConstMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneConstMean, fracOneConstSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", nGenPartReco/(1.*nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());
   // RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", (1./(1.*nGenSignal))*sqrt(nGenPartReco + (nGenPartReco*nGenPartReco/(1.*nGenSignal))) ) ;

    RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

    RooRealVar JpsiLeakMean("JpsiLeakMean", "JpsiLeakMean", nGenJpsiLeak);
    RooRealVar JpsiLeakSigma("JpsiLeakSigma", "JpsiLeakSigma", nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
    RooGaussian JpsiLeakConst("JpsiLeakConst", "JpsiLeakConst", nJpsiLeak, JpsiLeakMean, JpsiLeakSigma); 


   //**************** fit

    RooArgSet *par_set = totPdf.getParameters(dataGenTot);
    // initiateParams(par_set);
    initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
                   nGenPartReco, nGenComb, nGenJpsiLeak, trueExpConst, nSignal, nPartReco, nComb, fracZero, fracOne, *expoConst, nJpsiLeak,  constPartReco, fracPartRecoSigma);

    RooArgSet constraints(fracZeroConst, fracOneConst);
    if (constPartReco) constraints.add(fracPartRecoConst);
    if(nGenJpsiLeak>0) constraints.add(JpsiLeakConst);
   //nJpsiLeak.setConstant(true);

    RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(constraints));
    RooMinuit minuit(*nll);
    minuit.setStrategy(2);


    int migradRes(1);
    int hesseRes(4);
    double edm(10);
    int nrefit(0);

    RooFitResult* fitRes(0);


    for(int i(0); (i<10) && ( (migradRes != 0) || (hesseRes !=0) || (edm > 1e-4)); ++i)
    {
      // initiateParams(par_set);
      initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
                     nGenPartReco, nGenComb, nGenJpsiLeak, trueExpConst, nSignal, nPartReco, nComb, fracZero, fracOne, *expoConst, nJpsiLeak, constPartReco, fracPartRecoSigma);
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

   if(wantplot) plot_fit_result(plotsfile, totPdf, *dataGenTot);

   fw->Close();
   //delete and return
   delete nll;
   delete par_set;
   delete workspace;
   delete fw;
   delete combPDF;

}


void plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataSet dataGenTot)
{

  //**************Prepare TFile to save the plots

  TFile f2(plotsfile.c_str(), "UPDATE");
  //**************Plot the results of the fit

  RooArgSet *var_set = totPdf.getObservables(dataGenTot);
  TIterator *iter = var_set->createIterator();
  RooRealVar *var;
  
  std::vector<RooPlot*> plots;
  RooPlot* frame;

  while((var = (RooRealVar*) iter->Next()))
  {
  
    frame = var->frame();
    dataGenTot.plotOn(frame);
    totPdf.plotOn(frame, Components("histPdfPartReco"), LineColor(kBlue));
    totPdf.plotOn(frame, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
    totPdf.plotOn(frame, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
    totPdf.plotOn(frame, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
    totPdf.plotOn(frame, Components("histPdfJpsiLeak"), LineColor(14));
    totPdf.plotOn(frame, Components("combPDF"), LineColor(kBlack));
    totPdf.plotOn(frame, LineColor(kRed));
  
    plots.push_back(frame);
    
  }  
  
  if (!(plots.size())) return;

  TCanvas cFit("cFit", "cFit", 600, 800);
  cFit.Divide(1,2);
  cFit.cd(1); plots[0]->Draw();
  if (plots.size()>1){ 
    cFit.cd(2); plots[1]->Draw();
  }
  
  cFit.Write();
  f2.Close();
    
  
}

void PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& B_plus_M_corr, bool fit2D)
{
   if(fit2D) PlotShape2D(originDataSet, genDataSet, shape, plotsfile, canvName, B_plus_M, B_plus_M_corr);
   if(!fit2D) PlotShape1D(originDataSet, genDataSet, shape, plotsfile, canvName, B_plus_M);
}

void PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& B_plus_M_corr)
{
      //**************Prepare TFile to save the plots
      
      TFile f2(plotsfile.c_str(), "UPDATE");

      //**************Plot Signal Zero Gamma

      TH2F* th2fKey = (TH2F*)shape.createHistogram("th2Shape", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));
      cout<<genDataSet.sumEntries()<<endl;
      TH2F* th2fGen = (TH2F*)genDataSet.createHistogram("th2fGen", B_plus_M, Binning(20), YVar(B_plus_M_corr, Binning(20)));

      RooPlot* plotM = B_plus_M.frame();
      originDataSet.plotOn(plotM);
      shape.plotOn(plotM);

      RooPlot* plotMCorr = B_plus_M_corr.frame();
      originDataSet.plotOn(plotMCorr);
      shape.plotOn(plotMCorr);

      TCanvas canv(canvName.c_str(), canvName.c_str(), 800, 800);
      canv.Divide(2,2);
      canv.cd(1); th2fGen->Draw("lego");
      canv.cd(2); th2fKey->Draw("surf");
      canv.cd(3); plotM->Draw();
      canv.cd(4); plotMCorr->Draw();

      canv.Write();

      f2.Close();
}

void PlotShape1D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M)
{
   TFile f2(plotsfile.c_str(), "UPDATE");

   RooPlot* plotGen = B_plus_M.frame(Binning(20));
   genDataSet.plotOn(plotGen);

   RooPlot* plotM = B_plus_M.frame();
   originDataSet.plotOn(plotM);
   shape.plotOn(plotM);

   TCanvas canv(canvName.c_str(), canvName.c_str(), 800, 800);
   canv.Divide(1,2);
   canv.cd(1); plotGen->Draw();
   canv.cd(2); plotM->Draw();

   canv.Write();

   f2.Close();
}
