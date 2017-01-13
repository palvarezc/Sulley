#include "fitter_utils_HistFact.h"
#include "RooBinning.h"
#include "RooRandom.h"
#include "RooExpOfPolyTimesXBinned.h"
#include "RooExpOfPolyTimesX.h"
#include "RooExpOfPoly.h"

#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooRealSumPdf.h"
#include "RooProduct.h"

string measurementName = "my_measurement";
string channelName = "B2Kee";


FitterUtilsHistFact::FitterUtilsHistFact(int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_
                                         , double nGenFracOneGamma_,
                                         string workspacename_)
  :nGenSignal(nGenSignal_), nGenPartReco(nGenPartReco_), nGenComb(nGenComb_), nGenJpsiLeak(nGenJpsiLeak_), 
   nGenFracZeroGamma(nGenFracZeroGamma_), nGenFracOneGamma(nGenFracOneGamma_), 
   workspacename(workspacename_), 
   nMassBins(0), nmisPTBins(40), binWidthMass(40.), binWidthMassBroader(80.), binWidthMassExtended(80.)
{}


void FitterUtilsHistFact::initiateParams(RooArgSet* parset)
{
   RooRealVar *var;
   TIterator *iter = parset->createIterator();

   while((var = (RooRealVar*) iter->Next()))
   {
      if ( !var->isConstant() ) var->randomize();
   }
}


void FitterUtilsHistFact::initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, RooRealVar& nSignal, RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma, 
      RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee, 
      RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen )
{
   TRandom rand;
   rand.SetSeed();

   int nGenSignal = nGenSignalZeroGamma + nGenSignalOneGamma + nGenSignalTwoGamma;

   double nGenSignal2;
   double nGenPartReco2;
   if(!constPartReco)
   {
      nGenSignal2 = rand.Uniform(nGenSignal-5*sqrt(nGenSignal), nGenSignal+5*sqrt(nGenSignal));
      nGenPartReco2 = rand.Uniform(nGenPartReco-5*sqrt(nGenPartReco), nGenPartReco+5*sqrt(nGenPartReco));
   }
   if(constPartReco)
   { 
      double nGenSigPartReco( nGenSignal+nGenPartReco );
      double nGenSigPartReco2( rand.Uniform( nGenSigPartReco-5*sqrt(nGenSigPartReco), nGenSigPartReco+5*sqrt(nGenSigPartReco) ) );
      double fracPartReco1( nGenPartReco/(1.*nGenSignal));
      double fracPartReco2( rand.Uniform(fracPartReco1-5*fracPartRecoSigma.getVal(), fracPartReco1+5*fracPartRecoSigma.getVal()) ); 

      nGenPartReco2 = fracPartReco2*nGenSigPartReco2 / (1+fracPartReco2); 
      nGenSignal2 = nGenSigPartReco2 / (1+fracPartReco2); 
   }
   double nGenComb2 = rand.Uniform(nGenComb-5*sqrt(nGenComb), nGenComb+5*sqrt(nGenComb));
   double nGenJpsiLeak2 = rand.Uniform(nGenJpsiLeak-5*sqrt(nGenJpsiLeak), nGenJpsiLeak+5*sqrt(nGenJpsiLeak));


   nSignal.setVal(nGenSignal2);
   nSignal.setRange(TMath::Max(0.,nGenSignal2-10.*sqrt(nGenSignal)) , nGenSignal2+10*sqrt(nGenSignal));

   nPartReco.setVal(nGenPartReco2);
   nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(nGenPartReco)), nGenPartReco2+10*sqrt(nGenPartReco));


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

   l1Kee.setVal(rand.Uniform( l1KeeGen.getVal() - 5*l1KeeGen.getError(), l1KeeGen.getVal() + 5*l1KeeGen.getError() ) );
   l1Kee.setRange( l1KeeGen.getVal() - 10*l1KeeGen.getError(), l1KeeGen.getVal() + 10*l1KeeGen.getError() );

   l2Kee.setVal(rand.Uniform( l2KeeGen.getVal() - 5*l2KeeGen.getError(), l2KeeGen.getVal() + 5*l2KeeGen.getError() ) );
   l2Kee.setRange( l2KeeGen.getVal() - 10*l2KeeGen.getError(), l2KeeGen.getVal() + 10*l2KeeGen.getError() );

   l3Kee.setVal(rand.Uniform( l3KeeGen.getVal() - 5*l3KeeGen.getError(), l3KeeGen.getVal() + 5*l3KeeGen.getError() ) );
   l3Kee.setRange( l3KeeGen.getVal() - 10*l3KeeGen.getError(), l3KeeGen.getVal() + 10*l3KeeGen.getError() );

   l4Kee.setVal(rand.Uniform( l4KeeGen.getVal() - 5*l4KeeGen.getError(), l4KeeGen.getVal() + 5*l4KeeGen.getError() ) );
   l4Kee.setRange( l4KeeGen.getVal() - 10*l4KeeGen.getError(), l4KeeGen.getVal() + 10*l4KeeGen.getError() );

   l5Kee.setVal(rand.Uniform( l5KeeGen.getVal() - 5*l5KeeGen.getError(), l5KeeGen.getVal() + 5*l5KeeGen.getError() ) );
   l5Kee.setRange( l5KeeGen.getVal() - 10*l5KeeGen.getError(), l5KeeGen.getVal() + 10*l5KeeGen.getError() );
}


TH2D* FitterUtilsHistFact::make_base_histogram(string componentName, string componentFile, string componentTree,
                                              string cut, vector<string> varset, RooRealVar &x, RooRealVar &y, 
                                              string binningName)
{  


   //***********Get the datasets
   TFile* fcomponent = new TFile(componentFile.c_str());
   TTree* tcomponent = (TTree*)fcomponent->Get(componentTree.c_str());
   
   TH2D *hcomponent = new TH2D(componentName.c_str(),componentName.c_str(),
                               x.getBins(binningName.c_str()), x.getMin(binningName.c_str()), x.getMax(binningName.c_str()),
                               y.getBins(binningName.c_str()), y.getMin(binningName.c_str()), y.getMax(binningName.c_str()));

   // TH2D hcomponent(componentName.c_str(),componentName.c_str(),
   //                             x.getBins(binningName.c_str()), x.getMin(binningName.c_str()), x.getMax(binningName.c_str()),
   //                             y.getBins(binningName.c_str()), y.getMin(binningName.c_str()), y.getMax(binningName.c_str()));
   
                                        
   // cout<<x.GetName()<<"["<<x.getBins(binningName.c_str())<<","<<x.getMin(binningName.c_str())<<","<<x.getMax(binningName.c_str())<<endl;
   // cout<<y.GetName()<<"["<<y.getBins(binningName.c_str())<<","<<y.getMin(binningName.c_str())<<","<<y.getMax(binningName.c_str())<<endl;
  

   // TH2D* hcomponent = (TH2D*) x.createHistogram(componentName.c_str(),YVar(y),Binning(binningName.c_str()));
   hcomponent->Sumw2();

   //***********Set only variables needed

   tcomponent->SetBranchStatus("*",0);
   for(vector<string>::const_iterator i = varset.begin(); i != varset.end(); ++i) {
     tcomponent->SetBranchStatus((*i).c_str(), 1);
     // cout <<"Activating.. "<< *i << endl; 
   }

   string histcontent = y.GetName();
   histcontent += ":";
   histcontent += x.GetName();
   histcontent += ">>";
   histcontent += hcomponent->GetName();
   tcomponent->Draw(histcontent.c_str(),cut.c_str(),"qoff");
   
   // cout<<"Integral "<<componentName<<" : "<<tcomponent->GetEntries(cut.c_str())<<endl;
   // cout<<"Integral "<<componentName<<" : "<<hcomponent.GetSumOfWeights()<<endl;
   
 
   hcomponent->SetDirectory(0);
   delete tcomponent;
   delete fcomponent;

   return hcomponent;
   
   

}




void FitterUtilsHistFact::fit_combinatorial(string combfile, string combtree, string combcuts, vector<string> varset, 
                                            RooRealVar &B_plus_M, RooRealVar &misPT, RooWorkspace *workspace)
{

   //****************Get data

   TFile* fComb = new TFile(combfile.c_str());
   TTree* tComb = (TTree*)fComb->Get(combtree.c_str()); 

   RooRealVar *var;
   RooArgSet argset;
   
   tComb->SetBranchStatus("*",0);
   for(vector<string>::const_iterator i = varset.begin(); i != varset.end(); ++i) {
     tComb->SetBranchStatus((*i).c_str(), 1);
     if (*i=="B_plus_M" || *i=="misPT") continue;
     var = new RooRealVar((*i).c_str(),(*i).c_str(),-10,10);
     argset.add(*var);
   }
   argset.add(B_plus_M);
   argset.add(misPT);
   

   RooDataSet *dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, combcuts.c_str());


   //****************Define model
   RooRealVar l1Kee("l1Kee", "l1Kee", +2.3e-3, -1e-1, 1e-1);
   RooRealVar l2Kee("l2Kee", "l2Kee", +4.9e-3, -1e-1, 1e-1);
   RooRealVar l3Kee("l3Kee", "l3Kee", -4.1e-7, -1e-5, 1e-5);
   RooRealVar l4Kee("l4Kee", "l4Kee", -4e-7, -1e-5, 1e-5);
   RooRealVar l5Kee("l5Kee", "l5Kee", 1e-10, -1e-9, 1e-9);

   RooBinning m_binning(nMassBins,B_plus_M.getMin(),B_plus_M.getMax(),"m_binning");
   RooBinning pT_binning(nmisPTBins,misPT.getMin(),misPT.getMax(),"pT_binning");   

   RooAbsPdf *combPDF;
   // combPDF =  new RooExpOfPolyTimesXBinned("combPDF", "combPDF", B_plus_M, misPT, l1Kee, l2Kee, l3Kee, l4Kee, l5Kee, 
   //                                         m_binning, pT_binning);
   combPDF =  new RooExpOfPolyTimesX("combPDF", "combPDF", B_plus_M, misPT, l1Kee, l2Kee, l3Kee, l4Kee, l5Kee);

   combPDF->fitTo(*dataSetComb); // 


   //****************The following is just needed to estimate systematics

   
   RooRealVar l1KeeSyst("l1KeeSyst", "l1KeeSyst", +2.3e-3, -1e-1, 1e-1);
   RooRealVar l2KeeSyst("l2KeeSyst", "l2KeeSyst", +4.9e-3, -1e-1, 1e-1);
   RooRealVar l3KeeSyst("l3KeeSyst", "l3KeeSyst", -4.1e-7, -1e-5, 1e-5);
   RooRealVar l4KeeSyst("l4KeeSyst", "l4KeeSyst", -4e-7, -1e-5, 1e-5);
   RooRealVar l5KeeSyst("l5KeeSyst", "l5KeeSyst", 1e-10, -1e-9, 1e-9);
   RooExpOfPoly systPDF("systPDF", "systPDF", B_plus_M, misPT, l1KeeSyst, l2KeeSyst, l3KeeSyst, l4KeeSyst, l5KeeSyst); 
   systPDF.fitTo(*dataSetComb); // 

   l1KeeSyst.setConstant(true);
   l2KeeSyst.setConstant(true);
   l3KeeSyst.setConstant(true);
   l4KeeSyst.setConstant(true);
   l5KeeSyst.setConstant(true);



   //****************End of systematics


   RooRealVar l1KeeGen("l1KeeGen", "l1KeeGen", l1Kee.getVal(), -10, 10);
   l1KeeGen.setError(l1Kee.getError());
   l1KeeGen.setConstant(true);
   RooRealVar l2KeeGen("l2KeeGen", "l2KeeGen", l2Kee.getVal(), -10, 10);
   l2KeeGen.setError(l2Kee.getError());
   l2KeeGen.setConstant(true);
   RooRealVar l3KeeGen("l3KeeGen", "l3KeeGen", l3Kee.getVal(), -10, 10);
   l3KeeGen.setError(l3Kee.getError());
   l3KeeGen.setConstant(true);
   RooRealVar l4KeeGen("l4KeeGen", "l4KeeGen", l4Kee.getVal(), -10, 10);
   l4KeeGen.setError(l4Kee.getError());
   l4KeeGen.setConstant(true);
   RooRealVar l5KeeGen("l5KeeGen", "l5KeeGen", l5Kee.getVal(), -10, 10);
   l5KeeGen.setError(l5Kee.getError());
   l5KeeGen.setConstant(true);


   //***************Save everything on a workspace
   workspace->import(l1Kee);
   workspace->import(l2Kee);
   workspace->import(l3Kee);
   workspace->import(l4Kee);
   workspace->import(l5Kee);
   workspace->import(l1KeeGen);
   workspace->import(l2KeeGen);
   workspace->import(l3KeeGen);
   workspace->import(l4KeeGen);
   workspace->import(l5KeeGen);
   workspace->import(l1KeeSyst);
   workspace->import(l2KeeSyst);
   workspace->import(l3KeeSyst);
   workspace->import(l4KeeSyst);
   workspace->import(l5KeeSyst);

   
   delete fComb;
   delete combPDF;

   delete dataSetComb;



}


void FitterUtilsHistFact::prepare_PDFs(string trigStr, string weightStr, string BDTVar, double BDTcut,
      string signalfile, string partrecofile, string combfile, string JpsiLeakfile,
      double minBMass, double maxBMass,
      string signaltree, string partrecotree, string combtree, string JpsiLeaktree)
{

  bool templateStat = true;
  

  //**********Define variables
  
  // RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV"); 
  RooRealVar B_plus_M("B_plus_M", "M_{visible}", minBMass, maxBMass, "MeV");
  RooRealVar misPT("misPT", "p_{#perp}", 0, 5000, "MeV");


  //***********Set Binning
  nMassBins = floor((maxBMass-minBMass)/(binWidthMass));
  RooBinning defaultMBins(floor((maxBMass-minBMass)/(binWidthMass)), B_plus_M.getMin(), B_plus_M.getMax() ); 
  RooBinning defaultMisPTBins(nmisPTBins, misPT.getMin(), misPT.getMax()); 
  RooBinning broaderMBins(floor((maxBMass-minBMass)/(binWidthMassBroader)), B_plus_M.getMin(), B_plus_M.getMax()); 
  RooBinning broaderMisPTBins(nmisPTBins, misPT.getMin(), misPT.getMax()); 
  RooBinning extendedMBins(floor((maxBMass-minBMass)/(binWidthMassExtended)), B_plus_M.getMin(), B_plus_M.getMax()); 
  RooBinning extendedMisPTBins(nmisPTBins, -2000,5000); 

  B_plus_M.setBinning( defaultMBins);
  misPT.setBinning( defaultMisPTBins );
  B_plus_M.setBinning( defaultMBins, "default");
  misPT.setBinning( defaultMisPTBins, "default");
  B_plus_M.setBinning( broaderMBins, "broaderBins");
  misPT.setBinning( broaderMisPTBins, "broaderBins" );
  B_plus_M.setBinning( extendedMBins, "extended");
  misPT.setBinning( extendedMisPTBins, "extended" );

  // B_plus_DTFM_M_zero.setBins(100);
  // RooArgSet obsset(B_plus_M, misPT);
  
  
  //**********Define components

  cout<<"Creating histograms..."<<endl;

  string BDTCutString  = "1"; //"("+BDTVar+">"+d2s(BDTcut)+")";
  string TrigCutString = "("+trigStr+"  > 0.9)";
  string MassCutString = "B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass);
  
  string ZeroGamma = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5)";
  string OneGamma  = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5)";
  string TwoGamma  = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5)";
  
  string SignalZeroCuts = "("+BDTCutString+" && "+TrigCutString+" && "+MassCutString+" && "+ZeroGamma+")*"+weightStr;
  string SignalOneCuts = "("+BDTCutString+" && "+TrigCutString+" && "+MassCutString+" && "+OneGamma+")*"+weightStr;
  string SignalTwoCuts = "("+BDTCutString+" && "+TrigCutString+" && "+MassCutString+" && "+TwoGamma+")*"+weightStr;
  string PartRecoCuts = "("+BDTCutString+" && "+TrigCutString+" && "+MassCutString+")*"+weightStr;
  string JpsiLeakCuts = "("+BDTCutString+" && "+TrigCutString+" && "+MassCutString+")*"+weightStr;
  string CombCuts  = BDTCutString+" && "+TrigCutString+" && "+MassCutString+" && K_Kst_isMuon == 0";   
  
  vector<string> argsetComb = {"B_plus_DTFM_M_zero", "misPT", "B_plus_M", 
                               BDTVar, trigStr, "e_plus_BremMultiplicity", "e_minus_BremMultiplicity", "K_Kst_isMuon"};
  vector<string> argsetPartReco = {"B_plus_DTFM_M_zero", "misPT",  "B_plus_M", 
                                    BDTVar, trigStr, "e_plus_BremMultiplicity", "e_minus_BremMultiplicity", weightStr};
  vector<string> argsetJpsiLeak = {"B_plus_DTFM_M_zero", "misPT",  "B_plus_M", 
                                   BDTVar, trigStr, "e_plus_BremMultiplicity", "e_minus_BremMultiplicity", weightStr};
  vector<string> argsetSignal = {"misPT",  "B_plus_M", 
                                 BDTVar, trigStr, "e_plus_BremMultiplicity", "e_minus_BremMultiplicity", weightStr};
  
  // TH2D h_SignalZero = make_base_histogram("SignalZero", signalfile, signaltree, SignalZeroCuts, argsetSignal, B_plus_M, misPT, "default");
  // TH2D h_SignalOne = make_base_histogram("SignalOne", signalfile, signaltree, SignalOneCuts, argsetSignal, B_plus_M, misPT, "default");
  // TH2D h_SignalTwo = make_base_histogram("SignalTwo", signalfile, signaltree, SignalTwoCuts, argsetSignal, B_plus_M, misPT, "default");
  // TH2D h_PartReco = make_base_histogram("PartReco", partrecofile, partrecotree, PartRecoCuts, argsetPartReco, B_plus_M, misPT, "default");
  // TH2D h_JpsiLeak = make_base_histogram("JpsiLeak", JpsiLeakfile, JpsiLeaktree, JpsiLeakCuts, argsetJpsiLeak, B_plus_M, misPT, "default");
  // TH2D h_JpsiLeakbroader = make_base_histogram("JpsiLeakbroader", JpsiLeakfile, JpsiLeaktree, JpsiLeakCuts, argsetJpsiLeak, B_plus_M, misPT, "broaderBins");
  // TH2D h_Comb  = make_base_histogram("Comb", combfile, combtree, CombCuts, argsetComb, B_plus_M, misPT, "default");
  
    
  // misPT.setRange(-2000, 5000);
  // TH2D h_CombExt = make_base_histogram("CombExt", combfile, combtree, CombCuts, argsetComb,  B_plus_M, misPT, "extended");
  // misPT.setRange(0, 5000);


  TH2D* h_SignalZero = make_base_histogram("SignalZero", signalfile, signaltree, SignalZeroCuts, argsetSignal, B_plus_M, misPT, "default");
  TH2D* h_SignalOne = make_base_histogram("SignalOne", signalfile, signaltree, SignalOneCuts, argsetSignal, B_plus_M, misPT, "default");
  TH2D* h_SignalTwo = make_base_histogram("SignalTwo", signalfile, signaltree, SignalTwoCuts, argsetSignal, B_plus_M, misPT, "default");
  TH2D* h_PartReco = make_base_histogram("PartReco", partrecofile, partrecotree, PartRecoCuts, argsetPartReco, B_plus_M, misPT, "default");
  TH2D* h_JpsiLeak = make_base_histogram("JpsiLeak", JpsiLeakfile, JpsiLeaktree, JpsiLeakCuts, argsetJpsiLeak, B_plus_M, misPT, "default");
  TH2D* h_JpsiLeakbroader = make_base_histogram("JpsiLeakbroader", JpsiLeakfile, JpsiLeaktree, JpsiLeakCuts, argsetJpsiLeak, B_plus_M, misPT, "broaderBins");
  TH2D* h_Comb  = make_base_histogram("Comb", combfile, combtree, CombCuts, argsetComb, B_plus_M, misPT, "default");
  
  cout<<"h_SignalZero bins: "<<h_SignalZero->GetNbinsX()*h_SignalZero->GetNbinsY()<<endl;
  cout<<"h_SignalOne bins: "<<h_SignalOne->GetNbinsX()*h_SignalOne->GetNbinsY()<<endl;
  cout<<"h_SignalTwo bins: "<<h_SignalTwo->GetNbinsX()*h_SignalTwo->GetNbinsY()<<endl;
  cout<<"h_PartReco bins: "<<h_PartReco->GetNbinsX()*h_PartReco->GetNbinsY()<<endl;
  cout<<"h_JpsiLeak bins: "<<h_JpsiLeak->GetNbinsX()*h_JpsiLeak->GetNbinsY()<<endl;
  cout<<"h_Comb bins: "<<h_Comb->GetNbinsX()*h_Comb->GetNbinsY()<<endl;
  

  misPT.setRange(-2000, 5000);
  TH2D* h_CombExt = make_base_histogram("CombExt", combfile, combtree, CombCuts, argsetComb,  B_plus_M, misPT, "extended");
  misPT.setRange(0, 5000);
  
 double ErrorJpsi(0);
 cout<<"Number of JpsiLeak : "<< h_JpsiLeak->Integral() <<endl;
 cout<<"Number of SignalZero: "<< h_SignalZero->Integral() <<endl;
 cout<<"Number of SignalOne: "<< h_SignalOne->Integral() <<endl;
 cout<<"Number of SignalTwo: "<< h_SignalTwo->Integral() <<endl;
 cout<<"Number of PartReco: "<< h_PartReco->Integral() <<endl;
 cout<<"Number of Comb: "<< h_Comb->Integral() <<endl;
 if (h_JpsiLeak->Integral()>0) ErrorJpsi = 1./sqrt(h_JpsiLeak->Integral());
  
  
  //**************Create PDFs
  
  cout<<"Creating pdfs..."<<endl;

  double fzero = 0.368;
  double fone = 0.484;
  double ftwo = 1- fzero - fone;
  
  RooStats::HistFactory::Measurement meas(measurementName.c_str(),measurementName.c_str());
  string measPath = workspacename;
  measPath.resize(measPath.find("/"));
  meas.SetOutputFilePrefix(measPath+"/");
  meas.SetExportOnly(kTRUE);
  
  meas.SetPOI("nSignal");

  meas.SetLumi(1.0);
  meas.SetLumiRelErr(0.05);  
  
  RooStats::HistFactory::Channel chan(channelName.c_str());
  chan.SetStatErrorConfig(1e-5,"Poisson");

  //Setup the zero brem 
  RooStats::HistFactory::Sample sigzero("SignalZeroSample");
  sigzero.SetHisto(h_SignalZero);
  if(templateStat) sigzero.ActivateStatError();
  sigzero.SetNormalizeByTheory(kFALSE);
  sigzero.AddNormFactor("nSignal", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
  sigzero.AddNormFactor("fracZero", fzero,0.,1.);
  sigzero.AddNormFactor("mcNorm_sigzero", 1./h_SignalZero->Integral(), 1e-9, 1.);
  // sigzero.AddOverallSys("fracZeroConstraint",0.9,1.1); // Relative to fracZero
  chan.AddSample(sigzero);


  //Setup the one brem 
  RooStats::HistFactory::Sample sigone("SignalOneSample");
  if(templateStat) sigone.ActivateStatError();
  sigone.SetHisto(h_SignalOne);
  sigone.SetNormalizeByTheory(kFALSE);
  sigone.AddNormFactor("nSignal",  1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
  sigone.AddNormFactor("fracOne", fone,0.,1.);
  sigone.AddNormFactor("mcNorm_sigone", 1./h_SignalOne->Integral(), 1e-9, 1.);
  // sigone.AddOverallSys("fracOneConstraint",0.9,1.1); // Relative to fracOne
  chan.AddSample(sigone);


  //Setup the two brem 
  RooStats::HistFactory::Sample sigtwo("SignalTwoSample");
  if(templateStat) sigtwo.ActivateStatError();
  sigtwo.SetHisto(h_SignalTwo);
  sigtwo.SetNormalizeByTheory(kFALSE);
  sigtwo.AddNormFactor("nSignal",  1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
  sigtwo.AddNormFactor("fracTwo", ftwo,0.,1.);
  sigtwo.AddNormFactor("mcNorm_sigtwo", 1./h_SignalTwo->Integral(), 1e-9, 1.);
  // sigtwo.AddOverallSys("fracTwoConstraint",0.9,1.1); // Relative to fracTwo (need to figure out how to lock it to 1-zero-one)
  chan.AddSample(sigtwo);
  

  //Setup the part reco
  RooStats::HistFactory::Sample partreco("PartRecoSample");
  if(templateStat) partreco.ActivateStatError();
  partreco.SetHisto(h_PartReco);
  partreco.SetNormalizeByTheory(kFALSE);
  partreco.AddNormFactor("nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
  partreco.AddNormFactor("mcNorm_partreco", 1./h_PartReco->Integral(), 1e-9, 1.);
  // if (constPartReco) partreco.AddOverallSys("PartRecoConstraint",0.9,1.1); // Relative to nPartReco
  chan.AddSample(partreco);
  

  //Setup the jpsi leak
  RooStats::HistFactory::Sample jpsileak("JpsiLeakSample");
  if(templateStat) jpsileak.ActivateStatError();
  jpsileak.SetHisto(h_JpsiLeak);
  jpsileak.SetNormalizeByTheory(kFALSE);
  jpsileak.AddNormFactor("nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
  jpsileak.AddNormFactor("mcNorm_jpsileak", 1./h_JpsiLeak->Integral(), 1e-9, 1.);
  // if (constJpsiLeak) jpsileak.AddOverallSys("JpsiLeakConstraint",0.9,1.1); // Relative to nJpsiLeak
  chan.AddSample(jpsileak); 


  meas.AddChannel(chan);

  // meas.CollectHistograms();


  //Create the workspace that is going to save the config
  RooWorkspace *workspace;
  workspace = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
  // workspace = new RooWorkspace("combined");
  
  
  ModelConfig *mc = (ModelConfig*) workspace->obj("ModelConfig"); // Get model (noComb) manually
   
  // Lets tell roofit the right names for our histogram variables
  RooArgSet *obs = (RooArgSet*) mc->GetObservables();
  RooRealVar *B_plus_M_mod = (RooRealVar*) obs->find(("obs_x_"+channelName).c_str());
  RooRealVar *misPT_mod = (RooRealVar*) obs->find(("obs_y_"+channelName).c_str());
  // B_plus_M_mod->SetName("B_plus_M");
  B_plus_M_mod->SetTitle("B_plus_M");
  B_plus_M_mod->setUnit("MeV/c^{2}");
  misPT_mod->SetTitle("misPT");
  // misPT_mod->SetName("misPT");
  misPT_mod->setUnit("MeV/c");



  //***************Create combinatorial from fit to data
  // fit_combinatorial(combfile, combtree, CombCuts, argsetComb,
  //                   *B_plus_M_mod, *misPT_mod, workspace);

  fit_combinatorial(combfile, combtree, CombCuts, argsetComb,
                    B_plus_M, misPT, workspace);
  
  RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
  workspace->import(nComb);  


  //***************Save some extra factors
  workspace->factory(Form("nMC_SignalZero[%f]",h_SignalZero->Integral()));
  workspace->factory(Form("nMC_SignalOne[%f]",h_SignalOne->Integral()));
  workspace->factory(Form("nMC_SignalTwo[%f]",h_SignalTwo->Integral()));
  workspace->factory(Form("nMC_PartReco[%f]",h_PartReco->Integral()));
  workspace->factory(Form("nMC_JpsiLeak[%f]",h_JpsiLeak->Integral()));
  workspace->factory(Form("nMC_Comb[%f]",h_Comb->Integral()));
  workspace->factory(Form("nMC_CombExt[%f]",h_CombExt->Integral()));



  //***************Create RooHistPdfs for generation
  
  RooArgList arglist(*B_plus_M_mod, *misPT_mod);
  // RooArgList arglist(B_plus_M, misPT);

  RooDataHist dataHistSignalZeroGamma("dataHistSignalZeroGamma", "dataHistSignalZeroGamma", arglist, h_SignalZero); 
  RooDataHist dataHistSignalOneGamma("dataHistSignalOneGamma", "dataHistSignalOneGamma", arglist, h_SignalOne); 
  RooDataHist dataHistSignalTwoGamma("dataHistSignalTwoGamma", "dataHistSignalTwoGamma", arglist, h_SignalTwo); 
  RooDataHist dataHistComb("dataHistComb", "dataHistComb", arglist, h_Comb); 
  RooDataHist dataHistPartReco("dataHistPartReco", "dataHistPartReco", arglist, h_PartReco); 
  RooDataHist dataHistJpsiLeak("dataHistJpsiLeak", "dataHistJpsiLeak", arglist, h_JpsiLeakbroader); 

   cout<<"Preparing the 3 2D histPdf: 1";
   RooArgSet argset2(*B_plus_M_mod, *misPT_mod);
   // RooArgSet argset2(B_plus_M, misPT);
   RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma,2); cout<<" 2";
   RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma,2); cout<<" 3";
   RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma,2); cout<<" 4";
   RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco,2); cout<<" 5";
   RooHistPdf histPdfJpsiLeak("histPdfJpsiLeak", "histPdfJpsiLeak", argset2, dataHistJpsiLeak,2); cout<<" 6";
   RooHistPdf histPdfComb("histPdfComb", "histPdfComb", argset2, dataHistComb,2);

   workspace->import(histPdfSignalZeroGamma);
   workspace->import(histPdfSignalOneGamma);
   workspace->import(histPdfSignalTwoGamma);
   workspace->import(histPdfPartReco);
   workspace->import(histPdfJpsiLeak);
   workspace->import(histPdfComb);

   workspace->factory(Form("fractionalErrorJpsiLeak[%f]",ErrorJpsi));
   workspace->writeToFile(workspacename.c_str());


  // deleted by ~Sample
  // h_SignalZero->SetDirectory(0);
  // h_SignalOne->SetDirectory(0);
  // h_SignalTwo->SetDirectory(0);
  // h_PartReco->SetDirectory(0);
  // h_JpsiLeak->SetDirectory(0);
  // h_JpsiLeakbroader->SetDirectory(0);
  // h_Comb->SetDirectory(0);
  // h_CombExt->SetDirectory(0);


  // delete h_SignalZero;
  // delete h_SignalOne;
  // delete h_SignalTwo;
  // delete h_PartReco;
  // delete h_JpsiLeak;
  // delete h_JpsiLeakbroader;
  // delete h_Comb;
  // delete h_CombExt;

   workspace->Print();
   
   cout<<"PDFs prepared!"<<endl;
   delete workspace;  


}


void FitterUtilsHistFact::generate(bool wantPlots, string plotsfile)
{
   //***************Get the PDFs from the workspace

   // TFile fw(workspacename.c_str(), "UPDATE");   
   // RooWorkspace* workspace = (RooWorkspace*)fw.Get("combined");   
  TFile* fw = TFile::Open(workspacename.c_str(), "UPDATE");   
   RooWorkspace* workspace = (RooWorkspace*) fw->Get("combined");   
   RooRealVar *l1KeeGen = workspace->var("l1KeeGen");
   RooRealVar *l2KeeGen = workspace->var("l2KeeGen");
   RooRealVar *l3KeeGen = workspace->var("l3KeeGen");
   RooRealVar *l4KeeGen = workspace->var("l4KeeGen");
   RooRealVar *l5KeeGen = workspace->var("l5KeeGen");

   cout<<"HELLO0"<<endl;
   // Get observables
   ModelConfig *mc = (ModelConfig*) workspace->obj("ModelConfig");
   RooArgSet *obs = (RooArgSet*) mc->GetObservables();
   RooRealVar *B_plus_M = (RooRealVar*) obs->find(("obs_x_"+channelName).c_str());
   RooRealVar *misPT = (RooRealVar*) obs->find(("obs_y_"+channelName).c_str());
   // B_plus_M->SetName("B_plus_M");
   B_plus_M->SetTitle("B_plus_M");
   B_plus_M->setUnit("MeV/c^{2}");
   // misPT->SetName("misPT");
   misPT->SetTitle("misPT");
   misPT->setUnit("MeV/c");

   // Get generation pdfs
   RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
   RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
   RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
   RooHistPdf *histPdfComb = (RooHistPdf *) workspace->pdf("histPdfComb");
   RooHistPdf *histPdfJpsiLeak(0);
   cout<<"HELLO0a"<<endl;
   histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");


   cout<<"Here"<<endl;
   RooAbsPdf *combPDF;
   RooBinning m_binning(nMassBins,B_plus_M->getMin(),B_plus_M->getMax(),"m_binning");
   RooBinning pT_binning(nmisPTBins,misPT->getMin(),misPT->getMax(),"pT_binning");   
   // combPDF =  new RooExpOfPolyTimesXBinned("combPDF", "combPDF", *B_plus_M, *misPT, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen, m_binning, pT_binning);

   combPDF =  new RooExpOfPolyTimesX("combPDF", "combPDF", *B_plus_M, *misPT, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);


   //***************Prepare generation

   int nGenSignalZeroGamma(floor(nGenFracZeroGamma*nGenSignal));
   int nGenSignalOneGamma(floor(nGenFracOneGamma*nGenSignal));
   int nGenSignalTwoGamma(floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));
   cout<<"HELLO2"<<endl;


   RooArgSet argset2(*B_plus_M, *misPT);

   cout<<"Preparing the generation of events 1";

   RooRandom::randomGenerator()->SetSeed();
   RooAbsPdf::GenSpec* GenSpecSignalZeroGamma = histPdfSignalZeroGamma->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalZeroGamma)); cout<<" 2 ";
   RooAbsPdf::GenSpec* GenSpecSignalOneGamma = histPdfSignalOneGamma->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalOneGamma)); cout<<" 3 ";
   RooAbsPdf::GenSpec* GenSpecSignalTwoGamma = histPdfSignalTwoGamma->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenSignalTwoGamma)); cout<<" 4 ";
   RooAbsPdf::GenSpec* GenSpecPartReco =  histPdfPartReco->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenPartReco)); cout<<" 5 "<<endl;
   RooAbsPdf::GenSpec* GenSpecComb = combPDF->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenComb));
   RooAbsPdf::GenSpec* GenSpecJpsiLeak(0);
   if(nGenJpsiLeak>1) GenSpecJpsiLeak = histPdfJpsiLeak->prepareMultiGen(argset2, RooFit::Extended(1), NumEvents(nGenJpsiLeak));

   cout<<"HELLO3"<<endl;


   cout<<"Variable loaded:"<<endl;
   B_plus_M->Print();//B_plus_DTFM_M_zero->Print();
   misPT->Print(); 

   cout<<"HELLO4"<<endl;

   //***************Generate some datasets

   cout<<"Generating signal Zero Photon"<<endl;
   RooDataSet* dataGenSignalZeroGamma = histPdfSignalZeroGamma->generate(*GenSpecSignalZeroGamma);//(argset, 250, false, true, "", false, true);
   dataGenSignalZeroGamma->SetName("dataGenSignalZeroGamma"); dataGenSignalZeroGamma->SetTitle("dataGenSignalZeroGamma");
   cout<<"Generating signal One Photon"<<endl;
   RooDataSet* dataGenSignalOneGamma = histPdfSignalOneGamma->generate(*GenSpecSignalOneGamma);//(argset, 250, false, true, "", false, true);
   dataGenSignalOneGamma->SetName("dataGenSignalOneGamma"); dataGenSignalOneGamma->SetTitle("dataGenSignalOneGamma");
   cout<<"Generating signal two Photons"<<endl;
   RooDataSet* dataGenSignalTwoGamma = histPdfSignalTwoGamma->generate(*GenSpecSignalTwoGamma);//(argset, 250, false, true, "", false, true);
   dataGenSignalTwoGamma->SetName("dataGenSignalTwoGamma"); dataGenSignalTwoGamma->SetTitle("dataGenSignalTwoGamma");
   cout<<"Generating combinatorial"<<endl;
   RooDataSet* dataGenComb = combPDF->generate(*GenSpecComb);//(argset, 100, false, true, "", false, true);
   dataGenComb->SetName("dataGenComb"); dataGenComb->SetTitle("dataGenComb");
   cout<<"Generating PartReco"<<endl;
   RooDataSet* dataGenPartReco = histPdfPartReco->generate(*GenSpecPartReco);//argset, 160, false, true, "", false, true);
   dataGenPartReco->SetName("dataGenPartReco"); dataGenPartReco->SetTitle("dataGenPartReco");
   RooDataSet* dataGenJpsiLeak(0);
   cout<<"HELLO5"<<endl;
   if(nGenJpsiLeak>1)
   {
      cout<<"Generating Leaking JPsi"<<endl;
      dataGenJpsiLeak = histPdfJpsiLeak->generate(*GenSpecJpsiLeak);//argset, 160, false, true, "", false, true);
      dataGenJpsiLeak->SetName("dataGenJpsiLeak"); dataGenJpsiLeak->SetTitle("dataGenJpsiLeak");
   }

   // cout<<"Hello7: "<<dataSetComb<<" "<<dataGenComb<<" "<<combPDF<<endl;
   // if(wantPlots) PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cCombKeys", *B_plus_M, *misPT);
   if(wantPlots) PlotShape(*dataGenComb, *dataGenComb, *combPDF, plotsfile, "cCombKeys", *B_plus_M, *misPT);

   //*************Saving the generated datasets in a workspace

   cout<<"HELLO6"<<endl;
   RooWorkspace* workspaceGen = new RooWorkspace("workspaceGen", "workspaceGen");

   workspaceGen->import(*dataGenSignalZeroGamma);
   workspaceGen->import(*dataGenSignalOneGamma);
   workspaceGen->import(*dataGenSignalTwoGamma);
   workspaceGen->import(*dataGenComb);
   workspaceGen->import(*dataGenPartReco);
   if(nGenJpsiLeak>1) workspaceGen->import(*dataGenJpsiLeak);

   // RooWorkspace workspaceGen("workspaceGen", "workspaceGen");

   // workspaceGen.import(*dataGenSignalZeroGamma);
   // workspaceGen.import(*dataGenSignalOneGamma);
   // workspaceGen.import(*dataGenSignalTwoGamma);
   // workspaceGen.import(*dataGenComb);
   // workspaceGen.import(*dataGenPartReco);
   // if(nGenJpsiLeak>1) workspaceGen.import(*dataGenJpsiLeak);
    
   
   fw->cd();
   workspaceGen->Write("", TObject::kOverwrite);


   delete workspaceGen;
   fw->Close();   
   delete workspace;
   delete fw;
   
   gROOT->Reset();
   
   delete dataGenSignalZeroGamma;
   delete dataGenSignalOneGamma;
   delete dataGenSignalTwoGamma;
   delete dataGenComb;
   delete dataGenPartReco;
   delete dataGenJpsiLeak;

   delete GenSpecSignalZeroGamma;
   delete GenSpecSignalOneGamma;
   delete GenSpecSignalTwoGamma;
   delete GenSpecComb;
   delete GenSpecPartReco;
   delete combPDF;
   delete GenSpecJpsiLeak;

   cout<<"HELLO7"<<endl;
   // delete histPdfSignalZeroGamma; 
   // delete histPdfSignalOneGamma;
   // delete histPdfSignalTwoGamma;
   // delete histPdfPartReco;
   // delete histPdfComb;
   // delete histPdfJpsiLeak;






}




void FitterUtilsHistFact::fit(bool wantplot, bool constPartReco,
      double fracPartReco_const,
      ofstream& out, TTree* t, bool update, string plotsfile)
{


  //***************** Define parameters and observables

   // TFile fw(workspacename.c_str(), "UPDATE");   
   TFile *fw = new TFile(workspacename.c_str(), "UPDATE");   
   RooWorkspace* workspace = (RooWorkspace*)fw->Get("combined");   
   // RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   // RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   // RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   // RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
   // RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   // RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);
   // RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(nPartReco,nSignal));
   // RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));
   RooRealVar* nSignal = workspace->var("nSignal");
   RooRealVar* nPartReco = workspace->var("nPartReco");
   RooRealVar* nComb = workspace->var("nComb");
   RooRealVar* fracZero = workspace->var("fracZero");
   RooRealVar* fracOne = workspace->var("fracOne");
   RooRealVar* fracTwo = workspace->var("fracTwo");
   RooRealVar* nJpsiLeak = workspace->var("nJpsiLeak");
   RooRealVar *l1KeeGen = workspace->var("l1KeeGen");
   RooRealVar *l2KeeGen = workspace->var("l2KeeGen");
   RooRealVar *l3KeeGen = workspace->var("l3KeeGen");
   RooRealVar *l4KeeGen = workspace->var("l4KeeGen");
   RooRealVar *l5KeeGen = workspace->var("l5KeeGen");
   RooRealVar *l1Kee = workspace->var("l1Kee");
   RooRealVar *l2Kee = workspace->var("l2Kee");
   RooRealVar *l3Kee = workspace->var("l3Kee");
   RooRealVar *l4Kee = workspace->var("l4Kee");
   RooRealVar *l5Kee = workspace->var("l5Kee");
   RooRealVar *l1KeeSyst = workspace->var("l1KeeSyst");
   RooRealVar *l2KeeSyst = workspace->var("l2KeeSyst");
   RooRealVar *l3KeeSyst = workspace->var("l3KeeSyst");
   RooRealVar *l4KeeSyst = workspace->var("l4KeeSyst");
   RooRealVar *l5KeeSyst = workspace->var("l5KeeSyst");
   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");
   RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(*nPartReco,*nSignal));

   cout<<"VALUE OF l1Kee IN FIT: "<<l1Kee->getVal()<<" +- "<<l1Kee->getError()<<endl;
   cout<<"VALUE OF l2Kee IN FIT: "<<l2Kee->getVal()<<" +- "<<l2Kee->getError()<<endl;
   cout<<"VALUE OF l3Kee IN FIT: "<<l3Kee->getVal()<<" +- "<<l3Kee->getError()<<endl;
   cout<<"VALUE OF l4Kee IN FIT: "<<l4Kee->getVal()<<" +- "<<l4Kee->getError()<<endl;
   cout<<"VALUE OF l5Kee IN FIT: "<<l5Kee->getVal()<<" +- "<<l5Kee->getError()<<endl;

   // Get observables
   ModelConfig *mc = (ModelConfig*) workspace->obj("ModelConfig");
   RooArgSet *obs = (RooArgSet*) mc->GetObservables();
   RooRealVar *B_plus_M = (RooRealVar*) obs->find(("obs_x_"+channelName).c_str());
   RooRealVar *misPT = (RooRealVar*) obs->find(("obs_y_"+channelName).c_str());
   // B_plus_M->SetNameTitle("B_plus_M", "B_plus_M");
   B_plus_M->SetTitle("B_plus_M");
   B_plus_M->setUnit("MeV/c^{2}");
   // misPT->SetNameTitle("misPT", "misPT");
   misPT->SetTitle("misPT");
   misPT->setUnit("MeV/c");

   cout<<"B_plus_M nbins: "<<B_plus_M->getBins("default")<<endl;
   cout<<"misPT nbins: "<<misPT->getBins("default")<<endl;
   
  // //***********Set Binning
  
  // RooBinning defaultMBins(floor(16), B_plus_M->getMin(), B_plus_M->getMax() ); 
  // RooBinning defaultMisPTBins(floor(40), misPT->getMin(), misPT->getMax()); 

  // B_plus_M->setBinning( defaultMBins);
  // misPT->setBinning( defaultMisPTBins );


  //***************** Build model
   
   // Get HistFactory pdf (no combinatorial)
   RooSimultaneous *model_nocomb = (RooSimultaneous*)mc->GetPdf();
   RooArgList obsTerms;
   RooArgList pois_constraints;
   FactorizeHistFactoryPdf(*obs, *model_nocomb, obsTerms, pois_constraints);

   RooRealSumPdf* model_nocomb_pdf = ( RooRealSumPdf*) obsTerms.first();
   
   RooArgList *model_nocomb_pdf_coefList = (RooArgList*) model_nocomb_pdf->coefList().Clone();
   RooArgList *model_nocomb_pdf_funcList = (RooArgList*) model_nocomb_pdf->funcList().Clone();


   // Define combinatorial model
   RooBinning m_binning(nMassBins,B_plus_M->getMin(),B_plus_M->getMax(),"m_binning");
   RooBinning pT_binning(nmisPTBins,misPT->getMin(),misPT->getMax(),"pT_binning");
   RooExpOfPolyTimesXBinned *combPDF_unnorm =  new RooExpOfPolyTimesXBinned("combPDF_unnorm", "combPDF_unnorm", 
                                                                            *B_plus_M, *misPT, *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee,
                                                                            m_binning, pT_binning);

   RooAbsReal *combIntegral = combPDF_unnorm->createIntegral(RooArgSet(*B_plus_M, *misPT));
   combIntegral->SetNameTitle("combIntegral","combIntegral");
   // RooFormulaVar *combNorm = new RooFormulaVar("combNorm","1./combIntegral",RooArgList(*combIntegral));
   RooRealVar *combNorm = new RooRealVar("combNorm","combNorm",1./combIntegral->getVal());
   RooProduct *combPDF = new RooProduct("combPDF","combPDF",RooArgList(*combPDF_unnorm,*combNorm));  

   model_nocomb_pdf_coefList->add(*nComb);
   model_nocomb_pdf_funcList->add(*combPDF); //normalization

   // Combined model
   RooRealSumPdf *model_pdf = new RooRealSumPdf((channelName+"_comb").c_str(),(channelName+"_comb").c_str(),
                                                *model_nocomb_pdf_funcList,*model_nocomb_pdf_coefList,1);  
   model_pdf->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator") ;
   model_pdf->specialIntegratorConfig(kTRUE)->method2D().setLabel("RooBinIntegrator") ;
   model_pdf->specialIntegratorConfig(kTRUE)->methodND().setLabel("RooBinIntegrator") ;
   model_pdf->forceNumInt();
   model_pdf->setAttribute("GenerateBinned"); // for use with RooSimultaneous::generate in mixed mode

   // Add constraints
   RooArgList *model_and_constraints = (RooArgList*)pois_constraints.Clone();
   model_and_constraints->add(*model_pdf);
   RooProdPdf *model_pdf_constrain = new RooProdPdf("RK_kinematic_comb_constrain","RK_kinematic_comb_constrain",*model_and_constraints);

   // Total model
   RooCategory *idx = (RooCategory*) obs->find("channelCat");
   RooSimultaneous model("modified_pdf","modified_pdf",*idx);
   model.addPdf(*model_pdf_constrain,idx->getLabel());  
   // model.addPdf(*model_nocomb_pdf,idx->getLabel());  

   HistFactorySimultaneous totPdf( model );
   totPdf.SetNameTitle("totPdf", "totPdf");
   

   //***************** Get data

   cout<<"CACA2"<<endl;
   RooWorkspace* workspaceGen = (RooWorkspace*)fw->Get("workspaceGen");
   cout<<"CACA2a"<<endl;
   RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGamma");
   cout<<"CACA2b"<<endl;
   RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen->data("dataGenSignalOneGamma");
   RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGamma");
   RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen->data("dataGenPartReco");
   cout<<"CACA2c"<<endl;
   RooDataSet* dataGenComb = (RooDataSet*)workspaceGen->data("dataGenComb");
   cout<<"CACA2d"<<endl;
   RooDataSet* dataGenJpsiLeak(0);
   if(nGenJpsiLeak>0) dataGenJpsiLeak = (RooDataSet*)workspaceGen->data("dataGenJpsiLeak");


   cout<<"CACA3"<<endl;


   // if(wantplot)
   // {
   //    //**************Must get the datasets

   //    RooDataSet* dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
   //    RooDataSet* dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
   //    RooDataSet* dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
   //    RooDataSet* dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
   //    RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
   //    RooDataSet* dataSetJpsiLeak = (RooDataSet*)workspace->data("dataSetJpsiLeak");

   //    //**************Plot all the different components
   // cout<<"CACA4"<<endl;

   //    cout<<"dataGenSignalZeroGamma: "<<dataGenSignalZeroGamma<<endl;
   //    PlotShape(*dataSetSignalZeroGamma, *dataGenSignalZeroGamma, *histPdfSignalZeroGamma, plotsfile, "cSignalZeroGamma", *B_plus_M, *misPT);
   //    PlotShape(*dataSetSignalOneGamma, *dataGenSignalOneGamma, *histPdfSignalOneGamma, plotsfile, "cSignalOneGamma", *B_plus_M, *misPT);
   //    PlotShape(*dataSetSignalTwoGamma, *dataGenSignalTwoGamma, *histPdfSignalTwoGamma, plotsfile, "cSignalTwoGamma", *B_plus_M, *misPT);
   //    PlotShape(*dataSetPartReco, *dataGenPartReco, *histPdfPartReco, plotsfile, "cPartReco", *B_plus_M, *misPT);
   //    PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cComb", *B_plus_M, *misPT);
   //    PlotShape(*dataSetComb, *dataGenComb, *histPdfComb, plotsfile, "cCombTemplate", *B_plus_M, *misPT);

   // cout<<"CACA5"<<endl;
   //    if(nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeak, *dataGenJpsiLeak, *histPdfJpsiLeak, plotsfile, "cJpsiLeak", *B_plus_M, *misPT);
   // }

   //***************Merge datasets

   cout<<"CACA6"<<endl;
   RooDataSet* dataGenTot_chan(dataGenPartReco);
   dataGenTot_chan->append(*dataGenSignalZeroGamma);
   dataGenTot_chan->append(*dataGenSignalOneGamma);
   dataGenTot_chan->append(*dataGenSignalTwoGamma);
   dataGenTot_chan->append(*dataGenComb);
   if(nGenJpsiLeak>0) dataGenTot_chan->append(*dataGenJpsiLeak);

   RooDataHist *dataGenTot_chan_bin = dataGenTot_chan->binnedClone();
   RooDataHist *dataGenTot = new RooDataHist("dataGenTot","dataGenTot",RooArgSet(*B_plus_M, *misPT),
                                             Index(*idx),Import(idx->getLabel(),*dataGenTot_chan_bin));

   cout<<"Binned data number of bins: "<<dataGenTot->numEntries()<<endl;
   

   //**************Setup constant parameters
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("fracZero")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("fracOne")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("fracTwo")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_sigzero")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_sigone")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_sigtwo")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_partreco")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_jpsileak")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);
  // l1Kee->setConstant(kTRUE);
  l2Kee->setConstant(kTRUE);
  l3Kee->setConstant(kTRUE);
  l4Kee->setConstant(kTRUE);
  l5Kee->setConstant(kTRUE);
  
   cout<<"CACA7"<<endl;
  

   //**************** Constrain the fraction of zero and one photon


   int nGenSignalZeroGamma(floor(nGenFracZeroGamma*nGenSignal));
   int nGenSignalOneGamma(floor(nGenFracOneGamma*nGenSignal));
   int nGenSignalTwoGamma(floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));

   RooRealVar fracZeroConstMean("fracZeroConstMean", "fracZeroConstMean", nGenSignalZeroGamma*1./nGenSignal);
   RooRealVar fracZeroConstSigma("fracZeroConstSigma", "fracZeroConstSigma", sqrt(nGenSignalZeroGamma)/nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", *fracZero, fracZeroConstMean, fracZeroConstSigma); 

   RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", nGenSignalOneGamma*1./nGenSignal/(1-fracZeroConstMean.getVal()));
   RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", sqrt(nGenSignalOneGamma)/nGenSignal/(1-fracZeroConstMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", *fracOne, fracOneConstMean, fracOneConstSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", nGenPartReco/(1.*nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());

   RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

   RooRealVar JpsiLeakMean("JpsiLeakMean", "JpsiLeakMean", nGenJpsiLeak);
   RooRealVar JpsiLeakSigma("JpsiLeakSigma", "JpsiLeakSigma", nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
   RooGaussian JpsiLeakConst("JpsiLeakConst", "JpsiLeakConst", *nJpsiLeak, JpsiLeakMean, JpsiLeakSigma); 
   
   
   cout<<"CACA10"<<endl;
   //Extra TEST CONSTRAINT


   //RooRealVar combConstMean("combConstMean", "combConstMean", nGenComb);
   //RooRealVar combConstSigma("combConstSigma", "combConstSigma", 7.7);
   //RooGaussian combConst("combConst", "combConst", nComb, combConstMean, combConstSigma);

   //**************** fit
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;

   RooArgSet *par_set = totPdf.getParameters(dataGenTot);

   cout<<"CACA11"<<endl;

   initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
          *nSignal, *nPartReco, *nComb, *fracZero, *fracOne, *nJpsiLeak,  constPartReco, fracPartRecoSigma,
          *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);

   RooArgSet constraints(fracZeroConst, fracOneConst);
   if (constPartReco) constraints.add(fracPartRecoConst);
   if(nGenJpsiLeak>0) constraints.add(JpsiLeakConst);

   RooArgList* gammas = new RooArgList();
   ParamHistFunc* param_func=NULL;
   // bool hasStatUncert = getStatUncertaintyFromChannel(model_nocomb_pdf , param_func, gammas );
   bool hasStatUncert = getStatUncertaintyFromChannel(model_pdf_constrain , param_func, gammas );

   cout<< "Param_func nbins: "<<param_func->numBins()<<endl;
   

   RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Offset(kTRUE), Verbose(kTRUE));
   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, ExternalConstraints(constraints),Offset(kTRUE));
   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(constraints));
   RooMinuit minuit(*nll);
   minuit.setStrategy(2);
   minuit.setErrorLevel(0.5);
   minuit.setPrintLevel(1);   

   // std::cout << "Minimizing the Minuit (Migrad)" << std::endl;
   // w->saveSnapshot("TMCPARS",*allpars,kTRUE);
   // minuit_hf->setStrategy(2);
   // minuit_hf->fit("smh");
   // RooFitResult *tempResult=minuit_hf->save("TempResult","TempResult");
   // cout << tempResult->edm() << endl;
   // if (useMinos) minuit_hf->minos(RooArgSet(*poi));
   // result = minuit_hf->save("Result","Result");
 
   int migradRes(1);
   int hesseRes(4);

   vector<int> migradResVec;
   vector<int> hesseResVec;

   double edm(10);
   int nrefit(0);

   RooFitResult* fitRes(0);
   vector<RooFitResult*> fitResVec;

   bool hasConverged(false);


   for(int i(0); (i<15) && !hasConverged ; ++i)
   {
      initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
          *nSignal, *nPartReco, *nComb, *fracZero, *fracOne, *nJpsiLeak,  constPartReco, fracPartRecoSigma,
          *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);

      cout<<"FITTING: starting with nsignal = "<<nSignal->getValV()<<" refit nbr. "<<i<<endl;
      //if(fitRes != NULL && fitRes != 0) delete fitRes;

      // cout<<"Eval nll: "<<nll->getVal()<<endl;
      
      

      migradRes = minuit.migrad();
      hesseRes = minuit.hesse();

      fitRes = minuit.save();
      edm = fitRes->edm();

      fitResVec.push_back(fitRes); 
      migradResVec.push_back(migradRes);
      hesseResVec.push_back(hesseRes);

      if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;

      ++nrefit;


      cout<<"Fitting nbr "<<i<<" done. Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   if(!hasConverged)
   {
      double minNll(1e20);
      int minIndex(-1);
      for(unsigned int i(0); i<fitResVec.size(); ++i)
      {
         if( fitResVec.at(i)->minNll() < minNll)
         {
            minIndex = i;
            minNll = fitResVec[i]->minNll();
         }
      }
      
      migradRes = migradResVec.at(minIndex);
      hesseRes = hesseResVec.at(minIndex);
      cout<<"Fit not converged, choose fit "<<minIndex<<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(t, fitRes,  update, migradRes, hesseRes, hasConverged);
   
   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);
   //totPdf.fitTo(*dataGenTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(wantplot) plot_fit_result(plotsfile, totPdf, dataGenTot);


   // this should remain commented out
   // {
     
   //   //**************Prepare TFile to save the plots
     
   //   TFile f2(plotsfile.c_str(), "UPDATE");
   //   //**************Plot the results of the fit
     
   //   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   //   TIterator *iter = var_set->createIterator();
   //   RooRealVar *var;
     
   //   std::vector<RooPlot*> plots;
   //   RooPlot* frame;
     
     
     
   //   while((var = (RooRealVar*) iter->Next()))
   //   {
       
   //     frame = var->frame();
   //     dataGenTot->plotOn(frame,Name("Data"),DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.6),DrawOption("ZP"));
       
   //     cout<<"here bins"<<endl;
       
   //     totPdf.plotOn(frame, Components("histPdfPartReco"), LineColor(kBlue));
   //     totPdf.plotOn(frame, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
   //     totPdf.plotOn(frame, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
   //     totPdf.plotOn(frame, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
   //     totPdf.plotOn(frame, Components("histPdfJpsiLeak"), LineColor(14));
   //     totPdf.plotOn(frame, Components("combPDF"), LineColor(kBlack));
   //     totPdf.plotOn(frame, LineColor(kRed));
       
   //     plots.push_back(frame);
       
   //   }  
     
   //   if (!(plots.size())) return;
     
   //   TCanvas cFit("cFit", "cFit", 600, 800);
   //   cFit.Divide(1,2);
   //   cFit.cd(1); plots[0]->Draw();
   //   if (plots.size()>1){ 
   //     cFit.cd(2); plots[1]->Draw();
   //   }
     
   //   cFit.Write();
   //   f2.Close();
     
     
   // }
   
   fw->Close();
   delete fw;
   //delete and return
   delete nll;
   delete par_set;
   delete workspace;
   delete workspaceGen;
   cout<<"deleting binned combinatorial"<<endl;
   delete combPDF_unnorm;
   cout<<"deleted!"<<endl;
   delete combPDF;
   delete gammas;
   delete dataGenTot;
   delete model_pdf;
   delete model_pdf_constrain;
   delete combNorm;
   
   

}


void FitterUtilsHistFact::plot_fit_result(string plotsfile, RooSimultaneous &totPdf, RooDataHist *dataGenTot)
{

   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   RooCategory* idx = (RooCategory*) (&totPdf.indexCat());
   
   while((var = (RooRealVar*) iter->Next()))
   {

     if (strcmp(var->ClassName(),"RooCategory")==0) continue;
     
      frame = var->frame();
      dataGenTot->plotOn(frame,Name("Data"),DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.6),DrawOption("ZP"));


      totPdf.plotOn(frame,Name("model"), Slice(*idx),ProjWData(*idx,*dataGenTot),DrawOption("F"),FillColor(kRed),
                             LineWidth(0));
      // resids[i]=frame->pullHist();
      totPdf.plotOn(frame,Name("comb"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                             DrawOption("F"),FillColor(kCyan),Components("combPDF*,*Zero*,*One*,*Two*,*PartReco*,*Leak*"),LineWidth(0));
      dataGenTot->plotOn(frame,Name("Data"),DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.8),DrawOption("ZP"));

      totPdf.plotOn(frame,Name("partreco"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                             DrawOption("F"),FillColor(kMagenta),Components("*Zero*,*One*,*Two*,*PartReco*,*Leak*"),LineWidth(0));

      totPdf.plotOn(frame,Name("partreco"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                             DrawOption("F"),FillColor(kYellow),Components("*Zero*,*One*,*Two*,*Leak*"),LineWidth(0));
      
      totPdf.plotOn(frame,Name("signal"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                             DrawOption("F"),FillColor(kBlue),Components("*Zero*,*One*,*Two*"),LineWidth(0));
      totPdf.plotOn(frame,Name("sigzero"), Slice(*idx),ProjWData(*idx,*dataGenTot),Components("*Zero*"),
                             LineWidth(1),LineColor(kYellow));
      totPdf.plotOn(frame,Name("sigone"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                             LineColor(kViolet),Components("*One*"),LineWidth(1));
      totPdf.plotOn(frame,Name("sigtwo"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                             LineColor(kGreen),Components("*Two*"),LineWidth(1));
      
      
      cout<<"here bins"<<endl;



      // totPdf.plotOn(frame, Components("histPdfPartReco"), LineColor(kBlue));
      // totPdf.plotOn(frame, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
      // totPdf.plotOn(frame, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
      // totPdf.plotOn(frame, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
      // totPdf.plotOn(frame, Components("histPdfJpsiLeak"), LineColor(14));
      // totPdf.plotOn(frame, Components("combPDF"), LineColor(kBlack));
      // totPdf.plotOn(frame, LineColor(kRed));

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

void FitterUtilsHistFact::PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
{
   PlotShape2D(originDataSet, genDataSet, shape, plotsfile, canvName, B_plus_M, misPT);
}

void FitterUtilsHistFact::PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
{
   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");

   //**************Plot Signal Zero Gamma

   TH2F* th2fKey = (TH2F*)shape.createHistogram("th2Shape", B_plus_M, Binning(20), YVar(misPT, Binning(20)));
   cout<<genDataSet.sumEntries()<<endl;
   TH2F* th2fGen = (TH2F*)genDataSet.createHistogram("th2fGen", B_plus_M, Binning(20), YVar(misPT, Binning(20)));

   RooPlot* plotM = B_plus_M.frame();
   originDataSet.plotOn(plotM);
   shape.plotOn(plotM);

   RooPlot* plotMisPT = misPT.frame();
   originDataSet.plotOn(plotMisPT);
   shape.plotOn(plotMisPT);

   TCanvas canv(canvName.c_str(), canvName.c_str(), 800, 800);
   canv.Divide(2,2);
   canv.cd(1); th2fGen->Draw("lego");
   canv.cd(2); th2fKey->Draw("surf");
   canv.cd(3); plotM->Draw();
   canv.cd(4); plotMisPT->Draw();

   canv.Write();

   f2.Close();
}

