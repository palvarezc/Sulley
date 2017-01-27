#include "fitter_utils_ExpOfPolyTimesX.h"
#include "RooBinning.h"
#include "RooRandom.h"
#include "RooExpOfPolyTimesX.h"
#include "RooExpOfPoly.h"

FitterUtilsExpOfPolyTimesX::FitterUtilsExpOfPolyTimesX(int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_, double nGenFracOneGamma_,
    string workspacename_)
   :nGenSignal(nGenSignal_), nGenPartReco(nGenPartReco_), nGenComb(nGenComb_), nGenJpsiLeak(nGenJpsiLeak_), nGenFracZeroGamma(nGenFracZeroGamma_), nGenFracOneGamma(nGenFracOneGamma_), 
   workspacename(workspacename_)
   {}



void FitterUtilsExpOfPolyTimesX::initiateParams(RooArgSet* parset)
{
   RooRealVar *var;
   TIterator *iter = parset->createIterator();

   while((var = (RooRealVar*) iter->Next()))
   {
      if ( !var->isConstant() ) var->randomize();
   }
}


void FitterUtilsExpOfPolyTimesX::initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, RooRealVar& nSignal, RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma, 
      RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee, 
                                                RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen , bool constFracs, bool constComb)
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

   if (!constFracs)
   {
     
   fracZero.setVal(rand.Gaus(fracGenZero, sqrt(nGenSignalZeroGamma)/(1.*nGenSignal))) ;
   fracZero.setRange(0., 1.);
   fracOne.setVal(rand.Gaus(fracGenOne, sqrt(nGenSignalOneGamma)/(1.*nGenSignal))) ;
   fracOne.setRange(0., 1.);
   }


   if(!constComb)
   {
     
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
}


void FitterUtilsExpOfPolyTimesX::prepare_PDFs(string trigStr, string weightStr, string BDTVar, double BDTcut,
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
   RooRealVar weightVar(weightStr.c_str(), weightStr.c_str(), -10, 10);
   RooRealVar BDTRooRealVar(BDTVar.c_str(), BDTVar.c_str(), -1,1);
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", minBMass, maxBMass, "MeV");
   RooRealVar misPT("misPT", "p_{#perp}", 0, 5000, "MeV");
   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV"); 
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);
   RooRealVar K_Kst_isMuon("K_Kst_isMuon", "K_Kst_isMuon", 0, 1);

   RooRealVar weightPartReco("weightPartReco", "weightPartReco", 0, 10);
   RooRealVar weightLeakage("weightLeakage", "weightLeakage", 0, 10);
   RooRealVar dataMCWeightee("DataMCWeightee", "DataMCWeightee",0, 30);


   //***********Set only variables needed

   // tSignal->SetBranchStatus("*", 0); tSignal->SetBranchStatus("B_plus_M", 1); tSignal->SetBranchStatus("misPT", 1); tSignal->SetBranchStatus("B_plus_DTFM_M_zero", 1); tSignal->SetBranchStatus(BDTVar.c_str(),1);
   // tSignal->SetBranchStatus("e_plus_BremMultiplicity", 1); tSignal->SetBranchStatus("e_minus_BremMultiplicity", 1); tSignal->SetBranchStatus(trigStr.c_str()); tSignal->SetBranchStatus("DataMCWeightee",1); tSignal->SetBranchStatus("K_Kst_isMuon", 1);

   // tPartReco->SetBranchStatus("*", 0); tPartReco->SetBranchStatus("B_plus_M", 1); tPartReco->SetBranchStatus("misPT", 1); tPartReco->SetBranchStatus("B_plus_DTFM_M_zero", 1);tPartReco->SetBranchStatus(BDTVar.c_str(),1); tPartReco->SetBranchStatus("K_Kst_isMuon", 1);
   // tPartReco->SetBranchStatus("e_plus_BremMultiplicity", 1); tPartReco->SetBranchStatus("e_minus_BremMultiplicity", 1); tPartReco->SetBranchStatus(trigStr.c_str()); tPartReco->SetBranchStatus("weightPartReco",1);

   // tComb->SetBranchStatus("*", 0); tComb->SetBranchStatus("B_plus_M", 1); tComb->SetBranchStatus("misPT", 1); tComb->SetBranchStatus("B_plus_DTFM_M_zero", 1);tComb->SetBranchStatus(BDTVar.c_str(),1);
   // tComb->SetBranchStatus("e_plus_BremMultiplicity", 1); tComb->SetBranchStatus("e_minus_BremMultiplicity", 1); tComb->SetBranchStatus(trigStr.c_str()); tComb->SetBranchStatus("K_Kst_isMuon");

   // tJpsiLeak->SetBranchStatus("*", 0); tJpsiLeak->SetBranchStatus("B_plus_M", 1); tJpsiLeak->SetBranchStatus("misPT", 1); 
   // tJpsiLeak->SetBranchStatus("B_plus_DTFM_M_zero", 1);tJpsiLeak->SetBranchStatus(BDTVar.c_str(),1);
   // tJpsiLeak->SetBranchStatus("e_plus_BremMultiplicity", 1); tJpsiLeak->SetBranchStatus("e_minus_BremMultiplicity", 1); tJpsiLeak->SetBranchStatus(trigStr.c_str()); 
   // tJpsiLeak->SetBranchStatus("weightLeakage",1); tJpsiLeak->SetBranchStatus("K_Kst_isMuon");



   tSignal->SetBranchStatus("*", 0); tSignal->SetBranchStatus("B_plus_M", 1); tSignal->SetBranchStatus("misPT", 1); tSignal->SetBranchStatus("B_plus_DTFM_M_zero", 1); tSignal->SetBranchStatus(BDTVar.c_str(),1);
   tSignal->SetBranchStatus("e_plus_BremMultiplicity", 1); tSignal->SetBranchStatus("e_minus_BremMultiplicity", 1); tSignal->SetBranchStatus(trigStr.c_str()); tSignal->SetBranchStatus(weightStr.c_str(),1); tSignal->SetBranchStatus("K_Kst_isMuon", 1);

   tPartReco->SetBranchStatus("*", 0); tPartReco->SetBranchStatus("B_plus_M", 1); tPartReco->SetBranchStatus("misPT", 1); tPartReco->SetBranchStatus("B_plus_DTFM_M_zero", 1);tPartReco->SetBranchStatus(BDTVar.c_str(),1); tPartReco->SetBranchStatus("K_Kst_isMuon", 1);
   tPartReco->SetBranchStatus("e_plus_BremMultiplicity", 1); tPartReco->SetBranchStatus("e_minus_BremMultiplicity", 1); 
   tPartReco->SetBranchStatus(trigStr.c_str()); tPartReco->SetBranchStatus(weightStr.c_str(),1);

   tComb->SetBranchStatus("*", 0); tComb->SetBranchStatus("B_plus_M", 1); tComb->SetBranchStatus("misPT", 1); tComb->SetBranchStatus("B_plus_DTFM_M_zero", 1);tComb->SetBranchStatus(BDTVar.c_str(),1);
   tComb->SetBranchStatus("e_plus_BremMultiplicity", 1); tComb->SetBranchStatus("e_minus_BremMultiplicity", 1); tComb->SetBranchStatus(trigStr.c_str()); tComb->SetBranchStatus("K_Kst_isMuon");

   tJpsiLeak->SetBranchStatus("*", 0); tJpsiLeak->SetBranchStatus("B_plus_M", 1); tJpsiLeak->SetBranchStatus("misPT", 1); 
   tJpsiLeak->SetBranchStatus("B_plus_DTFM_M_zero", 1);tJpsiLeak->SetBranchStatus(BDTVar.c_str(),1);
   tJpsiLeak->SetBranchStatus("e_plus_BremMultiplicity", 1); tJpsiLeak->SetBranchStatus("e_minus_BremMultiplicity", 1); 
   tJpsiLeak->SetBranchStatus(trigStr.c_str()); 
   tJpsiLeak->SetBranchStatus(weightStr.c_str(),1); tJpsiLeak->SetBranchStatus("K_Kst_isMuon");


   //***********Set Binning


   RooBinning defaultMBins(floor((maxBMass-minBMass)/(40.)), B_plus_M.getMin(), B_plus_M.getMax() ); 
   RooBinning defaultMisPTBins(floor(40), misPT.getMin(), misPT.getMax()); 
   RooBinning broaderMBins(floor((maxBMass-minBMass)/(80.)), B_plus_M.getMin(), B_plus_M.getMax()); 
   RooBinning broaderMisPTBins(floor(40), misPT.getMin(), misPT.getMax()); 

   B_plus_M.setBinning( defaultMBins);
   misPT.setBinning( defaultMisPTBins );
   B_plus_M.setBinning( broaderMBins, "broaderBins");
   misPT.setBinning( broaderMisPTBins, "broaderBins" );

   B_plus_DTFM_M_zero.setBins(100);

   // RooArgSet argset(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, K_Kst_isMuon);
   // RooArgSet argsetPartReco(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightPartReco);
   // RooArgSet argsetLeakage(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightLeakage);
   // RooArgSet argsetSignal(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, dataMCWeightee);


   RooArgSet argset(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, K_Kst_isMuon);
   RooArgSet argsetPartReco(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightVar);
   RooArgSet argsetLeakage(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightVar);
   RooArgSet argsetSignal(BDTRooRealVar, B_plus_DTFM_M_zero, misPT,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightVar);



   cout<<"getting the datasets:"<<endl;

   RooDataSet* dataSetSignalZeroGamma;
   RooDataSet* dataSetSignalOneGamma;
   RooDataSet* dataSetSignalTwoGamma;
   RooDataSet* dataSetPartReco;
   RooDataSet* dataSetJpsiLeak;
   RooDataSet* dataSetComb;

   TFile* fw(NULL);

   // Uncomment to recover the BDT cut!!!
   // string BDTCutString( ("("+BDTVar+">"+d2s(BDTcut)+")").c_str()  );
   string BDTCutString = "1";

   // dataSetSignalZeroGamma = new RooDataSet("dataSetSignalZeroGamma", "dataSetSignalZeroGamma", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("DataMCWeightee")  );

   // dataSetSignalOneGamma = new RooDataSet("dataSetSignalOneGamma", "dataSetSignalOneGamma", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("DataMCWeightee")  );

   // dataSetSignalTwoGamma = new RooDataSet("dataSetSignalTwoGamma", "dataSetSignalTwoGamma", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("DataMCWeightee")  );


   // dataSetPartReco = new RooDataSet("dataSetPartReco", "dataSetPartReco",  argsetPartReco, Import(*tPartReco),Cut(("("+trigStr+"  > 0.9) && "+BDTCutString+ " && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("weightPartReco"));

   // dataSetJpsiLeak = new RooDataSet("dataSetJpsiLeak", "dataSetJpsiLeak",  argsetLeakage, Import(*tJpsiLeak),Cut(("B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("weightLeakage"));

   //  // dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("("+trigStr+"  > 0.9) && (UBDT3R > "+d2s(BDTcut-0.03)+")  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   // dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("K_Kst_isMuon == 0 && ("+trigStr+"  > 0.9) && "+BDTCutString+"  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());

   // misPT.setRange(-2000, 5000);
   // RooDataSet* dataSetCombExt = new RooDataSet("dataSetCombExt", "dataSetCombExt", tComb, argset, ("K_Kst_isMuon == 0 && ("+trigStr+"  > 0.9) && "+BDTCutString+"  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   // misPT.setRange(0, 5000);



   dataSetSignalZeroGamma = new RooDataSet("dataSetSignalZeroGamma", "dataSetSignalZeroGamma", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar(weightStr.c_str())  );

   dataSetSignalOneGamma = new RooDataSet("dataSetSignalOneGamma", "dataSetSignalOneGamma", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar(weightStr.c_str())  );

   dataSetSignalTwoGamma = new RooDataSet("dataSetSignalTwoGamma", "dataSetSignalTwoGamma", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar(weightStr.c_str())  );


   dataSetPartReco = new RooDataSet("dataSetPartReco", "dataSetPartReco",  argsetPartReco, Import(*tPartReco),Cut(("("+trigStr+"  > 0.9) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar("weightPartReco"));

   dataSetJpsiLeak = new RooDataSet("dataSetJpsiLeak", "dataSetJpsiLeak",  argsetLeakage, Import(*tJpsiLeak),Cut(("B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()), WeightVar(weightStr.c_str()));

    // dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("("+trigStr+"  > 0.9) && (UBDT3R > "+d2s(BDTcut-0.03)+")  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("K_Kst_isMuon == 0 && ("+trigStr+"  > 0.9) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());

   misPT.setRange(-2000, 5000);
   RooDataSet* dataSetCombExt = new RooDataSet("dataSetCombExt", "dataSetCombExt", tComb, argset, ("K_Kst_isMuon == 0 && ("+trigStr+"  > 0.9) && "+BDTCutString+"  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   misPT.setRange(0, 5000);


   cout<<"Number of zero: "<< dataSetSignalZeroGamma->sumEntries()<<endl;
   cout<<"Number of one: "<< dataSetSignalOneGamma->sumEntries()<<endl;
   cout<<"Number of two: "<< dataSetSignalTwoGamma->sumEntries()<<endl;
   cout<<"Number of PartReco: "<< dataSetPartReco->sumEntries()<<endl;
   cout<<"Number of Jpsi leaking:"<< dataSetJpsiLeak->sumEntries()<<endl;
   cout<<"Number of combinatorial events:"<< dataSetComb->sumEntries()<<endl;


   cout<<"binning the datasets:"<<endl;

   RooArgSet argset2(B_plus_M, misPT);

   RooDataHist dataHistSignalZeroGamma("dataHistSignalZeroGamma", "dataHistSignalZeroGamma", argset2, *dataSetSignalZeroGamma); 
   RooDataHist dataHistSignalOneGamma("dataHistSignalOneGamma", "dataHistSignalOneGamma", argset2, *dataSetSignalOneGamma); 
   RooDataHist dataHistSignalTwoGamma("dataHistSignalTwoGamma", "dataHistSignalTwoGamma", argset2, *dataSetSignalTwoGamma); 
   RooDataHist dataHistComb("dataHistComb", "dataHistComb", argset2, *dataSetComb); 
   RooDataHist dataHistPartReco("dataHistPartReco", "dataHistPartReco", argset2, *dataSetPartReco); 
   // RooDataHist dataHistJpsiLeak("dataHistJpsiLeak", "dataHistJpsiLeak", argset2, "broaderBins");
   // dataHistJpsiLeak.add(*dataSetJpsiLeak); 
   RooDataHist dataHistJpsiLeak("dataHistJpsiLeak", "dataHistJpsiLeak", argset2, *dataSetJpsiLeak); 

   //*************** Compute Error on J/psi leak

   // double ErrorJpsi(0);
   // if(dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && "+BDTCutString+" && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()) > 0) ErrorJpsi = 1./sqrt(dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && "+BDTCutString+" && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()));
   // RooRealVar fractionalErrorJpsiLeak("fractionalErrorJpsiLeak", "fractionalErrorJpsiLeak", ErrorJpsi);
   // cout<<"JPSI LEAK: "<<dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && "+BDTCutString+" && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   // cout<<"JPSI LEAK fractional Error: "<<ErrorJpsi<<endl;


   double ErrorJpsi(0);
   if(dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()) > 0) ErrorJpsi = 1./sqrt(dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()));
   RooRealVar fractionalErrorJpsiLeak("fractionalErrorJpsiLeak", "fractionalErrorJpsiLeak", ErrorJpsi);
   cout<<"JPSI LEAK: "<<dataSetJpsiLeak->sumEntries(("("+trigStr+"  > 0.9) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());
   cout<<"JPSI LEAK fractional Error: "<<ErrorJpsi<<endl;


   //***************Create 2D histogram estimates from data


   cout<<"Preparing the 3 2D histPdf: 1";
   //   RooArgSet argset2(B_plus_M);

   // RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma,2); cout<<" 2";
   // RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma,2); cout<<" 3";
   // RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma,2); cout<<" 4";
   // RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco,2); cout<<" 5";
   // RooHistPdf histPdfJpsiLeak("histPdfJpsiLeak", "histPdfJpsiLeak", argset2, dataHistJpsiLeak,2); cout<<" 6";

   // No interpolation
   RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma); cout<<" 2";
   RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma); cout<<" 3";
   RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma); cout<<" 4";
   RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco); cout<<" 5";
   RooHistPdf histPdfJpsiLeak("histPdfJpsiLeak", "histPdfJpsiLeak", argset2, dataHistJpsiLeak); cout<<" 6";



   //***************Create combinatorial from fit to data

   RooRealVar l1Kee("l1Kee", "l1Kee", +2.3e-3, -1e-1, 1e-1);
   RooRealVar l2Kee("l2Kee", "l2Kee", +4.9e-3, -1e-1, 1e-1);
   RooRealVar l3Kee("l3Kee", "l3Kee", -4.1e-7, -1e-5, 1e-5);
   RooRealVar l4Kee("l4Kee", "l4Kee", -4e-7, -1e-5, 1e-5);
   RooRealVar l5Kee("l5Kee", "l5Kee", 1e-10, -1e-9, 1e-9);

   RooAbsPdf *combPDF;
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


   RooHistPdf histPdfComb("histPdfComb", "histPdfComb", argset2, dataHistComb,2);



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
   RooWorkspace workspace("workspace", "workspace");
   workspace.import(B_plus_DTFM_M_zero);
   workspace.import(B_plus_M);
   workspace.import(misPT);
   workspace.import(l1Kee);
   workspace.import(l2Kee);
   workspace.import(l3Kee);
   workspace.import(l4Kee);
   workspace.import(l5Kee);
   workspace.import(l1KeeGen);
   workspace.import(l2KeeGen);
   workspace.import(l3KeeGen);
   workspace.import(l4KeeGen);
   workspace.import(l5KeeGen);
   workspace.import(l1KeeSyst);
   workspace.import(l2KeeSyst);
   workspace.import(l3KeeSyst);
   workspace.import(l4KeeSyst);
   workspace.import(l5KeeSyst);
   workspace.import(*dataSetSignalZeroGamma);
   workspace.import(*dataSetSignalOneGamma);
   workspace.import(*dataSetSignalTwoGamma);
   workspace.import(*dataSetPartReco);
   workspace.import(*dataSetComb);
   workspace.import(*dataSetCombExt);
   workspace.import(*dataSetJpsiLeak);
   workspace.import(histPdfSignalZeroGamma);
   workspace.import(histPdfSignalOneGamma);
   workspace.import(histPdfSignalTwoGamma);
   workspace.import(histPdfPartReco);
   workspace.import(histPdfJpsiLeak);
   workspace.import(histPdfComb);
   workspace.import(fractionalErrorJpsiLeak);
   workspace.writeToFile(workspacename.c_str());



   delete fComb;
   delete fSignal;
   delete fPartReco;
   if(fw!=0 && fw != NULL) delete fw;
   delete combPDF;

   delete dataSetSignalZeroGamma; 
   delete dataSetSignalOneGamma;
   delete dataSetSignalTwoGamma;
   delete dataSetPartReco;
   delete dataSetJpsiLeak;
   delete dataSetComb;
   delete dataSetCombExt;
}


void FitterUtilsExpOfPolyTimesX::generate(bool wantPlots, string plotsfile)
{
   //***************Get the PDFs from the workspace

   TFile fw(workspacename.c_str(), "UPDATE");   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooRealVar *l1KeeGen = workspace->var("l1KeeGen");
   RooRealVar *l2KeeGen = workspace->var("l2KeeGen");
   RooRealVar *l3KeeGen = workspace->var("l3KeeGen");
   RooRealVar *l4KeeGen = workspace->var("l4KeeGen");
   RooRealVar *l5KeeGen = workspace->var("l5KeeGen");
   RooRealVar *l1KeeSyst = workspace->var("l1KeeSyst");
   RooRealVar *l2KeeSyst = workspace->var("l2KeeSyst");
   RooRealVar *l3KeeSyst = workspace->var("l3KeeSyst");
   RooRealVar *l4KeeSyst = workspace->var("l4KeeSyst");
   RooRealVar *l5KeeSyst = workspace->var("l5KeeSyst");
   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");
   RooDataSet* dataSetCombExt = (RooDataSet*)workspace->data("dataSetCombExt");
   RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");

   cout<<"VALUE OF l1 IN GENERATE: "<<l1KeeGen->getVal()<<" +- "<<l1KeeGen->getError()<<endl;
   cout<<"VALUE OF l2 IN GENERATE: "<<l2KeeGen->getVal()<<" +- "<<l2KeeGen->getError()<<endl;
   cout<<"VALUE OF l3 IN GENERATE: "<<l3KeeGen->getVal()<<" +- "<<l3KeeGen->getError()<<endl;
   cout<<"VALUE OF l4 IN GENERATE: "<<l4KeeGen->getVal()<<" +- "<<l4KeeGen->getError()<<endl;
   cout<<"VALUE OF l5 IN GENERATE: "<<l5KeeGen->getVal()<<" +- "<<l5KeeGen->getError()<<endl;


   cout<<"HELLO0"<<endl;
   RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
   RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
   RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
   RooHistPdf *histPdfComb = (RooHistPdf *) workspace->pdf("histPdfComb");
   RooHistPdf *histPdfJpsiLeak(0);
   cout<<"HELLO0a"<<endl;
   if(nGenJpsiLeak>1) histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");


   RooAbsPdf *combPDF;
   combPDF =  new RooExpOfPolyTimesX("combPDF", "combPDF", *B_plus_M, *misPT, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen );

   //WARNING WARNING WARNING WARNING LITTLE FUCKER: USE THESE LINES ONLY IF YOU WANT TO ESTIMATE SYSTEMATICS

//   combPDF = histPdfComb;

//   RooAbsPdf *combPDF;
//   combPDF =  new RooExpOfPoly("combPDF", "combPDF", *B_plus_M, *misPT, *l1KeeSyst, *l2KeeSyst, *l3KeeSyst, *l4KeeSyst, *l5KeeSyst );
//   cout<<"WARNING WARNING WARNING: systematic mode, gen function different from fit function"<<endl;


   ///////////////create a keysPdf also for systematics
   // cout<<"HELLO0b"<<endl;
   // TVectorD rho(2);
   // rho[0] = 2.5;
   // rho[1] = 1.5;
   // cout<<"HELLO0c"<<endl;
   // misPT->setRange(-2000, 5000);
   // cout<<"HELLO0d "<<dataSetCombExt<<endl;
   // combPDF = new RooNDKeysPdf("keysPDFComb", "keysPDFComb", RooArgList(*B_plus_M, *misPT), *dataSetCombExt, rho, "ma",3, true);
   // cout<<"HELLO0e"<<endl;
   // misPT->setRange(0, 5000);

   cout<<"HELLO1"<<endl;

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

   cout<<"Hello7: "<<dataSetComb<<" "<<dataGenComb<<" "<<combPDF<<endl;
   if(wantPlots) PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cCombKeys", *B_plus_M, *misPT);

   //*************Saving the generated datasets in a workspace

   cout<<"HELLO6"<<endl;
   RooWorkspace workspaceGen("workspaceGen", "workspaceGen");

   workspaceGen.import(*dataGenSignalZeroGamma);
   workspaceGen.import(*dataGenSignalOneGamma);
   workspaceGen.import(*dataGenSignalTwoGamma);
   workspaceGen.import(*dataGenComb);
   workspaceGen.import(*dataGenPartReco);
   if(nGenJpsiLeak>1) workspaceGen.import(*dataGenJpsiLeak);
    
   
   fw.cd();
   workspaceGen.Write("", TObject::kOverwrite);


   //delete workspace;
   fw.Close();


   delete dataGenSignalZeroGamma;
   delete dataGenSignalOneGamma;
   delete dataGenSignalTwoGamma;
   delete dataGenComb;
   delete dataGenPartReco;
   if(nGenJpsiLeak>1) delete dataGenJpsiLeak;

   delete GenSpecSignalZeroGamma;
   delete GenSpecSignalOneGamma;
   delete GenSpecSignalTwoGamma;
   delete GenSpecComb;
   delete GenSpecPartReco;
   delete combPDF;
   if(nGenJpsiLeak>1) delete GenSpecJpsiLeak;

   delete histPdfSignalZeroGamma; 
   delete histPdfSignalOneGamma;
   delete histPdfSignalTwoGamma;
   delete histPdfPartReco;
   delete histPdfComb;
   if(nGenJpsiLeak>1) delete histPdfJpsiLeak;

   cout<<"HELLO7"<<endl;

}




void FitterUtilsExpOfPolyTimesX::fit(bool wantplot, bool constPartReco,
      double fracPartReco_const,
      ofstream& out, TTree* t, bool update, string plotsfile)
{


  bool constFracs =0;
  bool constComb =0;
  
   //***************Get the PDFs from the workspace

   TFile fw(workspacename.c_str());   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
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

   cout<<"VALUE OF l1Kee IN FIT: "<<l1Kee->getVal()<<" +- "<<l1Kee->getError()<<endl;
   cout<<"VALUE OF l2Kee IN FIT: "<<l2Kee->getVal()<<" +- "<<l2Kee->getError()<<endl;
   cout<<"VALUE OF l3Kee IN FIT: "<<l3Kee->getVal()<<" +- "<<l3Kee->getError()<<endl;
   cout<<"VALUE OF l4Kee IN FIT: "<<l4Kee->getVal()<<" +- "<<l4Kee->getError()<<endl;
   cout<<"VALUE OF l5Kee IN FIT: "<<l5Kee->getVal()<<" +- "<<l5Kee->getError()<<endl;

   RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
   RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
   RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
   RooHistPdf *histPdfComb = (RooHistPdf *) workspace->pdf("histPdfComb");
   RooHistPdf *histPdfJpsiLeak(0);
   if(nGenJpsiLeak>0) histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");

   RooAbsPdf *combPDF;

   cout<<"CACA1"<<endl;

   //Here, put "Gen" if the parameter remains fixed
   
   //************ONly for systematics
  // l1Kee->setVal(l1KeeSyst->getVal());
  // l2Kee->setVal(l2KeeSyst->getVal());
  // l3Kee->setVal(l3KeeSyst->getVal());
  // l4Kee->setVal(l4KeeSyst->getVal());
  // l5Kee->setVal(l5KeeSyst->getVal());
   //*********************************

   // combPDF =  new RooExpOfPolyTimesX("combPDF", "combPDF", *B_plus_M, *misPT, *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee);
   combPDF =  new RooExpOfPolyTimesX("combPDF", "combPDF", *B_plus_M, *misPT, *l1Kee, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);


   cout<<"CACA2"<<endl;
   RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
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
   cout<<"CACA4"<<endl;

      cout<<"dataGenSignalZeroGamma: "<<dataGenSignalZeroGamma<<endl;
      PlotShape(*dataSetSignalZeroGamma, *dataGenSignalZeroGamma, *histPdfSignalZeroGamma, plotsfile, "cSignalZeroGamma", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalOneGamma, *dataGenSignalOneGamma, *histPdfSignalOneGamma, plotsfile, "cSignalOneGamma", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalTwoGamma, *dataGenSignalTwoGamma, *histPdfSignalTwoGamma, plotsfile, "cSignalTwoGamma", *B_plus_M, *misPT);
      PlotShape(*dataSetPartReco, *dataGenPartReco, *histPdfPartReco, plotsfile, "cPartReco", *B_plus_M, *misPT);
      PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cComb", *B_plus_M, *misPT);
      PlotShape(*dataSetComb, *dataGenComb, *histPdfComb, plotsfile, "cCombTemplate", *B_plus_M, *misPT);

   cout<<"CACA5"<<endl;
      if(nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeak, *dataGenJpsiLeak, *histPdfJpsiLeak, plotsfile, "cJpsiLeak", *B_plus_M, *misPT);
   }

   //***************Merge datasets

   cout<<"CACA6"<<endl;
   RooDataSet* dataGenTot(dataGenPartReco);
   dataGenTot->append(*dataGenSignalZeroGamma);
   dataGenTot->append(*dataGenSignalOneGamma);
   dataGenTot->append(*dataGenSignalTwoGamma);
   dataGenTot->append(*dataGenComb);
   if(nGenJpsiLeak>0) dataGenTot->append(*dataGenJpsiLeak);

   RooDataHist *dataGenTotHist = dataGenTot->binnedClone("dataGenTotHist","dataGenTotHist");
   

   cout<<"CACA7"<<endl;
   //**************Prepare fitting function

   RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
   RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);
   RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(nPartReco,nSignal));
   RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));

   // nJpsiLeak.setConstant(kTRUE);
   

   cout<<"CACA8"<<endl;

   RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", RooArgList(*histPdfSignalZeroGamma, *histPdfSignalOneGamma, *histPdfSignalTwoGamma), RooArgList(fracZero, fracOneRec));

   RooArgList pdfList(histPdfSignal, *histPdfPartReco, *combPDF);
   RooArgList yieldList(nSignal, nPartReco, nComb);

   // RooArgList pdfList(histPdfSignal, *histPdfPartReco);
   // RooArgList yieldList(nSignal, nPartReco);

   if(nGenJpsiLeak>0)
   {
      pdfList.add(*histPdfJpsiLeak);
      yieldList.add(nJpsiLeak); 
   }
   RooAddPdf totPdf("totPdf", "totPdf", pdfList, yieldList);
   cout<<"CACA9"<<endl;

   //**************** Constrain the fraction of zero and one photon


   int nGenSignalZeroGamma(floor(nGenFracZeroGamma*nGenSignal));
   int nGenSignalOneGamma(floor(nGenFracOneGamma*nGenSignal));
   int nGenSignalTwoGamma(floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));


   //   This is the code from Thibaud
   RooRealVar fracZeroConstMean("fracZeroConstMean", "fracZeroConstMean", nGenSignalZeroGamma*1./nGenSignal);
   RooRealVar fracZeroConstSigma("fracZeroConstSigma", "fracZeroConstSigma", sqrt(nGenSignalZeroGamma)/nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroConstMean, fracZeroConstSigma); 

   RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", nGenSignalOneGamma*1./nGenSignal/(1-fracZeroConstMean.getVal()));
   RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", sqrt(nGenSignalOneGamma)/nGenSignal/(1-fracZeroConstMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneConstMean, fracOneConstSigma); 


   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", nGenPartReco/(1.*nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());

   RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

   RooRealVar JpsiLeakMean("JpsiLeakMean", "JpsiLeakMean", nGenJpsiLeak);
   RooRealVar JpsiLeakSigma("JpsiLeakSigma", "JpsiLeakSigma", nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
   RooGaussian JpsiLeakConst("JpsiLeakConst", "JpsiLeakConst", nJpsiLeak, JpsiLeakMean, JpsiLeakSigma); 
   cout<<"CACA10"<<endl;

   //Extra TEST CONSTRAINT   
   


   //RooRealVar combConstMean("combConstMean", "combConstMean", nGenComb);
   //RooRealVar combConstSigma("combConstSigma", "combConstSigma", 7.7);
   //RooGaussian combConst("combConst", "combConst", nComb, combConstMean, combConstSigma);

   //**************** fit
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;

   // RooArgSet *par_set = totPdf.getParameters(dataGenTot);
   RooArgSet *par_set = totPdf.getParameters(dataGenTotHist);

   cout<<"CACA11"<<endl;

   initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
          nSignal, nPartReco, nComb, fracZero, fracOne, nJpsiLeak,  constPartReco, fracPartRecoSigma,
                  *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen, constFracs, constComb);

   RooArgSet constraints(fracZeroConst, fracOneConst);
   if (constPartReco) constraints.add(fracPartRecoConst);
   if(nGenJpsiLeak>0) constraints.add(JpsiLeakConst);

   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(constraints));
   RooAbsReal* nll = totPdf.createNLL(*dataGenTotHist, Extended(), ExternalConstraints(constraints));
   RooMinuit minuit(*nll);
   minuit.setStrategy(2);


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
          nSignal, nPartReco, nComb, fracZero, fracOne, nJpsiLeak,  constPartReco, fracPartRecoSigma,
                     *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen, constFracs, constComb);

      cout<<"FITTING: starting with nsignal = "<<nSignal.getValV()<<" refit nbr. "<<i<<endl;
      //if(fitRes != NULL && fitRes != 0) delete fitRes;

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

   if(wantplot) plot_fit_result(plotsfile, totPdf, *dataGenTotHist);


   fw.Close();
   //delete and return
   delete nll;
   delete par_set;
   delete workspace;
   delete workspaceGen;
   delete combPDF;

}


void FitterUtilsExpOfPolyTimesX::plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataSet dataGenTot)
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



void FitterUtilsExpOfPolyTimesX::plot_fit_result(string plotsfile, RooAbsPdf &totPdf, RooDataHist dataGenTot)
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


void FitterUtilsExpOfPolyTimesX::PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
{
   PlotShape2D(originDataSet, genDataSet, shape, plotsfile, canvName, B_plus_M, misPT);
}

void FitterUtilsExpOfPolyTimesX::PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
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

