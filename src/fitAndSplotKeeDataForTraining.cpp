#include "fitAndSplotKeeDataForTraining.h"
#include "TRandom.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "TCanvas.h"
#include "TTree.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooStats/SPlot.h"
//#include "RooDoubleCrystalBall.hpp"
#include "RooGaussian.h"
#include "RooNumIntConfig.h"
#include "TTreeFormula.h"
#include "RooKeysPdf.h"
// #include "RooHypatia.hpp"
// #include "RooHypatia2.hpp"
// #include "RooIpatia2.hpp"
#include "RooChebychev.h"
#include<boost/algorithm/string/replace.hpp>

void FitAndSplotKeeDataForTraining::initiateHistoYield()
{
   TFile fHistoYield(fHistoYieldName.c_str(), "RECREATE");
   TH1D histoYield("histoYield", "histoYield", 100, 0, 100);
   histoYield.Write("", TObject::kOverwrite);
   fHistoYield.Close();
}

void FitAndSplotKeeDataForTraining::addMassShiftScale(int trigCat, string extracut, string label)
{

  string trigCatS("Trig"+i2s(trigCat));

  //get workspace
  TFile fw(workspaceFileName.c_str());

  RooWorkspace* workspaceFit = (RooWorkspace*)fw.Get(("workspaceDataFit"+trigCatS).c_str());

  workspaceFit->Print();
  if(!workspaceFit)
  {
    cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, no workspace found in file "<<workspaceFileName<<endl;
    return;
  }

  RooRealVar *sigmaScaleFactor = workspaceFit->var(("sigmaScaleFactor"+trigCatS).c_str());
  RooRealVar *meanShift = workspaceFit->var(("meanShift"+trigCatS).c_str());


  //get the trees
  string totFileName(tupleMCdir+"/"+tupleMCname);
  //if(trigCat == 10) totFileName = "/vols/lhcb/th1011/KeeTuples/data/presel/B2Kee_Strip21_data_JpsiKPresel.root";

  TFile f( totFileName.c_str());


  TTree* t = (TTree*)f.Get(treeDataName.c_str());
  if(!t)
  {
    cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, no tree "<<treeDataName<<" found in "<<tupleDataDir<<"/"<<tupleDataName<<endl;
    return;
  }

   double B_plus_M;
   t->SetBranchAddress("B_plus_M", &B_plus_M);

  //prepare the cut to apply to the tree

  //  string trigMassCut( ("(B_plus_DTFM_M_zero>"+d2s(B_plus_DTF_M_cut)+" && B_plus_M>"+d2s(minBMass_data)+" && B_plus_M <"+d2s(maxBMass_data)).c_str() );
  // trigMassCut += " && passTrigCat"+i2s(trigCat)+" > 0.5)";

  string trigMassCut( "passTrigCat"+i2s(trigCat)+" > 0.5");

   TTreeFormula ttfTrigMassCut("ttfTrigMassCut", trigMassCut.c_str(), t);

   //clone the tree

   string nameNewFile(totFileName);
   nameNewFile.insert(nameNewFile.size()-5, "_"+label);
   TFile f2(nameNewFile.c_str(), "RECREATE");
   TTree* t2 = t->CloneTree(0);

   //fill the cloned tree with sWeights
   double sigmascale;
   double meanshift;
   double B_plus_M_new;
   double B_plus_M_NoMassShiftScale;
   double massBplusMC = 5279.15;

   t2->SetBranchAddress("B_plus_M",&B_plus_M_new);

   t2->Branch(("MassShift_"+trigCatS).c_str(), &meanshift, ("MassShift_"+trigCatS+"/D").c_str());
   t2->Branch(("MassScale_"+trigCatS).c_str(), &sigmascale, ("MassScale_Trig"+trigCatS+"/D").c_str());
   t2->Branch("B_plus_M_NoMassShiftScale", &B_plus_M_NoMassShiftScale, "B_plus_M_NoMassShiftScale/D");

   cout<<"Putting the mass shift and scale in the tree... "<<endl;
   
   int nEntries(t->GetEntries());

   if(extracut == "") extracut = "1";
   TTreeFormula ttfExtraCut( "ttfExtraCut", extracut.c_str(), t);

   int j(0);

   for(int i(0); i<nEntries; ++i)
   {
      t->GetEntry(i);

      if(i % (nEntries/10) == 0) cout<<100*i/nEntries<<"\% "<<flush;

      if(ttfTrigMassCut.EvalInstance())
      {
        sigmascale = sigmaScaleFactor->getVal();
        meanshift = meanShift->getVal();        
        B_plus_M_NoMassShiftScale = B_plus_M;
        B_plus_M_new = B_plus_M_NoMassShiftScale*sigmascale + massBplusMC*(1-sigmascale) + meanshift;
        
         ++j;

         if(ttfExtraCut.EvalInstance()) t2->Fill();
      }
   }

   cout<<endl;

   f2.cd();
   t2->Write();

   fw.Close();
   f.Close();
   f2.Close();
}

void FitAndSplotKeeDataForTraining::makeSWeightedTree(int trigCat, string extracut, string label)
{
   string trigCatS("Trig"+i2s(trigCat));

   //Get the sWeighted data

   TFile fw(workspaceFileName.c_str());
   RooWorkspace* workspaceSPlot = (RooWorkspace*)fw.Get(("workspaceSPlot"+trigCatS).c_str());

   if(!workspaceSPlot)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, no workspace found"<<endl; 
      return;
   }
   
   RooDataSet* wDataSet = (RooDataSet*)workspaceSPlot->data(("wDataSet"+trigCatS).c_str());

   if(!wDataSet)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, error downloading stuff from workspace"<<endl;
      return;
   }

   //get the trees


   string totFileName(tupleDataDir+"/"+tupleDataName);
   //if(trigCat == 10) totFileName = "/vols/lhcb/th1011/KeeTuples/data/presel/B2Kee_Strip21_data_JpsiKPresel.root";

   TFile f( totFileName.c_str());


   TTree* t = (TTree*)f.Get(treeDataName.c_str());
   if(!t)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, no tree "<<treeDataName<<" found in "<<tupleDataDir<<"/"<<tupleDataName<<endl;
      return;
   }

   //prepare the cut to apply to the tree
   
   string trigMassCut( ("(B_plus_DTFM_M_zero>"+d2s(B_plus_M_min_data)+" && B_plus_DTFM_M_zero <"+d2s(B_plus_M_max_data)).c_str() ); 
   trigMassCut += " && passTrigCat"+i2s(trigCat)+" > 0.5)";
  // if(trigCat == 0) trigMassCut += " && L0ETOSOnly_d > 0.5)";
  // if(trigCat == 1) trigMassCut += " && L0HTOSOnly_d > 0.5)";
  // if(trigCat == 2) trigMassCut += " && (B_plus_L0MuonDecision_TIS || B_plus_L0HadronDecision_TIS || B_plus_L0ElectronDecision_TIS || B_plus_L0PhotonDecision_TIS) && !K_Kst_L0HadronDecision_TOS && L0ETOSOnly_d < 0.5)";
  // if(trigCat == 3) trigMassCut += " && (B_plus_L0MuonDecision_TIS || B_plus_L0HadronDecision_TIS || B_plus_L0ElectronDecision_TIS || B_plus_L0PhotonDecision_TIS))";
  // if(trigCat == 10) trigMassCut += " && (L0ETOSOnly_d == 1 && B_plus_Hlt1Global_TIS && B_plus_Hlt2Global_TIS))";

   TTreeFormula ttfTrigMassCut("ttfTrigMassCut", trigMassCut.c_str(), t);

   //clone the tree

   string nameNewFile(totFileName);
   nameNewFile.insert(nameNewFile.size()-5, label);
   TFile f2(nameNewFile.c_str(), "RECREATE");
   TTree* t2 = t->CloneTree(0);

   //fill the cloned tree with sWeights

   double sig_sw;
   double bkg_sw;
   double piee_sw;
   //double rare_sw;
   double prc_sw;

   t2->Branch("sig_sw", &sig_sw, "sig_sw/D");
   t2->Branch("bkg_sw", &bkg_sw, "bkg_sw/D");
   t2->Branch("piee_sw", &piee_sw, "piee_sw/D");
   //t2->Branch("rare_sw", &rare_sw, "rare_sw/D");
   t2->Branch("prc_sw", &prc_sw, "prc_sw/D");

   if( wDataSet->sumEntries() == t->GetEntries( trigMassCut.c_str()  ) )
   {
      cout<<"Putting the sweights in the tree... "<<endl;
   }
   else
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, mismatch in the number of entries"<<endl;
      return;
   }

   int nEntries(t->GetEntries());

   if(extracut == "") extracut = "1";
   TTreeFormula ttfExtraCut( "ttfExtraCut", extracut.c_str(), t);

   int j(0);

   for(int i(0); i<nEntries; ++i)
   {
      t->GetEntry(i);

      if(i % (nEntries/10) == 0) cout<<100*i/nEntries<<"\% "<<flush;

      if(ttfTrigMassCut.EvalInstance()) 
      {
         sig_sw = wDataSet->get(j)->getRealValue(("sig"+trigCatS+"_sw").c_str());
         bkg_sw = wDataSet->get(j)->getRealValue(("bkg"+trigCatS+"_sw").c_str());
         piee_sw = wDataSet->get(j)->getRealValue(("pieeVar"+trigCatS+"_sw").c_str());
         //rare_sw = wDataSet->get(j)->getRealValue(("rareVar"+trigCatS+"_sw").c_str());
         prc_sw = wDataSet->get(j)->getRealValue(("prc"+trigCatS+"_sw").c_str());
         ++j;

         if(ttfExtraCut.EvalInstance()) t2->Fill();
      }
   }

   cout<<endl;

   f2.cd();
   t2->Write();

   fw.Close();
   f.Close();
   f2.Close();
}

void FitAndSplotKeeDataForTraining::computeSWeight(int trigCat)
{
   string trigCatS("Trig"+i2s(trigCat));

   //get workspace, data, model

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspaceFit = (RooWorkspace*)fw.Get(("workspaceDataFit"+trigCatS).c_str());

   workspaceFit->Print();
   if(!workspaceFit)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::makeSWeightedTree, no workspace found in file "<<workspaceFileName<<endl;
      return;
   }

   RooDataSet* data = (RooDataSet*)workspaceFit->data(("data"+trigCatS).c_str());
   RooAbsPdf* model_total = workspaceFit->pdf(("modelTotForSWeight"+trigCatS).c_str());
   RooRealVar* sig = workspaceFit->var(("sig"+trigCatS).c_str());
   RooRealVar* bkg = workspaceFit->var(("bkg"+trigCatS).c_str());
   RooRealVar* prc = workspaceFit->var(("prc"+trigCatS).c_str());
   RooRealVar* pieeVar = workspaceFit->var(("pieeVar"+trigCatS).c_str());
   //RooRealVar* rareVar = workspaceFit->var(("rareVar"+trigCatS).c_str());

   model_total->graphVizTree( (plotdir+"/modeltree.dot").c_str() );

   if(!data || !model_total || !sig || !bkg || !prc || !pieeVar)// || !rareVar)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining::computeSWeight, error downloading stuff from workspace"<<endl;
      cout<<data<<" "<<model_total<<" "<<sig<<" "<<bkg<<" "<<prc<<" "<<pieeVar<<" "<<endl;
      return;
   }

   //compute the sWeights, get the sW data

   cout<<"Computing the sWeights..."<<endl;
   RooStats::SPlot wdata("wData", "wData", *data, model_total, RooArgList(*sig, *bkg, *prc, *pieeVar));

   cout<<"Recovering the sWeighted data set"<<endl;
   RooDataSet* wDataSet = wdata.GetSDataSet();

   wDataSet->SetName( ("wDataSet"+trigCatS).c_str() );
   wDataSet->SetTitle( ("wDataSet"+trigCatS).c_str() );

   //print some checks 

   cout<<"***********************************"<<endl<<"sWeight checks"<<endl<<"***********************************"<<endl;

   cout<<"sig from fit: "<<sig->getVal()<<", sig from sWeight: "<<wdata.GetYieldFromSWeight(("sig"+trigCatS).c_str())<<endl;
   cout<<"comb from fit: "<<bkg->getVal()<<", comb from sWeight: "<<wdata.GetYieldFromSWeight(("bkg"+trigCatS).c_str())<<endl;
   cout<<"prc from fit: "<<prc->getVal()<<", prc from sWeight: "<<wdata.GetYieldFromSWeight(("prc"+trigCatS).c_str())<<endl;
   cout<<"piee from fit: "<<pieeVar->getVal()<<", piee from sWeight: "<<wdata.GetYieldFromSWeight(("pieeVar"+trigCatS).c_str())<<endl;
   //cout<<"rare from fit: "<<rareVar->getVal()<<", rare from sWeight: "<<wdata.GetYieldFromSWeight(("rareVar"+trigCatS).c_str())<<endl;

   cout<<"***********************************"<<endl<<"check done"<<endl<<"***********************************"<<endl;

   //save in the workspace 

   RooWorkspace workspaceSPlot(("workspaceSPlot"+trigCatS).c_str(), ("workspaceSPlot"+trigCatS).c_str());
   workspaceSPlot.import(*wDataSet);

   fw.cd();
   workspaceSPlot.Write("", TObject::kOverwrite);

   fw.Close();
}


void FitAndSplotKeeDataForTraining::mergeAllTriggers()
{
   TFile fw(workspaceFileName.c_str(), "UPDATE");


   RooWorkspace* ws0 = (RooWorkspace*)fw.Get("workspaceDataFitTrig0");
   RooWorkspace* ws1 = (RooWorkspace*)fw.Get("workspaceDataFitTrig1");
   RooWorkspace* ws2 = (RooWorkspace*)fw.Get("workspaceDataFitTrig2");

   RooAbsPdf* modelBkg0 = ws0->pdf("modelBkgTotTrig0");
   RooAbsPdf* modelBkg1 = ws1->pdf("modelBkgTotTrig1");
   RooAbsPdf* modelBkg2 = ws2->pdf("modelBkgTotTrig2");

   RooRealVar* B_plus_DTFM_M_zero = ws0->var("B_plus_DTFM_M_zero");

   RooAbsPdf* modelSig0 = ws0->pdf("modelSigTotTrig0");
   RooAbsPdf* modelSig1 = ws1->pdf("modelSigTotTrig1");
   RooAbsPdf* modelSig2 = ws2->pdf("modelSigTotTrig2");

   RooRealVar* sig0 = (RooRealVar*)ws0->var("sigTrig0");
   RooRealVar* sig1 = (RooRealVar*)ws1->var("sigTrig1");
   RooRealVar* sig2 = (RooRealVar*)ws2->var("sigTrig2");

   RooRealVar* bkg0 = (RooRealVar*)ws0->var("bkgTrig0");
   RooRealVar* bkg1 = (RooRealVar*)ws1->var("bkgTrig1");
   RooRealVar* bkg2 = (RooRealVar*)ws2->var("bkgTrig2");

   RooDataSet* data = (RooDataSet*)ws0->data("dataTrig0");
   RooDataSet* data1 = (RooDataSet*)ws1->data("dataTrig1");
   RooDataSet* data2 = (RooDataSet*)ws2->data("dataTrig2");   

   RooRealVar* prc0 = (RooRealVar*)ws0->var(("prcTrig0"));
   RooRealVar* prc1 = (RooRealVar*)ws1->var(("prcTrig1"));
   RooRealVar* prc2 = (RooRealVar*)ws2->var(("prcTrig2"));

   RooRealVar* fracCharm0 = (RooRealVar*)ws0->var( ("fracCharmTrig0") );
   RooRealVar* fracCharm1 = (RooRealVar*)ws1->var( ("fracCharmTrig1") );
   RooRealVar* fracCharm2 = (RooRealVar*)ws2->var( ("fracCharmTrig2") );

   RooAbsPdf* modelCharm0 = ws0->pdf( ("modelCharmTrig0"));
   RooAbsPdf* modelCharm1 = ws1->pdf( ("modelCharmTrig1"));
   RooAbsPdf* modelCharm2 = ws2->pdf( ("modelCharmTrig2"));

   RooAbsPdf* modelRare0 = ws0->pdf( ("modelRareTrig0"));
   RooAbsPdf* modelRare1 = ws1->pdf( ("modelRareTrig1"));
   RooAbsPdf* modelRare2 = ws2->pdf( ("modelRareTrig2"));

   data->append(*data1); 
   data->append(*data2); 
   data->SetName("dataTot");


   RooFormulaVar fracSig0("fracSig0", "fracSig0", "@0/(@0+@1+@2)", RooArgList(*sig0, *sig1, *sig2));
   RooFormulaVar fracSig1("fracSig1", "fracSig1", "@1/(@0+@1+@2)", RooArgList(*sig0, *sig1, *sig2));
   RooFormulaVar fracSig2("fracSig2", "fracSig2", "@2/(@0+@1+@2)", RooArgList(*sig0, *sig1, *sig2));

   RooFormulaVar fracBkg0("fracBkg0", "fracBkg0", "@0/(@0+@1+@2)", RooArgList(*bkg0, *bkg1, *bkg2));
   RooFormulaVar fracBkg1("fracBkg1", "fracBkg1", "@1/(@0+@1+@2)", RooArgList(*bkg0, *bkg1, *bkg2));
   RooFormulaVar fracBkg2("fracBkg2", "fracBkg2", "@2/(@0+@1+@2)", RooArgList(*bkg0, *bkg1, *bkg2));


//   RooFormulaVar sig("sig", "sig", "@0+@1+@2", RooArgList(*sig0, *sig1, *sig2));
//   RooFormulaVar bkg("bkg", "bkg", "@0+@1+@2", RooArgList(*bkg0, *bkg1, *bkg2));
   RooRealVar sig("sig", "sig", sig0->getVal()+sig1->getVal()+sig2->getVal());
   RooRealVar bkg("bkg", "bkg", bkg0->getVal()+bkg1->getVal()+bkg2->getVal());

   RooRealVar nCharm0("nCharm0", "nCharm0", fracCharm0->getVal()*prc0->getVal() );
   RooRealVar nRare0("nRare0", "nRare0", prc0->getVal() - nCharm0.getVal() );
   RooRealVar nCharm1("nCharm1", "nCharm1", fracCharm1->getVal()*prc1->getVal() );
   RooRealVar nRare1("nRare1", "nRare1", prc1->getVal() - nCharm1.getVal() );
   RooRealVar nCharm2("nCharm2", "nCharm2", fracCharm2->getVal()*prc2->getVal() );
   RooRealVar nRare2("nRare2", "nRare2", prc2->getVal() - nCharm2.getVal() );


   RooFormulaVar nCharm("nCharm", "nCharm", "(@0+@1+@2)", RooArgList( nCharm0, nCharm1, nCharm2) );
   RooFormulaVar nRare("nRare", "nRare", "(@0+@1+@2)", RooArgList( nRare0, nRare1, nRare2) );

   RooFormulaVar fracCharmTrig0( "fracCharmTrig0", "fracCharmTrig0", "@0/(@0+@1+@2)", RooArgList( nCharm0, nCharm1, nCharm2) );
   RooFormulaVar fracCharmTrig1( "fracCharmTrig1", "fracCharmTrig1", "@1/(@0+@1+@2)", RooArgList( nCharm0, nCharm1, nCharm2) );


   RooFormulaVar fracRareTrig0( "fracRareTrig0", "fracRareTrig0", "@0/(@0+@1+@2)", RooArgList( nRare0, nRare1, nRare2) );
   RooFormulaVar fracRareTrig1( "fracRareTrig1", "fracRareTrig1", "@1/(@0+@1+@2)", RooArgList( nRare0, nRare1, nRare2) );

   RooAddPdf modelSigTot("modelSigTot", "modelSigTot", RooArgList(*modelSig0, *modelSig1, *modelSig2), RooArgList(fracSig0, fracSig1));
   RooAddPdf modelBkgTot("modelBkgTot", "modelBkgTot", RooArgList(*modelBkg0, *modelBkg1, *modelBkg2), RooArgList(fracBkg0, fracBkg1));

   RooAddPdf modelCharmTot("modelCharmTot", "modelCharmTot", RooArgList(*modelCharm0, *modelCharm1, *modelCharm2), RooArgList(fracCharmTrig0, fracCharmTrig1) );
   RooAddPdf modelRareTot("modelRareTot", "modelRareTot", RooArgList(*modelRare0, *modelRare1, *modelRare2), RooArgList(fracRareTrig0, fracRareTrig1) );



   RooAddPdf modelTot("modelTot", "modelTot", RooArgList(modelSigTot, modelBkgTot, modelCharmTot, modelRareTot), RooArgList(sig, bkg, nCharm, nRare));



   //plot


   RooPlot* frame = B_plus_DTFM_M_zero->frame();
   data->plotOn(frame);
   modelTot.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");

   ofstream out((plotdir+"plotKeeData.dat").c_str());
   saveFitInfo(out, frame, 6, &modelTot);
   out.close();

   modelTot.plotOn(frame, RooFit::LineColor(kGreen), RooFit::Components("modelBkgTot"));
   modelTot.plotOn(frame, RooFit::LineColor(kBlue), RooFit::Components("modelSigTot"));
   modelTot.plotOn(frame, RooFit::LineColor(kMagenta), RooFit::Components("modelCharmTot"));
   modelTot.plotOn(frame, RooFit::LineColor(kOrange), RooFit::Components("modelRareTot"));

   TCanvas canv("canv", "canv", 600, 600);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotKeeData.pdf").c_str());
   canvTot.Print((plotdir+"plotKeeData.root").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotKeeDataLogy.pdf").c_str());
   canvTot.Print((plotdir+"plotKeeDataLogy.root").c_str());

   //save the pdf in workspace 

   fw.cd();
   RooWorkspace ws3("workspaceDataFitMerged", "workspaceDataFitMerged");
   ws3.import(modelTot);
   ws3.import(*data);

   ws3.Write("", TObject::kOverwrite);
   

   fw.Close();
   fpull.Close();
}


void FitAndSplotKeeDataForTraining::fitData1PhotCat(int trigCat, int PhotCat, RooRealVar& yield)
{
   TFile fw(workspaceFileName.c_str(), "UPDATE");
   cout<<"fitting data for trig "<<trigCat<<" and phot cat "<<PhotCat<<endl;

   string cat("Trig"+i2s(trigCat)+"Phot"+i2s(PhotCat));
   string trigCatS("Trig"+i2s(trigCat));
   string trigCatSPrc("Trig"+i2s(trigCat)+"Phot-1");

   RooWorkspace* workspaceMCFit = (RooWorkspace*)fw.Get(("workspaceMCFit"+cat).c_str());
   cout<<"Hello0: "<<workspaceMCFit<<endl;

   string bremCut;
   if(PhotCat == 0) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) == 0)";
   if(PhotCat == 1) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) == 1)";
   if(PhotCat == 2) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) >= 2)";

   RooWorkspace* workspaceDataData = (RooWorkspace*)fw.Get( ("workspaceDataData"+trigCatS).c_str()  );
   if(!workspaceDataData)
   {
      cerr<<"ERROR: in function fitKeeData1PhotCat, no workspace found in file "<<workspaceFileName<<endl;
      return;
   }

   RooRealVar* B_plus_M = workspaceDataData->var("B_plus_M");
   RooDataSet* data = (RooDataSet*)workspaceDataData->data(("data"+trigCatS).c_str());

   // RooWorkspace* workspaceCharmPrc = (RooWorkspace*)fw.Get(("workspaceCharmPrcFit"+cat).c_str());
   // RooWorkspace* workspaceRarePrc = (RooWorkspace*)fw.Get(("workspaceRarePrcFit"+cat).c_str());
   RooWorkspace* workspaceLb = (RooWorkspace*)fw.Get(("workspaceLbFit"+trigCatS).c_str());
   RooWorkspace* workspaceCharmPrc = (RooWorkspace*)fw.Get(("workspaceCharmPrcFit"+trigCatSPrc).c_str());
   RooWorkspace* workspaceRarePrc = (RooWorkspace*)fw.Get(("workspaceRarePrcFit"+trigCatSPrc).c_str());
   RooWorkspace* workspacePiee = (RooWorkspace*)fw.Get(("workspacePieeFit"+cat).c_str());
   RooWorkspace* workspaceRare = (RooWorkspace*)fw.Get(("workspaceRareFit"+cat).c_str());

   if(!data || !B_plus_M)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining, error downloading stuff from workspace"<<endl;
      return;
   }

   data = (RooDataSet*)data->reduce( bremCut.c_str() );
   
   //get Kee MC fitvar

   RooRealVar* mean = (RooRealVar*)workspaceMCFit->var(("mean"+cat).c_str());
   cout<<"Hello3: "<<mean<<endl;
   mean->setConstant(true);
   RooRealVar* mean2 = (RooRealVar*)workspaceMCFit->var(("mean2"+cat).c_str());
   cout<<"Hello3: "<<mean2<<endl;
   mean2->setConstant(true);
   RooRealVar* mean3 = (RooRealVar*)workspaceMCFit->var(("mean3"+cat).c_str());
   cout<<"Hello3: "<<mean3<<endl;
   mean3->setConstant(true);
   RooRealVar* mean4 = (RooRealVar*)workspaceMCFit->var(("mean4"+cat).c_str());
   cout<<"Hello3: "<<mean4<<endl;
   mean4->setConstant(true);
   RooRealVar* sigma = (RooRealVar*)workspaceMCFit->var(("sigma"+cat).c_str());
   cout<<"Hello4: "<<sigma<<endl;
   sigma->setConstant(true);
   RooRealVar* sigma2 = (RooRealVar*)workspaceMCFit->var(("sigma2"+cat).c_str());
   cout<<"Hello4: "<<sigma2<<endl;
   sigma2->setConstant(true);
   RooRealVar* sigma3 = (RooRealVar*)workspaceMCFit->var(("sigma3"+cat).c_str());
   cout<<"Hello4: "<<sigma3<<endl;
   sigma3->setConstant(true);
   RooRealVar* sigma4 = (RooRealVar*)workspaceMCFit->var(("sigma4"+cat).c_str());
   cout<<"Hello4: "<<sigma4<<endl;
   sigma4->setConstant(true);
//   RooRealVar* fracSigma0 = (RooRealVar*)workspaceMCFit0->var(("fracSigma"+cat0).c_str());
//   cout<<"Hello5: "<<fracSigma0<<endl;
//   fracSigma0->setConstant(true);
   RooRealVar* al = (RooRealVar*)workspaceMCFit->var(("alpha_left"+cat).c_str());
   cout<<"Hello6: "<<al<<endl;
   al->setConstant(true);
   RooRealVar* ar = (RooRealVar*)workspaceMCFit->var(("alpha_right"+cat).c_str());
   cout<<"Hello7: "<<ar<<endl;
   ar->setConstant(true);
   RooRealVar* nl = (RooRealVar*)workspaceMCFit->var(("n_left"+cat).c_str());
   cout<<"Hello8: "<<nl<<endl;
   nl->setConstant(true);
   RooRealVar* nr = (RooRealVar*)workspaceMCFit->var(("n_right"+cat).c_str());
   cout<<"Hello9: "<<nr<<endl;
   nr->setConstant(true);
   RooRealVar* fracCB2 = (RooRealVar*)workspaceMCFit->var(("fracCB2"+cat).c_str());
   cout<<"Hello1: "<<fracCB2<<endl;
   fracCB2->setConstant(true);
   RooRealVar* fracGaus = (RooRealVar*)workspaceMCFit->var(("fracGaus"+cat).c_str());
   cout<<"Hello1: "<<fracGaus<<endl;
   fracGaus->setConstant(true);
   RooRealVar* fracGaus2 = (RooRealVar*)workspaceMCFit->var(("fracGaus2"+cat).c_str());
   cout<<"Hello1: "<<fracGaus2<<endl;
   fracGaus2->setConstant(true);
   RooRealVar* fracCB3 = (RooRealVar*)workspaceMCFit->var(("fracCB3"+cat).c_str());
   cout<<"Hello1: "<<fracCB3<<endl;
   fracCB3->setConstant(true);
   RooRealVar* nEvents = (RooRealVar*)workspaceMCFit->var(("nEvents"+cat).c_str());
   cout<<"Hello11: "<<nEvents<<endl;
   nEvents->setConstant(true);
   // RooKeysPdf* charmShape = (RooKeysPdf*)workspaceCharmPrc->pdf(("templateCharmPrc"+cat).c_str());
   // RooRealVar* nMCCharmPrc = (RooRealVar*)workspaceCharmPrc->var( ("nEventsInDataRangeCharmPrc"+cat).c_str() );
   // RooKeysPdf* rareShape = (RooKeysPdf*)workspaceRarePrc->pdf(("templateRarePrc"+cat).c_str());
   // RooRealVar* nMCRarePrc = (RooRealVar*)workspaceRarePrc->var( ("nEventsInDataRangeRarePrc"+cat).c_str() );
   // cout<<"Hello12: "<<charmShape<<" "<<nMCCharmPrc<<" "<<rareShape<<" "<<nMCRarePrc<<endl;

   RooKeysPdf* LbShape = (RooKeysPdf*)workspaceLb->pdf(("templateLb"+trigCatS).c_str());
   RooKeysPdf* charmShape = (RooKeysPdf*)workspaceCharmPrc->pdf(("templateCharmPrc"+trigCatSPrc).c_str());
   RooRealVar* nMCLb = (RooRealVar*)workspaceLb->var( ("nEventsInDataRangeLb"+trigCatS).c_str() );
   RooRealVar* nMCCharmPrc = (RooRealVar*)workspaceCharmPrc->var( ("nEventsInDataRangeCharmPrc"+trigCatSPrc).c_str() );
   RooKeysPdf* rareShape = (RooKeysPdf*)workspaceRarePrc->pdf(("templateRarePrc"+trigCatSPrc).c_str());
   RooRealVar* nMCRarePrc = (RooRealVar*)workspaceRarePrc->var( ("nEventsInDataRangeRarePrc"+trigCatSPrc).c_str() );
   cout<<"Hello12: "<<charmShape<<" "<<nMCCharmPrc<<" "<<LbShape<<" "<<nMCLb<<" "<<rareShape<<" "<<nMCRarePrc<<endl;

   //get Piee MC fit variables
   

   RooRealVar* meanPiee = (RooRealVar*)workspacePiee->var((("meanPiee")+cat).c_str());
   cout<<"Helloa"<<endl;
   meanPiee->setConstant(true);
   cout<<"Helloa"<<endl;
   RooRealVar* sigmaPiee = (RooRealVar*)workspacePiee->var((("sigmaPiee")+cat).c_str());
   sigmaPiee->setConstant(true);
   cout<<"Hellob"<<endl;
   RooRealVar* alPiee = (RooRealVar*)workspacePiee->var((("alpha_leftPiee")+cat).c_str());
   alPiee->setConstant(true);
   cout<<"Hellod"<<endl;
   RooRealVar* arPiee = (RooRealVar*)workspacePiee->var((("alpha_rightPiee")+cat).c_str());
   arPiee->setConstant(true);
   cout<<"Helloe"<<endl;
   RooRealVar* nlPiee = (RooRealVar*)workspacePiee->var((("n_leftPiee")+cat).c_str());
   nlPiee->setConstant(true);
   cout<<"Hellof"<<endl;
   RooRealVar* nrPiee = (RooRealVar*)workspacePiee->var((("n_rightPiee")+cat).c_str());
   nrPiee->setConstant(true);
   cout<<"Hellog"<<endl;
   RooRealVar* fracCB2Piee = (RooRealVar*)workspacePiee->var((("fracCB2Piee")+cat).c_str());
   fracCB2Piee->setConstant(true);
   cout<<"Helloh"<<endl;

   //get the rare model

   // RooChebychev*  rarePdf = (RooChebychev*)workspaceRare->pdf(("poly"+cat).c_str());

   //Construct the signal model
   

   RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+cat).c_str(), ("sigmaScaleFactor"+cat).c_str(), 1.12, 0.9, 1.5);
   RooRealVar meanShift(("meanShift"+cat).c_str(), ("meanShift"+cat).c_str(), 0, -100, 100);

   RooFormulaVar shiftedMean(("shiftedMean"+cat).c_str(), ("shiftedMean"+cat).c_str(), "@0+@1", RooArgList(*mean, meanShift));
   RooFormulaVar shiftedMean2(("shiftedMean2"+cat).c_str(), ("shiftedMean2"+cat).c_str(), "@0+@1", RooArgList(*mean2, meanShift));
   RooFormulaVar shiftedMean3(("shiftedMean3"+cat).c_str(), ("shiftedMean3"+cat).c_str(), "@0+@1", RooArgList(*mean3, meanShift));
   RooFormulaVar shiftedMean4(("shiftedMean4"+cat).c_str(), ("shiftedMean4"+cat).c_str(), "@0+@1", RooArgList(*mean4, meanShift));
   RooFormulaVar scaledSigma(("scaledSigma"+cat).c_str(), ("scaledSigma"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma));
   RooFormulaVar scaledSigma2(("scaledSigma2"+cat).c_str(), ("scaledSigma2"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma2));
   RooFormulaVar scaledSigma4(("scaledSigma4"+cat).c_str(), ("scaledSigma4"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma4));
   RooFormulaVar scaledSigma3(("scaledSigma3"+cat).c_str(), ("scaledSigma3"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma3));
   
   RooRealVar arScale("arScale", "arScale", 1,0,3);
   arScale.setConstant(true);
   RooFormulaVar scaledAr(("scaledAr"+cat).c_str(), ("scaledAr"+cat).c_str(), "@0*@1", RooArgList(*ar, arScale));

   RooCBShape cb1(("cb1"+cat).c_str(), ("cb1"+cat).c_str(), *B_plus_M, shiftedMean, scaledSigma, *al, *nl);
   // RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, shiftedMean, scaledSigma2, scaledAr, *nr);
   RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, shiftedMean2, scaledSigma2, scaledAr, *nr);
   
   RooGaussian gaus(("gaus"+cat).c_str(), ("gaus"+cat).c_str(), *B_plus_M, shiftedMean3, scaledSigma3 );
   RooGaussian gaus2(("gaus2"+cat).c_str(), ("gaus2"+cat).c_str(), *B_plus_M, shiftedMean4, scaledSigma4 );
  
   RooAbsPdf *modelSigTot;

   if (trigCat==0)
   {
     
     if (PhotCat==0)
     {
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), 
                                 RooArgList(cb2,cb1),RooArgList(*fracCB2),true); 
     
       // modelSigTot = new RooCBShape(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), *B_plus_M, shiftedMean, scaledSigma,*al, *nl);
     }
     if (PhotCat==1)
       // if (PhotCat==1 || PhotCat==2)
     {
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb2, gaus, gaus2, cb1) ,
                                   RooArgList(*fracCB2, *fracGaus, *fracGaus2), true);
     }
     if (PhotCat==2)
     {
       // modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb1, gaus_mean2) ,RooArgList(*fracCB2));
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb2, gaus, gaus2, cb1),
                                   RooArgList(*fracCB2, *fracGaus, *fracGaus2),true);
     }
   }

   if (trigCat==1)
   {  
     if (PhotCat==0)
     {
       // modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), 
                                 // RooArgList(cb1,gaus),RooArgList(*fracCB2),true);
       modelSigTot = new RooCBShape(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), 
                                    *B_plus_M, shiftedMean, scaledSigma,*al, *nl);
     }
     if (PhotCat==1)
       // if (PhotCat==1 || PhotCat==2)
     {
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb2, cb1) ,
                                   RooArgList(*fracCB2), true);
     }
     if (PhotCat==2)
     {
       // modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb1, gaus_mean2) ,RooArgList(*fracCB2));
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb2, gaus, cb1),
                                   RooArgList(*fracCB2, *fracGaus2),true);
     }
   }

   if (trigCat==2)
   {  
     if (PhotCat==0)
     {
       // modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), 
                                 // RooArgList(cb1,gaus),RooArgList(*fracCB2),true); 
     
       modelSigTot = new RooCBShape(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), 
                                    *B_plus_M, shiftedMean, scaledSigma,*al, *nl);
     
     }
     if (PhotCat==1)
       // if (PhotCat==1 || PhotCat==2)
     {
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb1, gaus) ,
                                   RooArgList(*fracCB2), true);
     }
     if (PhotCat==2)
     {
       // modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb1, gaus_mean2) ,RooArgList(*fracCB2));
       modelSigTot = new RooAddPdf(("modelSigTot"+cat).c_str(), ("modelSigTot"+cat).c_str(), RooArgList(cb2, gaus, cb1),
                                   RooArgList(*fracCB2, *fracGaus),true);
     }
   }
   
   //Construct the Piee model

   RooFormulaVar shiftedMeanPiee(("shiftedMeanPiee"+cat).c_str(), ("shiftedMeanPiee"+cat).c_str(), "@0+@1", 
                                 RooArgList(*meanPiee,meanShift));
   RooFormulaVar scaledSigmaPiee(("scaledSigmaPiee"+cat).c_str(), ("scaledSigmaPiee"+cat).c_str(), "@0*@1", 
                                 RooArgList(sigmaScaleFactor, *sigmaPiee));
   
   RooFormulaVar scaledArPiee(("scaledArPiee"+cat).c_str(), ("scaledArPiee"+cat).c_str(), "@0*@1", RooArgList(*arPiee, arScale));


   RooCBShape cb1Piee(("cb1Piee"+cat).c_str(), ("cb1Piee"+cat).c_str(), *B_plus_M, shiftedMeanPiee, scaledSigmaPiee, *alPiee, *nlPiee);
   RooCBShape cb2Piee(("cb2Piee"+cat).c_str(), ("cb2Piee"+cat).c_str(), *B_plus_M, shiftedMeanPiee, scaledSigmaPiee, 
                      scaledArPiee, *nrPiee);
   
   RooAddPdf modelPieeTot(("modelPieeTot"+cat).c_str(), ("modelPieeTot"+cat).c_str(), RooArgList(cb1Piee, cb2Piee ) ,
                          RooArgList(*fracCB2Piee));

   //construct the Prc model

   RooRealVar fracCharmPrc(("fracCharmPrc"+cat).c_str(), ("fracCharmPrc"+cat).c_str(), 
                           nMCCharmPrc->getVal()/(nMCCharmPrc->getVal()+nMCRarePrc->getVal()),   0,1);
   // fracCharmPrc.setVal(nMCCharmPrc->getVal()/(nMCCharmPrc->getVal()+nMCRarePrc->getVal()+0.25*nMCLb->getVal() ));
   // fracCharmPrc.setVal(0.);
   fracCharmPrc.setConstant(true);

   RooRealVar fracLb(("fracLb"+cat).c_str(), ("fracLb"+cat).c_str(), 
                     0.25*nMCLb->getVal()/(0.25*nMCLb->getVal()+nMCCharmPrc->getVal()+nMCRarePrc->getVal() ),   0,1);

   // RooAddPdf modelPrcTot(("modelPrcTot"+cat).c_str(), ("modelPrcTot"+cat).c_str(), RooArgList(*charmShape, *LbShape, *rareShape), RooArgList(fracCharmPrc, fracLb) );
   RooAddPdf modelPrcTot(("modelPrcTot"+cat).c_str(), ("modelPrcTot"+cat).c_str(), RooArgList(*charmShape, *rareShape), 
                         RooArgList(fracCharmPrc) );
   
   
   cout<<endl<<endl<<"MEUH MEUH MEUH "<<fracCharmPrc.getVal()<<" "<<data->sumEntries()<<endl<<endl;

   //construct the bkg model
   
   RooRealVar lambda(("lambda"+cat).c_str(),("lambda"+cat).c_str(),-5e-04,-0.01, 0.01);
   RooExponential expo(("modelBkgTot"+cat).c_str(), ("modelBkgTot"+cat).c_str(), *B_plus_M,  lambda);
   // // RooRealVar lambda(("lambda"+cat).c_str(),("lambda"+cat).c_str(),-1.,1.);
   // // RooChebychev comb(("modelComb"+cat).c_str(), ("modelComb"+cat).c_str(), *B_plus_M, RooArgList(lambda));
   // RooExponential comb(("modelComb"+cat).c_str(), ("modelComb"+cat).c_str(), *B_plus_M,  lambda);
   // // RooRealVar cheby1(("cheby1"+cat).c_str(),("cheby1"+cat).c_str(),-1.,1);
   // // RooRealVar cheby2(("cheby2"+cat).c_str(),("cheby2"+cat).c_str(),-1.,1.);
   // // RooChebychev Lb(("modelLb"+cat).c_str(), ("modelLb"+cat).c_str(), *B_plus_M, RooArgList(cheby1,cheby2));
   // RooRealVar cheby1(("cheby1"+trigCatS).c_str(),("cheby1"+trigCatS).c_str(),5200,5500);
   // RooRealVar cheby2(("cheby2"+trigCatS).c_str(),("cheby2"+trigCatS).c_str(),50,2000);
   // RooGaussian Lb(("modelLb"+trigCatS).c_str(), ("modelLb"+trigCatS).c_str(), *B_plus_M, cheby1,cheby2);
   // RooRealVar fracLbcomb(("fracLbcomb"+cat).c_str(),("fracLbcomb"+cat).c_str(),0,1);

   // RooAddPdf expo(("modelBkgTot"+cat).c_str(), ("modelBkgTot"+cat).c_str(), Lb, comb, fracLbcomb);


   //prepare the rare constraint

   RooRealVar fracRareMu("fracRareMu", "fracRareMu", 0, 1);

   fracRareMu.SetName(("fracRareMu"+cat).c_str()); fracRareMu.SetTitle(("fracRareMu"+cat).c_str());
   fracRareMu.setConstant(true);
   RooRealVar fracRareSigma( ("fracRareSigma"+cat).c_str(), ("fracRareSigma"+cat).c_str(),  fracRareMu.getError());
   fracRareSigma.setConstant(true);

   //construct the model total

   RooRealVar sig(("sig"+cat).c_str(), ("sig"+cat).c_str(),1.7e4, 0, 200000);
   RooRealVar bkg(("bkg"+cat).c_str(), ("bkg"+cat).c_str(),1e5, 0, 200000);
   RooRealVar prc(("prc"+cat).c_str(), ("prc"+cat).c_str(),1e3, 0, 200000);

   RooRealVar fracPiee(("fracpiee"+cat).c_str(), ("fracpiee"+cat).c_str(), 0, 1);
   RooFormulaVar piee(("piee"+cat).c_str(), ("piee"+cat).c_str(), "@0*@1", RooArgList(fracPiee,sig) ); 

   RooRealVar fracRare("fracRare", "fracRare", fracRareMu.getVal(), 0,  fracRareMu.getVal()+10*fracRareSigma.getVal());
   RooFormulaVar rare(("rare"+cat).c_str(), ("rare"+cat).c_str(), "@0*@1", RooArgList(fracRare,sig) ); 

   // RooAddPdf modelTot(("modelTot"+cat).c_str(),("modelTot"+cat).c_str(), RooArgList( *modelSigTot, expo, modelPrcTot, modelPieeTot, *rarePdf ) ,RooArgList(sig, bkg, prc, piee, rare ));

   RooArgList components(*modelSigTot, expo, modelPieeTot);
   RooArgList yields(sig, bkg, piee);

   // if (nMCCharmPrc->getVal()>0) 
   // {
   //   components.add(modelPrcTot);
   //   yields.add(prc);
   // }
   // else
   // {
   //   components.add(*rareShape);
   //   yields.add(prc);  
   // }

   RooAddPdf modelTot(("modelTot"+cat).c_str(),("modelTot"+cat).c_str(), components , yields);
   
   // RooAddPdf modelTot(("modelTot"+cat).c_str(),("modelTot"+cat).c_str(), RooArgList( *modelSigTot, expo, modelPrcTot, modelPieeTot ) ,
                      // RooArgList(sig, bkg, prc, piee ));
   // RooAddPdf modelTot(("modelTot"+cat).c_str(),("modelTot"+cat).c_str(), RooArgList( *modelSigTot, expo, modelPrcTot ) ,RooArgList(sig, bkg, prc ));
   //RooAddPdf modelTot(("modelTot"+cat).c_str(),("modelTot"+cat).c_str(), RooArgList( modelSigTot, expo ) ,RooArgList(sig, bkg));

   RooGaussian fracRareConstraint( ("fracRareConstraint"+cat).c_str(), ("fracRareConstraint"+cat).c_str(), 
                                   fracRare, fracRareMu, fracRareSigma);

   //prepare piee frac constraint

   RooRealVar fracPieeMu("fracPieeMu", "fracPieeMu", 0, 1);
   getFracPiee(trigCat, fracPieeMu);
   fracPieeMu.SetName(("fracPieeMu"+cat).c_str()); fracPieeMu.SetTitle(("fracPieeMu"+cat).c_str());
   fracPieeMu.setConstant(true);
   RooRealVar fracPieeSigma( ("fracPieeSigma"+cat).c_str(), ("fracPieeSigma"+cat).c_str(),  fracPieeMu.getError());
   fracPieeSigma.setConstant(true);

   RooGaussian fracPieeConstraint( ("fracPieeConstraint"+cat).c_str(), ("fracPieeConstraint"+cat).c_str(), fracPiee, fracPieeMu, 
                                   fracPieeSigma);

   fracPiee.setVal(fracPieeMu.getVal());
   fracPiee.setRange(fracPieeMu.getVal() - 10*fracPieeSigma.getVal(), fracPieeMu.getVal() + 10*fracPieeSigma.getVal());

   //fit

   // RooAbsReal* nll = modelTot.createNLL(*data, RooFit::NumCPU(8), ExternalConstraints(RooArgSet(fracPieeConstraint, fracRareConstraint) ));
   RooAbsReal* nll = modelTot.createNLL(*data, RooFit::NumCPU(8), ExternalConstraints(RooArgSet(fracPieeConstraint) ));
   // RooAbsReal* nll = modelTot.createNLL(*data, RooFit::NumCPU(8));

   RooMinuit m(*nll);
   int migradStatus, hesseStatus;


   migradStatus = m.migrad();
   hesseStatus = m.hesse();
   // m.minos();

   //plot

   B_plus_M->setBins(50);
   B_plus_M->SetTitle("M(Kee)");

   RooPlot* frame = B_plus_M->frame();
   data->plotOn(frame);
   modelTot.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");


   ofstream out((plotdir+"plotKeeData"+cat+".dat").c_str());
   saveFitInfo(out, frame, 6, &modelTot);

   modelTot.plotOn(frame, RooFit::LineColor(kGreen), RooFit::Components( ("modelBkgTot"+cat).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kBlue), RooFit::Components( ("modelSigTot"+cat).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kCyan), RooFit::Components( ("templateLb"+trigCatS).c_str()) );
   modelTot.plotOn(frame, RooFit::LineColor(kMagenta), RooFit::Components( ("templateCharmPrc"+trigCatSPrc).c_str()) );
   modelTot.plotOn(frame, RooFit::LineColor(kOrange), RooFit::Components( ("templateRarePrc"+trigCatSPrc).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kGray+2), RooFit::Components( ("modelPieeTot"+cat).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kViolet-6), RooFit::Components( ("poly"+cat).c_str() ));

   TCanvas canv("canv", "canv", 600, 600);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotKeeData"+cat+".pdf").c_str());
   canvTot.Print((plotdir+"plotKeeData"+cat+".root").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotKeeDataLogy"+cat+".pdf").c_str());
   canvTot.Print((plotdir+"plotKeeDataLogy"+cat+".root").c_str());

   //save infos

   yield.setVal(sig.getVal()); 
   yield.setError(sig.getError());

   RooFitResult* fitRes = m.save();
   double edm (fitRes->edm());
   out<<"Migrad status: "<<migradStatus<<endl;
   out<<"Hesse status: "<<hesseStatus<<endl;
   out<<"EDM: "<<edm<<endl;
   out.close();

   fw.cd();
   RooWorkspace ws2(("workspaceDataFit"+cat).c_str(), ("workspaceDataFit"+cat).c_str());
   ws2.import(modelTot);
   ws2.import(*data);

   ws2.Write("", TObject::kOverwrite);

   fw.Close();
   fpull.Close();



}


void FitAndSplotKeeDataForTraining::fitData(int trigCat, RooRealVar& yield, bool fast, bool saveYield, double shiftFix, double scaleFix)
{
   cout<<"fitting data"<<endl;


   string trigCatSData("Trig"+i2s(trigCat));

   //if(trigCat == 10) trigCat = 0;

   string cat0("Trig"+i2s(trigCat)+"Phot"+i2s(0));
   string cat1("Trig"+i2s(trigCat)+"Phot"+i2s(1));
   string cat2("Trig"+i2s(trigCat)+"Phot"+i2s(2));
   string trigCatS("Trig"+i2s(trigCat));
   string trigCatPrc("Trig"+i2s(trigCat)+"Phot-1");

   //Get the dataset

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspaceDataData = (RooWorkspace*)fw.Get( ("workspaceDataData"+trigCatSData).c_str()  );
   if(!workspaceDataData)
   {
      cerr<<"ERROR: in function fitKeeData, no workspace "<<("workspaceDataData"+trigCatSData).c_str()<< "found in file "<<workspaceFileName<<endl;
      return;
   }

   RooWorkspace* workspacePiee = (RooWorkspace*)fw.Get(("workspacePieeFit"+trigCatPrc).c_str());
   RooWorkspace* workspaceLb = (RooWorkspace*)fw.Get(("workspaceLbFit"+trigCatS).c_str());
   RooRealVar* B_plus_M = workspaceDataData->var("B_plus_M");
   RooDataSet* data = (RooDataSet*)workspaceDataData->data(("data"+trigCatSData).c_str());

   if(!data || !B_plus_M)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining, error downloading stuff from workspace"<<endl;
      return;
   }

   //Get the variables from the MC and PRC

   RooWorkspace* workspaceMCFit0 = (RooWorkspace*)fw.Get(("workspaceMCFit"+cat0).c_str());
   cout<<"Hello0: "<<workspaceMCFit0<<endl;
   RooWorkspace* workspaceMCFit1 = (RooWorkspace*)fw.Get(("workspaceMCFit"+cat1).c_str());
   cout<<"Hello1: "<<workspaceMCFit1<<endl;
   RooWorkspace* workspaceMCFit2 = (RooWorkspace*)fw.Get(("workspaceMCFit"+cat2).c_str());
   cout<<"Hello2: "<<workspaceMCFit2<<endl;

   
   RooWorkspace* workspaceCharmPrc = (RooWorkspace*)fw.Get(("workspaceCharmPrcFit"+trigCatPrc).c_str());
   RooWorkspace* workspaceRarePrc = (RooWorkspace*)fw.Get(("workspaceRarePrcFit"+trigCatPrc).c_str());
   //RooWorkspace* workspaceRare = (RooWorkspace*)fw.Get(("workspaceRareFit"+trigCatPrc).c_str());

   // Parameters trig0
   RooRealVar* mean0 = (RooRealVar*)workspaceMCFit0->var(("mean"+cat0).c_str());
   cout<<"Hello3: "<<mean0<<endl;
   mean0->setConstant(true);
   RooRealVar* mean20 = (RooRealVar*)workspaceMCFit0->var(("mean2"+cat0).c_str());
   cout<<"Hello3: "<<mean20<<endl;
   mean20->setConstant(true);
   RooRealVar* mean30 = (RooRealVar*)workspaceMCFit0->var(("mean3"+cat0).c_str());
   cout<<"Hello3: "<<mean30<<endl;
   mean30->setConstant(true);
   RooRealVar* mean40 = (RooRealVar*)workspaceMCFit0->var(("mean4"+cat0).c_str());
   cout<<"Hello3: "<<mean40<<endl;
   mean40->setConstant(true);

   RooRealVar* sigma0 = (RooRealVar*)workspaceMCFit0->var(("sigma"+cat0).c_str());
   cout<<"Hello4: "<<sigma0<<endl;
   sigma0->setConstant(true);
   RooRealVar* sigma20 = (RooRealVar*)workspaceMCFit0->var(("sigma2"+cat0).c_str());
   cout<<"Hello4: "<<sigma20<<endl;
   sigma20->setConstant(true);
   RooRealVar* sigma30 = (RooRealVar*)workspaceMCFit0->var(("sigma3"+cat0).c_str());
   cout<<"Hello4: "<<sigma30<<endl;
   sigma30->setConstant(true);
   RooRealVar* sigma40 = (RooRealVar*)workspaceMCFit0->var(("sigma4"+cat0).c_str());
   cout<<"Hello4: "<<sigma40<<endl;
   sigma40->setConstant(true);

   RooRealVar* al0 = (RooRealVar*)workspaceMCFit0->var(("alpha_left"+cat0).c_str());
   cout<<"Hello6: "<<al0<<endl;
   al0->setConstant(true);
   RooRealVar* ar0 = (RooRealVar*)workspaceMCFit0->var(("alpha_right"+cat0).c_str());
   cout<<"Hello7: "<<ar0<<endl;
   ar0->setConstant(true);
   RooRealVar* nl0 = (RooRealVar*)workspaceMCFit0->var(("n_left"+cat0).c_str());
   cout<<"Hello8: "<<nl0<<endl;
   nl0->setConstant(true);
   RooRealVar* nr0 = (RooRealVar*)workspaceMCFit0->var(("n_right"+cat0).c_str());
   cout<<"Hello9: "<<nr0<<endl;
   nr0->setConstant(true);

   RooRealVar* fracCB20 = (RooRealVar*)workspaceMCFit0->var(("fracCB2"+cat0).c_str());
   cout<<"Hello10: "<<fracCB20<<endl;
   fracCB20->setConstant(true);
   RooRealVar* fracGaus0 = (RooRealVar*)workspaceMCFit0->var(("fracGaus"+cat0).c_str());
   cout<<"Hello10: "<<fracGaus0<<endl;
   fracGaus0->setConstant(true);
   RooRealVar* fracGaus20 = (RooRealVar*)workspaceMCFit0->var(("fracGaus2"+cat0).c_str());
   cout<<"Hello10: "<<fracGaus20<<endl;
   fracGaus20->setConstant(true);
   RooRealVar* nEvents0 = (RooRealVar*)workspaceMCFit0->var(("nEvents"+cat0).c_str());
   cout<<"Hello11: "<<nEvents0<<endl;
   nEvents0->setConstant(true);
   RooRealVar* nEventsInDataRange0 = (RooRealVar*)workspaceMCFit0->var(("nEventsInDataRange"+cat0).c_str());
   cout<<"Hello12: "<<nEvents0<<endl;
   nEventsInDataRange0->setConstant(true);


   // Parameters trig1
   RooRealVar* mean1 = (RooRealVar*)workspaceMCFit1->var(("mean"+cat1).c_str());
   cout<<"Hello3: "<<mean1<<endl;
   mean1->setConstant(true);
   RooRealVar* mean21 = (RooRealVar*)workspaceMCFit1->var(("mean2"+cat1).c_str());
   cout<<"Hello3: "<<mean21<<endl;
   mean21->setConstant(true);
   RooRealVar* mean31 = (RooRealVar*)workspaceMCFit1->var(("mean3"+cat1).c_str());
   cout<<"Hello3: "<<mean31<<endl;
   mean31->setConstant(true);
   RooRealVar* mean41 = (RooRealVar*)workspaceMCFit1->var(("mean4"+cat1).c_str());
   cout<<"Hello3: "<<mean41<<endl;
   mean41->setConstant(true);

   RooRealVar* sigma1 = (RooRealVar*)workspaceMCFit1->var(("sigma"+cat1).c_str());
   cout<<"Hello4: "<<sigma1<<endl;
   sigma1->setConstant(true);
   RooRealVar* sigma21 = (RooRealVar*)workspaceMCFit1->var(("sigma2"+cat1).c_str());
   cout<<"Hello4: "<<sigma21<<endl;
   sigma21->setConstant(true);
   RooRealVar* sigma31 = (RooRealVar*)workspaceMCFit1->var(("sigma3"+cat1).c_str());
   cout<<"Hello4: "<<sigma31<<endl;
   sigma31->setConstant(true);
   RooRealVar* sigma41 = (RooRealVar*)workspaceMCFit1->var(("sigma4"+cat1).c_str());
   cout<<"Hello4: "<<sigma41<<endl;
   sigma41->setConstant(true);

   RooRealVar* al1 = (RooRealVar*)workspaceMCFit1->var(("alpha_left"+cat1).c_str());
   cout<<"Hello6: "<<al1<<endl;
   al1->setConstant(true);
   RooRealVar* ar1 = (RooRealVar*)workspaceMCFit1->var(("alpha_right"+cat1).c_str());
   cout<<"Hello7: "<<ar1<<endl;
   ar1->setConstant(true);
   RooRealVar* nl1 = (RooRealVar*)workspaceMCFit1->var(("n_left"+cat1).c_str());
   cout<<"Hello8: "<<nl1<<endl;
   nl1->setConstant(true);
   RooRealVar* nr1 = (RooRealVar*)workspaceMCFit1->var(("n_right"+cat1).c_str());
   cout<<"Hello9: "<<nr1<<endl;
   nr1->setConstant(true);

   RooRealVar* fracCB21 = (RooRealVar*)workspaceMCFit1->var(("fracCB2"+cat1).c_str());
   cout<<"Hello11: "<<fracCB21<<endl;
   fracCB21->setConstant(true);
   RooRealVar* fracGaus1 = (RooRealVar*)workspaceMCFit1->var(("fracGaus"+cat1).c_str());
   cout<<"Hello11: "<<fracGaus1<<endl;
   fracGaus1->setConstant(true);
   RooRealVar* fracGaus21 = (RooRealVar*)workspaceMCFit1->var(("fracGaus2"+cat1).c_str());
   cout<<"Hello11: "<<fracGaus21<<endl;
   fracGaus21->setConstant(true);
   RooRealVar* nEvents1 = (RooRealVar*)workspaceMCFit1->var(("nEvents"+cat1).c_str());
   cout<<"Hello11: "<<nEvents1<<endl;
   nEvents1->setConstant(true);
   RooRealVar* nEventsInDataRange1 = (RooRealVar*)workspaceMCFit1->var(("nEventsInDataRange"+cat1).c_str());
   cout<<"Hello12: "<<nEvents1<<endl;
   nEventsInDataRange1->setConstant(true);


   // Parameters trig2
   RooRealVar* mean2 = (RooRealVar*)workspaceMCFit2->var(("mean"+cat2).c_str());
   cout<<"Hello3: "<<mean2<<endl;
   mean2->setConstant(true);
   RooRealVar* mean22 = (RooRealVar*)workspaceMCFit2->var(("mean2"+cat2).c_str());
   cout<<"Hello3: "<<mean22<<endl;
   mean22->setConstant(true);
   RooRealVar* mean32 = (RooRealVar*)workspaceMCFit2->var(("mean3"+cat2).c_str());
   cout<<"Hello3: "<<mean32<<endl;
   mean32->setConstant(true);
   RooRealVar* mean42 = (RooRealVar*)workspaceMCFit2->var(("mean4"+cat2).c_str());
   cout<<"Hello3: "<<mean42<<endl;
   mean42->setConstant(true);

   RooRealVar* sigma2 = (RooRealVar*)workspaceMCFit2->var(("sigma"+cat2).c_str());
   cout<<"Hello4: "<<sigma2<<endl;
   sigma2->setConstant(true);
   RooRealVar* sigma22 = (RooRealVar*)workspaceMCFit2->var(("sigma2"+cat2).c_str());
   cout<<"Hello4: "<<sigma22<<endl;
   sigma22->setConstant(true);
   RooRealVar* sigma32 = (RooRealVar*)workspaceMCFit2->var(("sigma3"+cat2).c_str());
   cout<<"Hello4: "<<sigma32<<endl;
   sigma32->setConstant(true);
   RooRealVar* sigma42 = (RooRealVar*)workspaceMCFit2->var(("sigma4"+cat2).c_str());
   cout<<"Hello4: "<<sigma42<<endl;
   sigma42->setConstant(true);

   RooRealVar* al2 = (RooRealVar*)workspaceMCFit2->var(("alpha_left"+cat2).c_str());
   cout<<"Hello6: "<<al2<<endl;
   al2->setConstant(true);
   RooRealVar* ar2 = (RooRealVar*)workspaceMCFit2->var(("alpha_right"+cat2).c_str());
   cout<<"Hello7: "<<ar2<<endl;
   ar2->setConstant(true);
   RooRealVar* nl2 = (RooRealVar*)workspaceMCFit2->var(("n_left"+cat2).c_str());
   cout<<"Hello8: "<<nl2<<endl;
   nl2->setConstant(true);
   RooRealVar* nr2 = (RooRealVar*)workspaceMCFit2->var(("n_right"+cat2).c_str());
   cout<<"Hello9: "<<nr2<<endl;
   nr2->setConstant(true);

   RooRealVar* fracCB22 = (RooRealVar*)workspaceMCFit2->var(("fracCB2"+cat2).c_str());
   cout<<"Hello12: "<<fracCB22<<endl;
   fracCB22->setConstant(true);
   RooRealVar* fracGaus2 = (RooRealVar*)workspaceMCFit2->var(("fracGaus"+cat2).c_str());
   cout<<"Hello12: "<<fracGaus2<<endl;
   fracGaus2->setConstant(true);
   RooRealVar* fracGaus22 = (RooRealVar*)workspaceMCFit2->var(("fracGaus2"+cat2).c_str());
   cout<<"Hello12: "<<fracGaus22<<endl;
   fracGaus22->setConstant(true);
   RooRealVar* nEvents2 = (RooRealVar*)workspaceMCFit2->var(("nEvents"+cat2).c_str());
   cout<<"Hello11: "<<nEvents2<<endl;
   nEvents2->setConstant(true);
   RooRealVar* nEventsInDataRange2 = (RooRealVar*)workspaceMCFit2->var(("nEventsInDataRange"+cat2).c_str());
   cout<<"Hello12: "<<nEvents2<<endl;
   nEventsInDataRange2->setConstant(true);


   //get Prc

   RooKeysPdf* charmShape = (RooKeysPdf*)workspaceCharmPrc->pdf(("templateCharmPrc"+trigCatPrc).c_str());
   RooRealVar* nMCCharmPrc = (RooRealVar*)workspaceCharmPrc->var( ("nEventsInDataRangeCharmPrc"+trigCatPrc).c_str() );
   RooKeysPdf* rareShape = (RooKeysPdf*)workspaceRarePrc->pdf(("templateRarePrc"+trigCatPrc).c_str());
   RooRealVar* nMCRarePrc = (RooRealVar*)workspaceRarePrc->var( ("nEventsInDataRangeRarePrc"+trigCatPrc).c_str() );
   cout<<"Hello12: "<<charmShape<<" "<<nMCCharmPrc<<" "<<rareShape<<" "<<nMCRarePrc<<endl;
   nMCCharmPrc->setConstant(true);
   nMCRarePrc->setConstant(true);

   //get rare PDF

   //RooChebychev* rarePdf = 0;//(RooChebychev*)workspaceRare->pdf(("poly"+trigCatPrc).c_str());

   //get Piee MC fit variables
   

   RooRealVar* meanPiee = (RooRealVar*)workspacePiee->var(("meanPiee"+trigCatPrc).c_str());
   meanPiee->setConstant(true);
   RooRealVar* sigmaPiee = (RooRealVar*)workspacePiee->var(("sigmaPiee"+trigCatPrc).c_str());
   sigmaPiee->setConstant(true);
   RooRealVar* alPiee = (RooRealVar*)workspacePiee->var(("alpha_leftPiee"+trigCatPrc).c_str());
   alPiee->setConstant(true);
   RooRealVar* arPiee = (RooRealVar*)workspacePiee->var(("alpha_rightPiee"+trigCatPrc).c_str());
   arPiee->setConstant(true);
   RooRealVar* nlPiee = (RooRealVar*)workspacePiee->var(("n_leftPiee"+trigCatPrc).c_str());
   nlPiee->setConstant(true);
   RooRealVar* nrPiee = (RooRealVar*)workspacePiee->var(("n_rightPiee"+trigCatPrc).c_str());
   nrPiee->setConstant(true);
   RooRealVar* fracCB2Piee = (RooRealVar*)workspacePiee->var(("fracCB2Piee"+trigCatPrc).c_str());
   fracCB2Piee->setConstant(true);


   //get Lb MC fit variables
   

   RooKeysPdf* LbShape = (RooKeysPdf*)workspaceLb->pdf(("templateLb"+trigCatS).c_str());
   RooRealVar* nMCLb = (RooRealVar*)workspaceLb->var( ("nEventsInDataRangeLb"+trigCatS).c_str() );
   cout<<"Hello22: "<<LbShape<<" "<<nMCLb<<endl;
   nMCLb->setConstant(true);

   // RooRealVar* meanLb = (RooRealVar*)workspaceLb->var(("meanLb"+trigCatS).c_str());
   // meanLb->setConstant(true);
   // RooRealVar* sigmaLb = (RooRealVar*)workspaceLb->var(("sigmaLb"+trigCatS).c_str());
   // sigmaLb->setConstant(true);
   // RooRealVar* alLb = (RooRealVar*)workspaceLb->var(("alpha_leftLb"+trigCatS).c_str());
   // alLb->setConstant(true);
   // RooRealVar* arLb = (RooRealVar*)workspaceLb->var(("alpha_rightLb"+trigCatS).c_str());
   // arLb->setConstant(true);
   // RooRealVar* nlLb = (RooRealVar*)workspaceLb->var(("n_leftLb"+trigCatS).c_str());
   // nlLb->setConstant(true);
   // RooRealVar* nrLb = (RooRealVar*)workspaceLb->var(("n_rightLb"+trigCatS).c_str());
   // nrLb->setConstant(true);
   // RooRealVar* fracCB2Lb = (RooRealVar*)workspaceLb->var(("fracCB2Lb"+trigCatS).c_str());
   // // fracCB2Lb->setConstant(true);


   //Construct the signal model
   // RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+trigCatS).c_str(), ("sigmaScaleFactor"+trigCatS).c_str(), 1.12, 1, 1.5);
   // RooRealVar meanShift(("meanShift"+trigCatS).c_str(), ("meanShift"+trigCatS).c_str(), 1.5, -5, 5);
   RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+trigCatS).c_str(), ("sigmaScaleFactor"+trigCatS).c_str(), 1.12, 0.9, 1.5);
   RooRealVar meanShift(("meanShift"+trigCatS).c_str(), ("meanShift"+trigCatS).c_str(), 1.5, -20, 20);

   if(scaleFix > -1000)
   {
      sigmaScaleFactor.setVal(scaleFix);
      sigmaScaleFactor.setConstant(true);
   }

   if(shiftFix > -1000)
   {
      meanShift.setVal(shiftFix);
      meanShift.setConstant(true);
   }

   RooFormulaVar shiftedMean0(("shiftedMean"+cat0).c_str(), ("shiftedMean"+cat0).c_str(), "@0+@1", RooArgList(*mean0, meanShift));
   RooFormulaVar shiftedMean1(("shiftedMean"+cat1).c_str(), ("shiftedMean"+cat1).c_str(), "@0+@1", RooArgList(*mean1, meanShift));
   RooFormulaVar shiftedMean2(("shiftedMean"+cat2).c_str(), ("shiftedMean"+cat2).c_str(), "@0+@1", RooArgList(*mean2, meanShift));
   
   RooFormulaVar shiftedMean20(("shiftedMean2"+cat0).c_str(), ("shiftedMean2"+cat0).c_str(), "@0+@1", RooArgList(*mean20, meanShift));
   RooFormulaVar shiftedMean21(("shiftedMean2"+cat1).c_str(), ("shiftedMean2"+cat1).c_str(), "@0+@1", RooArgList(*mean21, meanShift));
   RooFormulaVar shiftedMean22(("shiftedMean2"+cat2).c_str(), ("shiftedMean2"+cat2).c_str(), "@0+@1", RooArgList(*mean22, meanShift));

   RooFormulaVar shiftedMean30(("shiftedMean3"+cat0).c_str(), ("shiftedMean3"+cat0).c_str(), "@0+@1", RooArgList(*mean30, meanShift));
   RooFormulaVar shiftedMean31(("shiftedMean3"+cat1).c_str(), ("shiftedMean3"+cat1).c_str(), "@0+@1", RooArgList(*mean31, meanShift));
   RooFormulaVar shiftedMean32(("shiftedMean3"+cat2).c_str(), ("shiftedMean3"+cat2).c_str(), "@0+@1", RooArgList(*mean32, meanShift));

   RooFormulaVar shiftedMean40(("shiftedMean4"+cat0).c_str(), ("shiftedMean4"+cat0).c_str(), "@0+@1", RooArgList(*mean40, meanShift));
   RooFormulaVar shiftedMean41(("shiftedMean4"+cat1).c_str(), ("shiftedMean4"+cat1).c_str(), "@0+@1", RooArgList(*mean41, meanShift));
   RooFormulaVar shiftedMean42(("shiftedMean4"+cat2).c_str(), ("shiftedMean4"+cat2).c_str(), "@0+@1", RooArgList(*mean42, meanShift));

   RooFormulaVar scaledSigma0(("scaledSigma"+cat0).c_str(), ("scaledSigma"+cat0).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma0));
   RooFormulaVar scaledSigma1(("scaledSigma"+cat1).c_str(), ("scaledSigma"+cat1).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma1));
   RooFormulaVar scaledSigma2(("scaledSigma"+cat2).c_str(), ("scaledSigma"+cat2).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma2));

   RooFormulaVar scaledSigma20(("scaledSigma2"+cat0).c_str(), ("scaledSigma2"+cat0).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma20));
   RooFormulaVar scaledSigma21(("scaledSigma2"+cat1).c_str(), ("scaledSigma2"+cat1).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma21));
   RooFormulaVar scaledSigma22(("scaledSigma2"+cat2).c_str(), ("scaledSigma2"+cat2).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma22));

   RooFormulaVar scaledSigma30(("scaledSigma3"+cat0).c_str(), ("scaledSigma3"+cat0).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma30));
   RooFormulaVar scaledSigma31(("scaledSigma3"+cat1).c_str(), ("scaledSigma3"+cat1).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma31));
   RooFormulaVar scaledSigma32(("scaledSigma3"+cat2).c_str(), ("scaledSigma3"+cat2).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma32));

   RooFormulaVar scaledSigma40(("scaledSigma4"+cat0).c_str(), ("scaledSigma4"+cat0).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma40));
   RooFormulaVar scaledSigma41(("scaledSigma4"+cat1).c_str(), ("scaledSigma4"+cat1).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma41));
   RooFormulaVar scaledSigma42(("scaledSigma4"+cat2).c_str(), ("scaledSigma4"+cat2).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma42));

   
   RooRealVar arScale("arScale", "arScale", 1,0,3);
   arScale.setConstant(true);
   RooFormulaVar scaledAr0(("scaledAr"+cat0).c_str(), ("scaledAr"+cat0).c_str(), "@0*@1", RooArgList(*ar0, arScale));
   RooFormulaVar scaledAr1(("scaledAr"+cat1).c_str(), ("scaledAr"+cat1).c_str(), "@0*@1", RooArgList(*ar1, arScale));
   RooFormulaVar scaledAr2(("scaledAr"+cat2).c_str(), ("scaledAr"+cat2).c_str(), "@0*@1", RooArgList(*ar2, arScale));


   RooCBShape cb10(("cb1"+cat0).c_str(), ("cb1"+cat0).c_str(), *B_plus_M, shiftedMean0, scaledSigma0, *al0, *nl0);
   RooCBShape cb20(("cb2"+cat0).c_str(), ("cb2"+cat0).c_str(), *B_plus_M, shiftedMean20, scaledSigma20, scaledAr0, *nr0);
   RooGaussian gaus10(("gaus1"+cat0).c_str(), ("gaus1"+cat0).c_str(), *B_plus_M, shiftedMean30, scaledSigma30);
   

   RooCBShape cb11(("cb1"+cat1).c_str(), ("cb1"+cat1).c_str(), *B_plus_M, shiftedMean1, scaledSigma1, *al1, *nl1);
   RooCBShape cb21(("cb2"+cat1).c_str(), ("cb2"+cat1).c_str(), *B_plus_M, shiftedMean21, scaledSigma21, scaledAr1, *nr1);
   RooGaussian gaus11(("gaus1"+cat1).c_str(), ("gaus1"+cat1).c_str(), *B_plus_M, shiftedMean31, scaledSigma31);
   RooGaussian gaus21(("gaus2"+cat1).c_str(), ("gaus2"+cat1).c_str(), *B_plus_M, shiftedMean41, scaledSigma41);


   RooCBShape cb12(("cb1"+cat2).c_str(), ("cb1"+cat2).c_str(), *B_plus_M, shiftedMean2, scaledSigma2, *al2, *nl2);
   RooCBShape cb22(("cb2"+cat2).c_str(), ("cb2"+cat2).c_str(), *B_plus_M, shiftedMean22, scaledSigma22, scaledAr2, *nr2);
   RooGaussian gaus12(("gaus1"+cat2).c_str(), ("gaus1"+cat2).c_str(), *B_plus_M, shiftedMean32, scaledSigma32);
   RooGaussian gaus22(("gaus2"+cat2).c_str(), ("gaus2"+cat2).c_str(), *B_plus_M, shiftedMean42, scaledSigma42);




   RooAbsPdf *modelSig0;
   RooAbsPdf *modelSig1;
   RooAbsPdf *modelSig2;
   
   if (trigCat==0)
   {
     modelSig0 = new RooAddPdf(("modelSig"+cat0).c_str(), ("modelSig"+cat0).c_str(), RooArgList(cb20, cb10 ) ,RooArgList(*fracCB20));
     // RooCBShape modelSig0(("modelSig"+cat0).c_str(), ("modelSig"+cat0).c_str(), *B_plus_M, shiftedMean0, scaledSigma0, *al0, *nl0);
     modelSig1 = new RooAddPdf(("modelSig"+cat1).c_str(), ("modelSig"+cat1).c_str(), RooArgList(cb21, gaus11, gaus21, cb11 ) ,
                         RooArgList(*fracCB21, *fracGaus1, *fracGaus21),true);
     modelSig2 = new RooAddPdf(("modelSig"+cat2).c_str(), ("modelSig"+cat2).c_str(), RooArgList(cb22, gaus12, gaus22, cb12 ) ,
                       RooArgList(*fracCB22, *fracGaus2, *fracGaus22),true);
   }
   if (trigCat==1)
   {
     // modelSig0 = new RooAddPdf(("modelSig"+cat0).c_str(), ("modelSig"+cat0).c_str(), RooArgList(cb10, gaus10 ) ,RooArgList(*fracCB20));
     modelSig0 = new RooCBShape(("modelSig"+cat0).c_str(), ("modelSig"+cat0).c_str(), *B_plus_M, shiftedMean0, scaledSigma0, *al0, *nl0);
     modelSig1 = new RooAddPdf(("modelSig"+cat1).c_str(), ("modelSig"+cat1).c_str(), RooArgList(cb21, cb11 ) ,
                         RooArgList(*fracCB21),true);
     modelSig2 = new RooAddPdf(("modelSig"+cat2).c_str(), ("modelSig"+cat2).c_str(), RooArgList(cb22, gaus12, cb12 ) ,
                       RooArgList(*fracCB22, *fracGaus2),true);     
   }
   if (trigCat==2)
   {
     // modelSig0 = new RooAddPdf(("modelSig"+cat0).c_str(), ("modelSig"+cat0).c_str(), RooArgList(cb10, gaus10 ) ,RooArgList(*fracCB20));
     modelSig0 = new RooCBShape(("modelSig"+cat0).c_str(), ("modelSig"+cat0).c_str(), *B_plus_M, shiftedMean0, scaledSigma0, *al0, *nl0);
     modelSig1 = new RooAddPdf(("modelSig"+cat1).c_str(), ("modelSig"+cat1).c_str(), RooArgList(cb11, gaus11 ) ,
                         RooArgList(*fracCB21),true);
     modelSig2 = new RooAddPdf(("modelSig"+cat2).c_str(), ("modelSig"+cat2).c_str(), RooArgList(cb22, gaus12, cb12 ) ,
                       RooArgList(*fracCB22, *fracGaus2),true);     
   }
   
   RooFormulaVar frac0gammaMu(("frac0gammaMu"+trigCatS).c_str(), ("frac0gammaMu"+trigCatS).c_str(), "@0/(@0+@1+@2)", 
                              RooArgList(*nEventsInDataRange0, *nEventsInDataRange1, *nEventsInDataRange2));
   // RooFormulaVar frac0gammaSigma(("frac0gammaSigma"+trigCatS).c_str(), ("frac0gammaSigma"+trigCatS).c_str(), "0.01*@0/(@0+@1+@2)", RooArgList(*nEventsInDataRange0, *nEventsInDataRange1, *nEventsInDataRange2));
   RooRealVar frac0gammaSigma(("frac0gammaSigma"+trigCatS).c_str(), ("frac0gammaSigma"+trigCatS).c_str(), 0.001);
   RooFormulaVar frac1gammaMu(("frac1gammaMu"+trigCatS).c_str(), ("frac1gammaMu"+trigCatS).c_str(), "@1/(@0+@1+@2)", 
                              RooArgList(*nEventsInDataRange0, *nEventsInDataRange1, *nEventsInDataRange2));
   // RooFormulaVar frac1gammaSigma(("frac1gammaSigma"+trigCatS).c_str(), ("frac1gammaSigma"+trigCatS).c_str(), "0.01*@1/(@0+@1+@2)", RooArgList(*nEventsInDataRange0, *nEventsInDataRange1, *nEventsInDataRange2));
   RooRealVar frac1gammaSigma(("frac1gammaSigma"+trigCatS).c_str(), ("frac1gammaSigma"+trigCatS).c_str(), 0.001);
   cout<<"Frac0 val from MC: "<<frac0gammaMu.getVal()<<endl;
   cout<<"Frac1 val from MC: "<<frac1gammaMu.getVal()<<endl;
   


   RooRealVar frac0gamma(("frac0gamma"+trigCatS).c_str(), ("frac0gamma"+trigCatS).c_str(), frac0gammaMu.getVal(), 
                         frac0gammaMu.getVal()-10*frac0gammaSigma.getVal(), frac0gammaMu.getVal()+10*frac0gammaSigma.getVal());
   RooRealVar frac1gamma(("frac1gamma"+trigCatS).c_str(), ("frac1gamma"+trigCatS).c_str(), frac1gammaMu.getVal(), 
                         frac1gammaMu.getVal()-10*frac1gammaSigma.getVal(), frac1gammaMu.getVal()+10*frac1gammaSigma.getVal());

   RooGaussian frac0gammaConstraint(("frac0gammaConstraint"+trigCatS).c_str(), ("frac0gammaConstraint"+trigCatS).c_str(), 
                                    frac0gamma, frac0gammaMu, frac0gammaSigma);
   RooGaussian frac1gammaConstraint(("frac1gammaConstraint"+trigCatS).c_str(), ("frac1gammaConstraint"+trigCatS).c_str(), 
                                    frac1gamma, frac1gammaMu, frac1gammaSigma);

   RooAddPdf modelSigTot(("modelSigTot"+trigCatS).c_str(), ("modelSigTot"+trigCatS).c_str(), 
                         RooArgList(*modelSig0, *modelSig1, *modelSig2), RooArgList(frac0gamma, frac1gamma));

   
   meanShift.setVal(0.);
   sigmaScaleFactor.setVal(1.);
   plotOnMC(&modelSigTot,B_plus_M,trigCat);
   
   
   if(fast) frac0gamma.setConstant(true);
   if(fast) frac1gamma.setConstant(true);

   if(nEventsInDataRange0->getVal() < 1.) 
   {
      frac0gamma.setVal(0.);
      frac0gamma.setConstant(true);
   }

   if(nEventsInDataRange1->getVal() < 1.) 
   {
      frac1gamma.setVal(0.);
      frac1gamma.setConstant(true);
   }

   //construct the bkg model
   
   RooRealVar lambda(("lambda"+trigCatS).c_str(),("lambda"+trigCatS).c_str(),-5e-04,-0.01, 0.01);
   RooExponential expo(("modelBkgTot"+trigCatS).c_str(), ("modelBkgTot"+trigCatS).c_str(), *B_plus_M,  lambda);
   // // RooRealVar lambda(("lambda"+trigCatS).c_str(),("lambda"+trigCatS).c_str(),-1.,1.);
   // // RooChebychev comb(("modelComb"+trigCatS).c_str(), ("modelComb"+trigCatS).c_str(), *B_plus_M, RooArgList(lambda));
   // RooExponential comb(("modelComb"+trigCatS).c_str(), ("modelComb"+trigCatS).c_str(), *B_plus_M,  lambda);
   // // RooRealVar cheby1(("cheby1"+trigCatS).c_str(),("cheby1"+trigCatS).c_str(),-1.,1);
   // // RooRealVar cheby2(("cheby2"+trigCatS).c_str(),("cheby2"+trigCatS).c_str(),-1.,1.);
   // // RooChebychev Lb(("modelLb"+trigCatS).c_str(), ("modelLb"+trigCatS).c_str(), *B_plus_M, RooArgList(cheby1,cheby2));
   // RooRealVar cheby1(("cheby1"+trigCatS).c_str(),("cheby1"+trigCatS).c_str(),5200,5500);
   // RooRealVar cheby2(("cheby2"+trigCatS).c_str(),("cheby2"+trigCatS).c_str(),50,2000);
   // RooGaussian Lb(("modelLb"+trigCatS).c_str(), ("modelLb"+trigCatS).c_str(), *B_plus_M, cheby1,cheby2);
   // RooRealVar fracLbcomb(("fracLbcomb"+trigCatS).c_str(),("fracLbcomb"+trigCatS).c_str(),0,1);

   // RooAddPdf expo(("modelBkgTot"+trigCatS).c_str(), ("modelBkgTot"+trigCatS).c_str(), Lb, comb, fracLbcomb);

   // RooRealVar bkg1(("bkg1"+trigCatS).c_str(),("bkg1"+trigCatS).c_str(),0.,1.);
   // RooRealVar bkg2(("bkg2"+trigCatS).c_str(),("bkg2"+trigCatS).c_str(),0.,1.);
   // RooRealVar bkg3(("bkg3"+trigCatS).c_str(),("bkg3"+trigCatS).c_str(),0.,1.);
   
   // RooChebychev expo(("modelBkgTot"+trigCatS).c_str(), ("modelBkgTot"+trigCatS).c_str(), *B_plus_M, RooArgList(bkg1, bkg2, bkg3));
   

   //Construct the Piee model

   RooFormulaVar shiftedMeanPiee(("shiftedMeanPiee"+trigCatS).c_str(), ("shiftedMeanPiee"+trigCatS).c_str(), "@0+@1", RooArgList(*meanPiee, meanShift));
   RooFormulaVar scaledSigmaPiee(("scaledSigmaPiee"+trigCatS).c_str(), ("scaledSigmaPiee"+trigCatS).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigmaPiee));
   
   RooFormulaVar scaledArPiee(("scaledArPiee"+trigCatS).c_str(), ("scaledArPiee"+trigCatS).c_str(), "@0*@1", RooArgList(*arPiee, arScale));


   RooCBShape cb1Piee(("cb1Piee"+trigCatS).c_str(), ("cb1Piee"+trigCatS).c_str(), *B_plus_M, shiftedMeanPiee, scaledSigmaPiee, *alPiee, *nlPiee);
   RooCBShape cb2Piee(("cb2Piee"+trigCatS).c_str(), ("cb2Piee"+trigCatS).c_str(), *B_plus_M, shiftedMeanPiee, scaledSigmaPiee, scaledArPiee, *nrPiee);
   
   RooAddPdf modelPieeTot(("modelPieeTot"+trigCatS).c_str(), ("modelPieeTot"+trigCatS).c_str(), RooArgList(cb1Piee, cb2Piee ) ,RooArgList(*fracCB2Piee));


   //Construct the Lb model


   // RooFormulaVar shiftedMeanLb(("shiftedMeanLb"+trigCatS).c_str(), ("shiftedMeanLb"+trigCatS).c_str(), "@0+@1", RooArgList(*meanLb, meanShift));
   // RooFormulaVar scaledSigmaLb(("scaledSigmaLb"+trigCatS).c_str(), ("scaledSigmaLb"+trigCatS).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigmaLb));
   
   // RooFormulaVar scaledArLb(("scaledArLb"+trigCatS).c_str(), ("scaledArLb"+trigCatS).c_str(), "@0*@1", RooArgList(*arLb, arScale));


   // RooCBShape cb1Lb(("cb1Lb"+trigCatS).c_str(), ("cb1Lb"+trigCatS).c_str(), *B_plus_M, shiftedMeanLb, scaledSigmaLb, *alLb, *nlLb);
   // RooCBShape cb2Lb(("cb2Lb"+trigCatS).c_str(), ("cb2Lb"+trigCatS).c_str(), *B_plus_M, shiftedMeanLb, scaledSigmaLb, scaledArLb, *nrLb);

   // RooAddPdf modelLbTot(("modelLbTot"+trigCatS).c_str(), ("modelLbTot"+trigCatS).c_str(), RooArgList(cb1Lb, cb2Lb ) ,RooArgList(*fracCB2Lb));

   



   //construct the PRC model

   RooRealVar fracCharmPrc(("fracCharmPrc"+trigCatS).c_str(), ("fracCharmPrc"+trigCatS).c_str(), nMCCharmPrc->getVal()/(nMCCharmPrc->getVal()+nMCRarePrc->getVal() ),   0,1);
   // fracCharmPrc.setVal(nMCCharmPrc->getVal()/(nMCCharmPrc->getVal()+nMCRarePrc->getVal()+0.25*nMCLb->getVal()));
   // fracCharmPrc.setVal(0.);
   fracCharmPrc.setConstant(true);

   RooRealVar fracLbPrc(("fracLbPrc"+trigCatS).c_str(), ("fracLbPrc"+trigCatS).c_str(), 0.25*nMCLb->getVal()/(0.25*nMCLb->getVal()+nMCRarePrc->getVal()+nMCCharmPrc->getVal() ),   0,1);
   fracLbPrc.setVal(0.25*nMCLb->getVal()/(0.25*nMCLb->getVal()+nMCRarePrc->getVal() ));
   

   // RooAddPdf modelPrcTot(("modelPrcTot"+trigCatS).c_str(), ("modelPrcTot"+trigCatS).c_str(), RooArgList(*LbShape, *rareShape), RooArgList(fracLbPrc) );
   // RooAddPdf modelPrcTot(("modelPrcTot"+trigCatS).c_str(), ("modelPrcTot"+trigCatS).c_str(), RooArgList(*charmShape, *LbShape, *rareShape), RooArgList(fracCharmPrc, fracLbPrc) );

   RooAbsPdf *modelPrcTot;
   
   if (nMCCharmPrc->getVal())
   {
   modelPrcTot = new RooAddPdf(("modelPrcTot"+trigCatS).c_str(), ("modelPrcTot"+trigCatS).c_str(), RooArgList(*charmShape, *rareShape), RooArgList(fracCharmPrc) );
   }
   else
   {
     modelPrcTot = rareShape;
   }
   
   

   cout<<endl<<endl<<"MEUH MEUH MEUH "<<fracLbPrc.getVal()<<" "<<data->sumEntries();
   cout<<endl<<endl<<"MEUH MEUH MEUH "<<fracCharmPrc.getVal()<<" "<<data->sumEntries()<<endl<<endl;

   //prepare initial values for yields

   double ntot0(data->sumEntries());
   double nsig0( 0.87*ntot0);
   double nprc0( 0.09*ntot0);
   double nbkg0( 0.04*ntot0);

   //construct the model total
   RooRealVar sig(("sig"+trigCatSData).c_str(), ("sig"+trigCatSData).c_str(),nsig0, 0, 2000000);
   RooRealVar bkg(("bkg"+trigCatSData).c_str(), ("bkg"+trigCatSData).c_str(),nbkg0, 0, 800000);
   RooRealVar prc(("prc"+trigCatSData).c_str(), ("prc"+trigCatSData).c_str(),nprc0, 0, 20000);
   RooRealVar fracPiee(("fracPiee"+trigCatSData).c_str(), ("fracPiee"+trigCatSData).c_str(), 0,0,1);
   RooFormulaVar piee(("piee"+trigCatSData).c_str(), ("piee"+trigCatSData).c_str(), "@0*@1", RooArgList(sig, fracPiee) ); 
   RooRealVar fracLb(("fracLb"+trigCatSData).c_str(), ("fracLb"+trigCatSData).c_str(), 0,0,0.3);
   RooFormulaVar lb(("lb"+trigCatSData).c_str(), ("lb"+trigCatSData).c_str(), "@0*@1", RooArgList(sig, fracLb) ); 

   RooRealVar fracRare(("fracRare"+trigCatSData).c_str(), ("fracRare"+trigCatSData).c_str(), 0,0,1);
   RooFormulaVar rare(("rare"+trigCatSData).c_str(), ("rare"+trigCatSData).c_str(), "@0*@1", RooArgList(sig, fracRare) ); 

   // RooArgList models( modelSigTot, expo, *modelPrcTot,*LbShape);
   // RooArgList yields( sig, bkg, prc,lb);
   RooArgList models( modelSigTot, expo, *modelPrcTot);
   RooArgList yields( sig, bkg, prc);

   if(!fast)
   {
      models.add(modelPieeTot);
      // models.add(modelLbTot);
      //models.add( *rarePdf);

      yields.add(piee);
      // yields.add(lb);
      //yields.add(rare);
   }

  // RooArgList models( modelSigTot, expo);
  // RooArgList yields( sig, bkg);

   RooAddPdf modelTot(("modelTot"+i2s(trigCat)).c_str(),("modelTot"+i2s(trigCat)).c_str(), models, yields);


   //prepare piee frac constraint

   RooRealVar fracPieeMu("fracPieeMu", "fracPieeMu", 0.5);
   if(!fast) getFracPiee(trigCat, fracPieeMu);
   fracPieeMu.SetName(("fracPieeMu"+trigCatS).c_str()); fracPieeMu.SetTitle(("fracPieeMu"+trigCatS).c_str());
   fracPieeMu.setConstant(true);
   RooRealVar fracPieeSigma( ("fracPieeSigma"+trigCatS).c_str(), ("fracPieeSigma"+trigCatS).c_str(),  fracPieeMu.getError());
   // RooRealVar fracPieeSigma( ("fracPieeSigma"+trigCatS).c_str(), ("fracPieeSigma"+trigCatS).c_str(),  0.02);
   cout<<"Error in fracPieeMu: "<<fracPieeMu.getError()<<endl;
   
   // fracPieeSigma.setConstant(true);

   RooGaussian fracPieeConstraint( ("fracPieeConstraint"+trigCatS).c_str(), ("fracPieeConstraint"+trigCatS).c_str(), fracPiee, fracPieeMu, fracPieeSigma);

   fracPiee.setVal(fracPieeMu.getVal());
   fracPiee.setRange(fracPieeMu.getVal() - 10*fracPieeSigma.getVal(), fracPieeMu.getVal() + 10*fracPieeSigma.getVal());
   
   // fracPiee.setConstant(true);
   


   //prepare the rare constraint

   RooRealVar fracRareMu("fracRareMu", "fracRareMu", 0., 1);
   fracRareMu.SetName(("fracRareMu"+trigCatS).c_str()); fracRareMu.SetTitle(("fracRareMu"+trigCatS).c_str());
   fracRareMu.setConstant(true);
   RooRealVar fracRareSigma( ("fracRareSigma"+trigCatS).c_str(), ("fracRareSigma"+trigCatS).c_str(),  fracRareMu.getError());
   fracRareSigma.setConstant(true);

   RooGaussian fracRareConstraint( ("fracRareConstraint"+trigCatS).c_str(), ("fracRareConstraint"+trigCatS).c_str(), fracRare, fracRareMu, fracRareSigma);

   fracRare.setVal(fracRareMu.getVal());
   fracRare.setRange(0, fracRareMu.getVal() + 20*fracRareSigma.getVal());


   //fit

   RooAbsReal* nll(NULL);

   if(!fast) 
   {
      if(nEventsInDataRange0->getVal() < 1.) nll = modelTot.createNLL(*data, RooFit::NumCPU(8), ExternalConstraints(RooArgSet(fracPieeConstraint, frac1gammaConstraint)));
      if(nEventsInDataRange1->getVal() < 1.) nll = modelTot.createNLL(*data, RooFit::NumCPU(8), ExternalConstraints(RooArgSet(fracPieeConstraint, frac0gammaConstraint)));
      if(nEventsInDataRange0->getVal() >= 1. && nEventsInDataRange1->getVal() >= 1.) nll = modelTot.createNLL(*data, RooFit::NumCPU(8), ExternalConstraints(RooArgSet(fracPieeConstraint, frac0gammaConstraint, frac1gammaConstraint)));
   }


   if(fast) nll = modelTot.createNLL(*data, RooFit::NumCPU(8));

   RooMinuit m(*nll);
   m.setStrategy(2);
   int migradStatus, hesseStatus;



   migradStatus = m.migrad();

   if(migradStatus == 4)
   {
      TRandom rand;
      cout<<endl<<"WARNING: fit has failed, going to retry..."<<endl;
      for(int i(0); i<4 && migradStatus != 0; ++i)
      {
         cout<<"Refit "<<i+1<<endl;
         sig.setVal(rand.Uniform(ntot0));
         prc.setVal(rand.Uniform(ntot0));
         bkg.setVal(rand.Uniform(ntot0));
         migradStatus = m.migrad();
      }
   }

   hesseStatus = m.hesse();
   if(!fast) m.minos();



   //plot

   B_plus_M->setBins(50);
   B_plus_M->SetTitle("M(Kee)");
   RooPlot* frame = B_plus_M->frame();
   data->plotOn(frame);
   modelTot.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");


   ofstream out((plotdir+"plotKeeData"+trigCatSData+".dat").c_str());
   saveFitInfo(out, frame, 6, &modelTot);

   modelTot.plotOn(frame, RooFit::LineColor(kGreen), RooFit::Components( ("modelBkgTot"+trigCatS).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kBlue), RooFit::Components( ("modelSigTot"+trigCatS).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kMagenta), RooFit::Components(("templateCharmPrc"+trigCatPrc).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kOrange), RooFit::Components( ("templateRarePrc"+trigCatPrc).c_str()));
   modelTot.plotOn(frame, RooFit::LineColor(kCyan), RooFit::Components( ("templateLb"+trigCatS).c_str()) );
   if (!fast) modelTot.plotOn(frame, RooFit::LineColor(kGray+2), RooFit::Components( ("modelPieeTot"+trigCatS).c_str() ));
   if (!fast) modelTot.plotOn(frame, RooFit::LineColor(kCyan+2), RooFit::Components( ("modelLbTot"+trigCatS).c_str() ));
   modelTot.plotOn(frame, RooFit::LineColor(kViolet+5), RooFit::Components( ("poly"+trigCatPrc).c_str() ));

   TCanvas canv("canv", "canv", 600, 600);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotKeeData"+trigCatSData+".pdf").c_str());
   canvTot.Print((plotdir+"plotKeeData"+trigCatSData+".root").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotKeeDataLogy"+trigCatSData+".pdf").c_str());
   canvTot.Print((plotdir+"plotKeeDataLogy"+trigCatSData+".root").c_str());

   //print infos 

   RooFitResult* fitRes = m.save();
   double edm (fitRes->edm());
   out<<"Migrad status: "<<migradStatus<<endl;
   out<<"Hesse status: "<<hesseStatus<<endl;
   out<<"EDM: "<<edm<<endl;
   out.close();

   //prepare what is needed for sWeights (fix everything but yields, turn formulae into variables)

   cout<<"HELLO1"<<endl;
   RooRealVar pieeVar(("pieeVar"+trigCatSData).c_str(), ("pieeVar"+trigCatSData).c_str(), piee.getVal());
   pieeVar.setError(piee.getPropagatedError(*fitRes));

   cout<<"HELLO1.5"<<endl;
   RooRealVar lbVar(("lbVar"+trigCatSData).c_str(), ("lbVar"+trigCatSData).c_str(), lb.getVal());
   lbVar.setError(lb.getPropagatedError(*fitRes));

   cout<<"HELLO2"<<endl;
   RooRealVar rareVar(("rareVar"+trigCatSData).c_str(), ("rareVar"+trigCatSData).c_str(), rare.getVal());
   rareVar.setError(rare.getPropagatedError(*fitRes));

   cout<<"HELLO3"<<endl;
   RooAddPdf modelTotForSWeight(("modelTotForSWeight"+trigCatSData).c_str(),("modelTotForSWeight"+trigCatSData).c_str(), RooArgList( modelSigTot, expo, *modelPrcTot, modelPieeTot),
      RooArgList(sig, bkg, prc, pieeVar));

   cout<<"HELLO4"<<endl;
   sigmaScaleFactor.setConstant(true);
   meanShift.setConstant(true);

   cout<<"HELLO5"<<endl;

   //save infos





   yield.setVal(sig.getVal());
   yield.setError(sig.getError());


   if(saveYield)
   {
      TFile fHistoYield(fHistoYieldName.c_str(), "UPDATE");
      TH1D* histoYield = (TH1D*)fHistoYield.Get("histoYield");
      if(histoYield)
      {
         histoYield->SetBinContent(trigCat, sig.getVal()); 
         histoYield->SetBinError(trigCat, sig.getError()); 
         histoYield->Write("", TObject::kOverwrite);
      }
      fHistoYield.Close();
   }

   fw.cd();
   RooWorkspace ws2(("workspaceDataFit"+trigCatSData).c_str(), ("workspaceDataFit"+trigCatSData).c_str());
   //ws2.import(modelTot);
   //
   cout<<"HELLO6"<<endl;
   ws2.import(*data);
   cout<<"HELLO7"<<endl;
   ws2.import(modelTotForSWeight);
   cout<<"HELLO8"<<endl;

   ws2.Write("", TObject::kOverwrite);
   cout<<"HELLO9"<<endl;

   fw.Close();
   cout<<"HELLO10"<<endl;
   fpull.Close();
   cout<<"HELLO11"<<endl;

   if(nll != NULL) delete nll;
   cout<<"HELLO12"<<endl;
   if(fitRes != NULL) delete fitRes;
   cout<<"HELLO13"<<endl;


}

void FitAndSplotKeeDataForTraining::fitAllPrc()
{
   for(int i(0); i<10; ++i)
   {
      for(int j(-1); j<3; ++j)
      {
         fitPrc(i,j,"Charm");
         fitPrc(i,j,"Rare");
      }
   }
}

void FitAndSplotKeeDataForTraining::fitPrc(int trigCat, int nPhotons, string whichOne)
{
   cout<<"fitting Prc"<<endl;

   if( whichOne != "Rare" && whichOne != "Charm")
   {
      cerr<<"ERROR: in function fitPrc, whichOne set to "<<whichOne<<", not valid."<<endl; 
      return;
   }

   string cat("Trig"+i2s(trigCat)+"Phot"+i2s(nPhotons));

   //Get the dataset

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspacePrcData = (RooWorkspace*)fw.Get( ("workspace"+whichOne+"PrcData"+cat).c_str()  );
   if(!workspacePrcData)
   {
      cerr<<"ERROR: in function fitPrc, no workspace "<<"workspace"+whichOne+"PrcData"+cat<<" found in file "<<workspaceFileName<<endl;
      return;
   }

   RooRealVar* B_plus_M = workspacePrcData->var("B_plus_M");
   B_plus_M->setBins(50);
   RooDataSet* data = (RooDataSet*)workspacePrcData->data("data");

   if(!data || !B_plus_M)
   {
      cerr<<"ERROR: in function fitPrc, error downloading stuff from workspace"<<endl;
      return;
   }

   //prepare the fit function

   RooKeysPdf model( ("template"+whichOne+"Prc"+cat).c_str(), ("template"+whichOne+"Prc"+cat).c_str(), *B_plus_M, *data, RooKeysPdf::MirrorBoth, 1.0); 

   //plot and save infos

   ofstream out((plotdir+"plotKee"+whichOne+"Prc"+cat+".dat").c_str());

   RooPlot* frame = B_plus_M->frame();
   data->plotOn(frame);
   model.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");




   TCanvas canv("canv", "canv", 600, 600);
   frame->SetMinimum(0.5);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotKee"+whichOne+"Prc"+cat+".pdf").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotKee"+whichOne+"PrcLogy"+cat+".pdf").c_str());

   //save the fit

   RooWorkspace workspacePrcFit( ("workspace"+whichOne+"PrcFit"+cat).c_str(), ("workspace"+whichOne+"PrcFit"+cat).c_str()  );

   workspacePrcFit.import(*data);
   workspacePrcFit.import(*B_plus_M);
   workspacePrcFit.import(model);

   //add a variable containing the number of events
   RooRealVar nEvents(("nEvents"+whichOne+"Prc"+cat).c_str(), ("nEvents"+whichOne+"Prc"+cat).c_str(), data->sumEntries());
   workspacePrcFit.import(nEvents);
   out<<"nEvents: ";
   nEvents.writeToStream(out, false);
   out<<endl;

   RooRealVar nEventsInDataRange(("nEventsInDataRange"+whichOne+"Prc"+cat).c_str(), ("nEventsInDataRange"+whichOne+"Prc"+cat).c_str(), data->sumEntries(),
         data->sumEntries( ("B_plus_M > "+d2s(minBMass_data)+" && B_plus_M < "+d2s(maxBMass_data) ).c_str()) );
   workspacePrcFit.import(nEventsInDataRange);
   nEventsInDataRange.writeToStream(out, false);

   out.close();

   fw.cd();
   workspacePrcFit.Write("", TObject::kOverwrite);

   cout<<"Workspace for fit has been saved:"<<endl;
   workspacePrcFit.Print();

   fw.Close();
   fpull.Close();

   //if(workspacePrcData != NULL) delete workspacePrcData;
   //if(B_plus_DTFM_M_zero != NULL) delete B_plus_DTFM_M_zero;
   //if(data != NULL) delete data;
   //if(frame != NULL) delete frame;
   //if(cpull != NULL) delete cpull;

}


void FitAndSplotKeeDataForTraining::fitRare(int trigCat, int nPhotons)
{
   cout<<"fitting Rare mode"<<endl;

   string cat("Trig"+i2s(trigCat)+"Phot"+i2s(nPhotons));

   //Get the dataset

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspaceRareData = (RooWorkspace*)fw.Get( ("workspaceRareData"+cat).c_str()  );
   if(!workspaceRareData)
   {
      cerr<<"ERROR: in function fitRare, no workspace found in file "<<workspaceFileName<<endl;
      return;
   }

   RooRealVar* B_plus_DTFM_M_zero = workspaceRareData->var("B_plus_DTFM_M_zero");
   B_plus_DTFM_M_zero->setBins(50);
   RooDataSet* data = (RooDataSet*)workspaceRareData->data("data");

   if(!data || !B_plus_DTFM_M_zero)
   {
      cerr<<"ERROR: in function fitRare, error downloading stuff from workspace"<<endl;
      return;
   }

   //prepare the fit function

   RooRealVar c1("c1", "c1", 0, -1, 1);
   RooRealVar c2("c2", "c2", 0, -1, 1);
   RooRealVar c3("c3", "c3", 0, -1, 1);
   RooRealVar c4("c4", "c4", 0, -1, 1);


   RooChebychev model( ("poly"+cat).c_str(), ("poly"+cat).c_str(), *B_plus_DTFM_M_zero, RooArgList(c1,c2,c3,c4)); 

   //fit
   //
   RooAbsReal* nll = model.createNLL(*data);

   RooMinuit m(*nll);
   int migradStatus, hesseStatus;
   m.minos();
   hesseStatus = m.hesse();
   migradStatus = m.migrad();

   //plot and save infos

   c1.setConstant(true);
   c2.setConstant(true);
   c3.setConstant(true);
   c4.setConstant(true);

   ofstream out((plotdir+"plotKeeRare"+cat+".dat").c_str());

   RooPlot* frame = B_plus_DTFM_M_zero->frame();
   data->plotOn(frame);
   model.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");




   TCanvas canv("canv", "canv", 600, 600);
   frame->SetMinimum(0.5);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotKeeRare"+cat+".pdf").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotKeeRareLogy"+cat+".pdf").c_str());

   //save the fit

   RooWorkspace workspaceRareFit( ("workspaceRareFit"+cat).c_str(), ("workspaceRareFit"+cat).c_str()  );

   workspaceRareFit.import(*data);
   workspaceRareFit.import(*B_plus_DTFM_M_zero);
   workspaceRareFit.import(model);

   saveFitInfo(out, frame, 3, &model);
   out<<"Migrad Status: "<<migradStatus<<endl;
   out<<"Hesse Status: "<<hesseStatus<<endl;


   fw.cd();
   workspaceRareFit.Write("", TObject::kOverwrite);

   cout<<"Workspace for rare fit has been saved:"<<endl;
   workspaceRareFit.Print();

   fw.Close();
   fpull.Close();

   if(nll != NULL) delete nll;
}

void FitAndSplotKeeDataForTraining::fitMCAuto(int trigCat, int nPhotons, bool fast)
{
   if(trigCat == 0 || trigCat == 5)
   {
      if(nPhotons == 0) fitMC(trigCat, 0, 5280, 15, 15, 0.7, -0.3, 4.7, 4.28, 0.8, fast);
      if(nPhotons == 1) fitMC(trigCat, 1, 5280, 15, 15, 0.7, -0.3, 2.85, 2.58, 0.6, fast);
      if(nPhotons == 2) fitMC(trigCat, 2, 5280, 15, 11, 0.7, -0.9, 2.7, 3.1, 0.6, fast);
   }

   if(trigCat == 1 || trigCat == 4 || trigCat == 6)
   {
      if(nPhotons == 0) fitMC(trigCat, 0, 5281, 12, 1.3, 0.9, -0.5, 4.3, 13, 0.6, fast );
      if(nPhotons == 1) fitMC(trigCat, 1, 5281, 20, 1.66, 0.6, -0.11, 3.6, 2.6, 0.6, fast );
      if(nPhotons == 2) fitMC(trigCat, 2, 5280, 15, 1.2,  0.4, -0.4, 6.26, 2.9, 0.8, fast );
   }

   if(trigCat == 2)
   {
      if(nPhotons == 0) fitMC(trigCat, 0, 5280, 15, 1.2, 0.4, -0.4, 4.15, 5.3, 0.6,  fast );
      if(nPhotons == 1) fitMC(trigCat, 1, 5280, 15, 1.2, 0.4, -0.4, 3.3, 2.45,0.6,  fast );
      if(nPhotons == 2) fitMC(trigCat, 2, 5280, 15.64, 1.02, 0.5, -0.5, 4.5, 3.4, 0.8, fast );
   }

   if(trigCat == 3 || trigCat == 7)
   {
      if(nPhotons == 0) fitMC(trigCat, 0, 5280, 15, 1.2, 0.4, -0.4, 4.15, 5.3, 0.6,  fast );
      if(nPhotons == 1) fitMC(trigCat, 1, 5280, 15, 1.2, 0.4, -0.4, 3.3, 2.45,0.6,  fast );
      if(nPhotons == 2) fitMC(trigCat, 2, 5280, 15.64, 1.02, 0.5, -0.5, 4.5, 3.4, 0.8, fast );
   }
}

void FitAndSplotKeeDataForTraining::fitMC(int trigCat, int nPhotons, double mean0, double sigma0, double fracSigma0, double al0, double ar0, double nl0, double nr0, double fracCB20, bool fast )
{
   cout<<"fitting"<<endl;

   string cat("Trig"+i2s(trigCat)+"Phot"+i2s(nPhotons));

   //Get the dataset

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspaceMCData = (RooWorkspace*)fw.Get( ("workspaceMCData"+cat).c_str()  );
   if(!workspaceMCData)
   {
      cerr<<"ERROR: in function fitKeeMC, no workspace found in file "<<workspaceFileName<<endl;
      return;
   }

   RooRealVar* B_plus_M = workspaceMCData->var("B_plus_M");
   B_plus_M->setBins(50);
   RooDataSet* data = (RooDataSet*)workspaceMCData->data("data");

   if(!data || !B_plus_M)
   {
      cerr<<"ERROR: in function FitAndSplotKeeDataForTraining, error downloading stuff from workspace"<<endl;
      return;
   }

   //prepare the fit function

   //CRYSTALBALL

   RooRealVar mean(("mean"+cat).c_str() , ("mean"+cat).c_str() , mean0, 5000,5800,"MeV");
   RooRealVar mean2(("mean2"+cat).c_str() , ("mean2"+cat).c_str() , mean0, 5000,5800,"MeV");
   RooRealVar mean3(("mean3"+cat).c_str() , ("mean3"+cat).c_str() , mean0, 5000,5800,"MeV");
   RooRealVar mean4(("mean4"+cat).c_str() , ("mean4"+cat).c_str() , mean0, 5000,5800,"MeV");
   RooRealVar sigma(("sigma"+cat).c_str(), ("sigma"+cat).c_str(), sigma0, 10, 2000, "MeV");
   RooRealVar fracSigma(("fracSigma"+cat).c_str(), ("fracSigma"+cat).c_str(), fracSigma0,1,2);
   RooRealVar sigma2(("sigma2"+cat).c_str(), ("sigma2"+cat).c_str(), sigma0, 10, 2000, "MeV");
   RooRealVar al(("alpha_left"+cat).c_str(), ("alpha_left"+cat).c_str(),al0,0.01,2);
   //RooRealVar ar(("alpha_right"+cat).c_str(), ("alpha_right"+cat).c_str(), ar0,-2,-0.01);
   RooRealVar ar(("alpha_right"+cat).c_str(), ("alpha_right"+cat).c_str(), ar0,-2,-0.01);
   RooRealVar nl(("n_left"+cat).c_str(), ("n_left"+cat).c_str(),nl0,1,50);
   RooRealVar nr(("n_right"+cat).c_str(), ("n_right"+cat).c_str(), nr0,1,15);
   RooRealVar al2(("alpha_left2"+cat).c_str(), ("alpha_left2"+cat).c_str(),al0,0.01,2);


   if(fast)
   {
      nl.setConstant(true);
      nr.setConstant(true);
   }

   RooRealVar sigma3(("sigma3"+cat).c_str(), ("sigma3"+cat).c_str(), 100, 15, 2000, "MeV");
   RooRealVar sigma4(("sigma4"+cat).c_str(), ("sigma4"+cat).c_str(), 100, 15, 2000, "MeV");
   RooRealVar sigma5(("sigma5"+cat).c_str(), ("sigma5"+cat).c_str(), 100, 15, 2000, "MeV");

   RooRealVar fracGaus(("fracGaus"+cat).c_str(), ("fracGaus"+cat).c_str(), 0.01, 0, 1);
   RooRealVar fracGaus2(("fracGaus2"+cat).c_str(), ("fracGaus2"+cat).c_str(), 0.01, 0, 1);
   RooRealVar fracCB2(("fracCB2"+cat).c_str(), ("fracCB2"+cat).c_str(), fracCB20, 0.01,1);
   RooRealVar fracCB3(("fracCB3"+cat).c_str(), ("fracCB3"+cat).c_str(), fracCB20, 0.01,1);

   // nl.setMax(100);
   // nl.setVal(100);
   // nl.setConstant(1);

   RooCBShape cb1(("cb1"+cat).c_str(), ("cb1"+cat).c_str(), *B_plus_M, mean, sigma, al, nl);
   // RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, mean, sigma2, ar, nr);
   RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, mean2, sigma2, ar, nr);

   RooGaussian gaus(("gaus"+cat).c_str(), ("gaus"+cat).c_str(), *B_plus_M, mean3, sigma3 );
   RooGaussian gaus2(("gaus2"+cat).c_str(), ("gaus2"+cat).c_str(), *B_plus_M, mean4, sigma4 );


   RooAbsPdf *model;

     if (nPhotons==0)
     {
       // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, cb1 ) ,RooArgList(fracCB2), true);
       model = new RooCBShape(("model"+cat).c_str(), ("model"+cat).c_str(), *B_plus_M, mean, sigma, al, nl);
     }
     if (nPhotons==1)
     {
       model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,RooArgList(fracCB2), true);
       // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, gaus, gaus2, cb1) ,
                           // RooArgList(fracCB2,fracGaus, fracGaus2), true);
     }
     if (nPhotons==2)
     {
       // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, gaus, gaus2, cb1 ) ,
                             // RooArgList(fracCB2,fracGaus, fracGaus2), true);
       model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,
                             RooArgList(fracCB2), true);
       // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus_mean2, gaus ) ,
       // RooArgList(fracCB2,fracGaus), true);
     }

   // if (trigCat==0)
   // {
     
   //   if (nPhotons==0)
   //   {
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, cb1 ) ,RooArgList(fracCB2), true);
   //     // model = new RooCBShape(("model"+cat).c_str(), ("model"+cat).c_str(), *B_plus_M, mean, sigma, al, nl);
   //   }
   //   if (nPhotons==1)
   //   {
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,RooArgList(fracCB2), true);
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, gaus, gaus2, cb1) ,
   //                         RooArgList(fracCB2,fracGaus, fracGaus2), true);
   //   }
   //   if (nPhotons==2)
   //   {
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, gaus, gaus2, cb1 ) ,
   //                           RooArgList(fracCB2,fracGaus, fracGaus2), true);
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,
   //     //                       RooArgList(fracCB2), true);
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus_mean2, gaus ) ,
   //     // RooArgList(fracCB2,fracGaus), true);
   //   }
   // }
   // if (trigCat==1)
   // {
   //   if (nPhotons==0)
   //   {
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,RooArgList(fracCB2), true);
   //     model = new RooCBShape(("model"+cat).c_str(), ("model"+cat).c_str(), *B_plus_M, mean, sigma, al, nl);
   //   }
   //   if (nPhotons==1)
   //   {
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,RooArgList(fracCB2), true);
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, gaus, cb1) ,
   //                         RooArgList(fracCB2,fracGaus), true);
   //   }
   //   if (nPhotons==2)
   //   {
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, cb1 ) ,
   //                           RooArgList(fracCB2), true);
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,
   //     //                       RooArgList(fracCB2), true);
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus_mean2, gaus ) ,
   //     // RooArgList(fracCB2,fracGaus), true);
   //   }
   // }
   // if (trigCat==2)
   // {
   //   if (nPhotons==0)
   //   {
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,RooArgList(fracCB2), true);
   //     model = new RooCBShape(("model"+cat).c_str(), ("model"+cat).c_str(), *B_plus_M, mean, sigma, al, nl);
   //   }
   //   if (nPhotons==1)
   //   {
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,RooArgList(fracCB2), true);
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus) ,
   //                         RooArgList(fracCB2), true);
   //   }
   //   if (nPhotons==2)
   //   {
   //     model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb2, gaus, cb1 ) ,
   //                           RooArgList(fracCB2,fracGaus), true);
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus ) ,
   //     //                       RooArgList(fracCB2), true);
   //     // model = new RooAddPdf(("model"+cat).c_str(), ("model"+cat).c_str(), RooArgList(cb1, gaus_mean2, gaus ) ,
   //     // RooArgList(fracCB2,fracGaus), true);
   //   }
   // }

   // IPATIA
   //   RooRealVar a2("a2","a2", TMath::Abs(ar0), 0.1,10);
   //   RooRealVar a("a","a",TMath::Abs(al0), .01,10);
   //   RooRealVar n("n","n",nl0,0.,40);
   //   RooRealVar n2("n2","n2",nr0,.01,40);
   //   RooRealVar ipa_s("ipa_s","ipa_s", sigma0, 6, 15);
   //   RooRealVar ipa_m("ipa_m","ipa_m",mean0, 5275,5285);//#5363,5373)
   //   RooRealVar beta("beta","beta",0);//#-0.1,0.1)
   //   RooRealVar zeta("zeta","zeta",5e-03);//#, 1e-03,.1)
   //   RooRealVar l("l","l",-2,-6,6);
   //   RooIpatia2 model(("model"+cat).c_str(), ("model"+cat).c_str(), *B_plus_M,l,zeta,beta,ipa_s,ipa_m,a,n,a2,n2);

   //check if there are events to fit

   bool fitIt(true);

   if(data->sumEntries() < 10) fitIt = false;

   //fit 

   RooAbsReal* nll = model->createNLL(*data, RooFit::NumCPU(8));

   RooMinuit m(*nll);
   int migradStatus(-1000);
   int hesseStatus(-1000);

   if(fitIt)
   {
      migradStatus = m.migrad();
      hesseStatus = m.hesse();
      if(!fast) m.minos();
   }

   //plot and save infos

   ofstream out((plotdir+"plotKeeMCLogy"+cat+".dat").c_str());

   RooPlot* frame = B_plus_M->frame();
   // data->plotOn(frame, RooFit::Binning(33));
   data->plotOn(frame);
   model->plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");

   saveFitInfo(out, frame, 6, model);
   RooFitResult* fitRes = m.save();
   double edm (fitRes->edm());

   out<<"Migrad status: "<<migradStatus<<endl;
   out<<"Hesse status: "<<hesseStatus<<endl;
   out<<"EDM: "<<edm<<endl;


   TCanvas canv("canv", "canv", 600, 600);
   frame->SetMinimum(0.5);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotKeeMC"+cat+".pdf").c_str());
   canvTot.Print((plotdir+"plotKeeMC"+cat+".root").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotKeeMCLogy"+cat+".pdf").c_str());
   canvTot.Print((plotdir+"plotKeeMCLogy"+cat+".root").c_str());

   //save the fit

   RooWorkspace workspaceMCFit(("workspaceMCFit"+cat).c_str(), ("workspaceMCFit"+cat).c_str());

   workspaceMCFit.import(*data);
   workspaceMCFit.import(*B_plus_M);
   workspaceMCFit.import(mean);
   workspaceMCFit.import(mean2);
   workspaceMCFit.import(mean3);
   workspaceMCFit.import(mean4);
   workspaceMCFit.import(sigma);
   workspaceMCFit.import(al);
   workspaceMCFit.import(nl);
   workspaceMCFit.import(ar);
   workspaceMCFit.import(nr);
   workspaceMCFit.import(fracSigma);
   workspaceMCFit.import(sigma2);
   workspaceMCFit.import(sigma3);
   workspaceMCFit.import(sigma4);
   workspaceMCFit.import(fracCB2);
   workspaceMCFit.import(fracCB3);
   workspaceMCFit.import(fracGaus);
   workspaceMCFit.import(fracGaus2);

   //add a variable containing the number of events
   RooRealVar nEvents(("nEvents"+cat).c_str(), ("nEvents"+cat).c_str(), data->sumEntries());
   workspaceMCFit.import(nEvents);
   nEvents.writeToStream(out, false);

   RooRealVar nEventsInDataRange(("nEventsInDataRange"+cat).c_str(), ("nEventsInDataRange"+cat).c_str(),
         data->sumEntries( ("B_plus_M > "+d2s(minBMass_data)+" && B_plus_M < "+d2s(maxBMass_data) ).c_str()) );
   workspaceMCFit.import(nEventsInDataRange);
   nEventsInDataRange.writeToStream(out, false);

   out.close();

   fw.cd();
   workspaceMCFit.Write("", TObject::kOverwrite);

   cout<<"Workspace for fit has been saved:"<<endl;
   workspaceMCFit.Print();

   fw.Close();
   fpull.Close();

   if(nll != NULL) delete nll;
   if(fitRes != NULL) delete fitRes;
   if(model != NULL) delete model;
}





void FitAndSplotKeeDataForTraining::fitLb(int trigCat, bool fast)
{
   cout<<"fitting Lb"<<endl;

   string trigCatS("Trig"+i2s(trigCat));

   //Get the dataset

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspaceMCData = (RooWorkspace*)fw.Get( ("workspaceLbData"+trigCatS).c_str()  );
   if(!workspaceMCData)
   {
      cerr<<"ERROR: in function fitLb, no workspace found in file "<<workspaceFileName<<endl;
      return;
   }

   RooRealVar* B_plus_M = workspaceMCData->var("B_plus_M");
   B_plus_M->setBins(50);
   RooDataSet* data = (RooDataSet*)workspaceMCData->data("data");

   if(!data || !B_plus_M)
   {
      cerr<<"ERROR: in function FitLb, error downloading stuff from workspace"<<endl;
      return;
   }

   //prepare the fit function

   RooKeysPdf model( ("templateLb"+trigCatS).c_str(), ("templateLb"+trigCatS).c_str(), *B_plus_M, *data, RooKeysPdf::MirrorBoth, 1.0); 

   //fit 

   //      RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   //         RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;
   RooAbsReal* nll = model.createNLL(*data, RooFit::NumCPU(8));

   RooMinuit m(*nll);
   m.setStrategy(2);
   int migradStatus(-1000), hesseStatus(-1000);

   if(!fast)
   {
      migradStatus = m.migrad();
      hesseStatus = m.hesse();
      m.minos();
   }

   //plot and save infos

   ofstream out((plotdir+"plotLbMCLogy"+i2s(trigCat)+".dat").c_str());

   RooPlot* frame = B_plus_M->frame();
   data->plotOn(frame);
   model.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");

   saveFitInfo(out, frame, 6, &model);
   RooFitResult* fitRes = m.save();
   double edm (fitRes->edm());

   out<<"Migrad status: "<<migradStatus<<endl;
   out<<"Hesse status: "<<hesseStatus<<endl;
   out<<"EDM: "<<edm<<endl;


   TCanvas canv("canv", "canv", 600, 600);
   frame->SetMinimum(0.5);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotLbMC"+i2s(trigCat)+".pdf").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotLbMCLogy"+i2s(trigCat)+".pdf").c_str());


   //save the fit

   RooWorkspace workspaceMCFit(("workspaceLbFit"+trigCatS).c_str(), ("workspaceLbFit"+trigCatS).c_str());

   workspaceMCFit.import(*data);
   workspaceMCFit.import(*B_plus_M);
   workspaceMCFit.import(model);

   //add a variable containing the number of events
   RooRealVar nEvents(("nEventsLb"), ("nEventsLb"), data->sumEntries());
   workspaceMCFit.import(nEvents);
   nEvents.writeToStream(out, false);

   RooRealVar nEventsInDataRange(("nEventsInDataRangeLb"+trigCatS).c_str(), ("nEventsInDataRangeLb"+trigCatS).c_str(),
         data->sumEntries( ("B_plus_M > "+d2s(minBMass_data)+" && B_plus_M < "+d2s(maxBMass_data) ).c_str()) );
   workspaceMCFit.import(nEventsInDataRange);
   nEventsInDataRange.writeToStream(out, false);

   out.close();

   fw.cd();
   workspaceMCFit.Write("", TObject::kOverwrite);

   cout<<"Workspace for fit has been saved:"<<endl;
   workspaceMCFit.Print();

   fw.Close();
   fpull.Close();

   if(nll != NULL) delete nll;
   if(fitRes != NULL) delete fitRes;
}


void FitAndSplotKeeDataForTraining::fitPiee(int trigCat, int PhotCat, bool fast)
{
   cout<<"fitting Piee"<<endl;

   string cat("Trig"+i2s(trigCat)+"Phot"+i2s(PhotCat));
   string trigCatS("Trig"+i2s(trigCat));


   //Get the dataset

   TFile fw(workspaceFileName.c_str(), "UPDATE");

   RooWorkspace* workspaceMCData = (RooWorkspace*)fw.Get( ("workspacePieeData"+cat).c_str()  );
   if(!workspaceMCData)
   {
      cerr<<"ERROR: in function fitPiee, no workspace found in file "<<workspaceFileName<<endl;
      return;
   }

   RooRealVar* B_plus_M = workspaceMCData->var("B_plus_M");
   B_plus_M->setBins(50);
   RooDataSet* data = (RooDataSet*)workspaceMCData->data("data");

   if(!data || !B_plus_M)
   {
      cerr<<"ERROR: in function FitPiee, error downloading stuff from workspace"<<endl;
      return;
   }

   //prepare the fit function

   RooRealVar mean((("meanPiee")+cat).c_str() , (("meanPiee")+cat).c_str() , 5316, 5000,5500,"MeV");
   RooRealVar sigma((("sigmaPiee")+cat).c_str(), (("sigmaPiee")+cat).c_str(), 19, 10, 20000, "MeV");
   //   RooRealVar sigma2((("sigma2Piee")+cat).c_str(), (("sigma2Piee")+cat).c_str(), 1, 20, "MeV");
   RooRealVar al((("alpha_leftPiee")+cat).c_str(), (("alpha_leftPiee")+cat).c_str(),1.057, 0.01,2);
   RooRealVar ar((("alpha_rightPiee")+cat).c_str(), (("alpha_rightPiee")+cat).c_str(),-0.44, -2,-0.01);
   RooRealVar nl((("n_leftPiee")+cat).c_str(), (("n_leftPiee")+cat).c_str(), 3, 1,15);
   RooRealVar nr((("n_rightPiee")+cat).c_str(), (("n_rightPiee")+cat).c_str(), 3 ,1,15);

   //// nl.setConstant(true);
   //// nr.setConstant(true);

   RooRealVar fracCB2(("fracCB2Piee"+cat).c_str(), ("fracCB2Piee"+cat).c_str(), 0.5, 0.01,0.99);

   RooCBShape cb1(("cb1Piee"), ("cb1Piee"), *B_plus_M, mean, sigma, al, nl);
   RooCBShape cb2(("cb2Piee"), ("cb2Piee"), *B_plus_M, mean, sigma, ar, nr);

   // fracCB2.setVal(1);
   // fracCB2.setConstant(1);
   

   RooAddPdf model(("modelPiee"), ("modelPiee"), RooArgList(cb1, cb2 ) ,RooArgList(fracCB2));
   //fit 

   //      RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   //         RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;
   RooAbsReal* nll = model.createNLL(*data, RooFit::NumCPU(8));

   RooMinuit m(*nll);
   m.setStrategy(2);
   int migradStatus(-1000), hesseStatus(-1000);

   if(!fast)
   {
      migradStatus = m.migrad();
      hesseStatus = m.hesse();
      m.minos();
   }

   //plot and save infos

   ofstream out((plotdir+"plotPieeMCLogy"+cat+".dat").c_str());

   RooPlot* frame = B_plus_M->frame();
   data->plotOn(frame);
   model.plotOn(frame, RooFit::LineColor(kRed) );

   savePullPlot(*frame, plotdir+"pullPlot.root");
   TFile fpull((plotdir+"pullPlot.root").c_str());
   TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");

   saveFitInfo(out, frame, 6, &model);
   RooFitResult* fitRes = m.save();
   double edm (fitRes->edm());

   out<<"Migrad status: "<<migradStatus<<endl;
   out<<"Hesse status: "<<hesseStatus<<endl;
   out<<"EDM: "<<edm<<endl;


   TCanvas canv("canv", "canv", 600, 600);
   frame->SetMinimum(0.5);
   frame->Draw();

   TCanvas canvTot("canvTot", "canvTot", 600, 600);

   canvTot.Divide(1,2);
   canvTot.cd(1);
   gPad->SetPad(0.005, 0.205, 0.995, 0.995);
   canv.DrawClonePad();
   canvTot.cd(2);
   gPad->SetPad(0.005, 0.005, 0.995, 0.2);
   cpull->DrawClonePad();

   canvTot.Print((plotdir+"plotPieeMC"+cat+".pdf").c_str());

   canv.SetLogy();
   canvTot.cd(1);

   canv.DrawClonePad();
   canvTot.Print((plotdir+"plotPieeMCLogy"+cat+".pdf").c_str());


   //save the fit

   RooWorkspace workspaceMCFit(("workspacePieeFit"+cat).c_str(), ("workspacePieeFit"+cat).c_str());

   workspaceMCFit.import(*data);
   workspaceMCFit.import(*B_plus_M);
   workspaceMCFit.import(mean);
   workspaceMCFit.import(sigma);
   workspaceMCFit.import(al);
   workspaceMCFit.import(nl);
   workspaceMCFit.import(ar);
   workspaceMCFit.import(nr);
   //   workspaceMCFit.import(sigma2);
   workspaceMCFit.import(fracCB2);

   //add a variable containing the number of events
   RooRealVar nEvents(("nEventsPiee"), ("nEventsPiee"), data->sumEntries());
   workspaceMCFit.import(nEvents);
   nEvents.writeToStream(out, false);

   RooRealVar nEventsInDataRange(("nEventsInDataRangePiee"), ("nEventsInDataRangePiee"),
         data->sumEntries( ("B_plus_M > "+d2s(minBMass_data)+" && B_plus_M < "+d2s(maxBMass_data) ).c_str()) );
   workspaceMCFit.import(nEventsInDataRange);
   nEventsInDataRange.writeToStream(out, false);


   //add fraction wrt Kaons
   RooRealVar fracPieeKee(("fracPiee"+trigCatS).c_str(), ("fracPiee"+trigCatS).c_str(),0.5);
   getFracPiee(trigCat, fracPieeKee);
   fracPieeKee.setConstant(1);
   workspaceMCFit.import(fracPieeKee);

   out.close();

   fw.cd();
   workspaceMCFit.Write("", TObject::kOverwrite);

   cout<<"Workspace for fit has been saved:"<<endl;
   workspaceMCFit.Print();

   fw.Close();
   fpull.Close();

   if(nll != NULL) delete nll;
   if(fitRes != NULL) delete fitRes;
}


void FitAndSplotKeeDataForTraining::prepareAllPrcMCWorkspaces()
{
   for(int i(0); i<10; ++i)
   {
      for(int j(-1); j<3; ++j)
      {
         prepareRarePrcWorkspace(i,j);
         prepareCharmPrcWorkspace(i,j);
         //prepareRareWorkspace(i,j);
         preparePieeWorkspace(i,j);
         if(j>=0) prepareMCWorkspace(i,j);
      }
   }
}





//void FitAndSplotKeeDataForTraining::prepareRareWorkspace(int trigCat, int nPhotons, string forcedWeight)
//{
//   TFile f( (tupleRareDir+"/"+tupleRareName).c_str() );
//   TTree* t1 = (TTree*)f.Get( (treeRareName).c_str() );
//
//   if(!t1)
//   {
//      cerr<<"ERROR: in function prepareRareWorkspace, no tree  "<<treeRareName<<" found in "<<tupleRareDir+"/"+tupleRareName<<endl;
//      return;
//   }
//
//   t1->SetBranchStatus("*", 0);
//   t1->SetBranchStatus("B_plus_DTFM_M_zero", 1);
//   t1->SetBranchStatus("e_plus_BremMultiplicity", 1);
//   t1->SetBranchStatus("e_minus_BremMultiplicity", 1);
//   t1->SetBranchStatus("PIDTrigDataCondSPDWeight_ETOSOnly", 1);
//   t1->SetBranchStatus("PIDTrigDataCondSPDWeight_HTOSOnly", 1);
//   t1->SetBranchStatus("PIDTrigDataCondSPDWeight_TISOnly", 1);
//   t1->SetBranchStatus("PIDTrigDataCondSPDWeight_TISINCL", 1);
//   if(trigCat >= 4) t1->SetBranchStatus( ("PIDTrigDataCondSPDWeight_TrigCat"+i2s(trigCat)).c_str() );
//
//
//   setBranchStatusTTF(t1, cutString);
//
//   TFile fTrash("/vols/lhcb/th1011/toErase/toErase.root", "RECREATE");
//   TTree* t;
//   if(cutString == "" && forcedWeight == "") t = t1;
//   if(cutString != "" && forcedWeight == "")  t = t1->CopyTree(cutString.c_str());
//   if(forcedWeight != "") copyTreeWithNewVar(t, t1, cutString, forcedWeight, "forcedWeightNameForFit");
//
//
//   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M^{DTF}(B^{+})",B_plus_M_min_MC-100,B_plus_M_max_MC+100, "MeV");
//   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity", "e_plus_BremMultiplicity", -1, 3);
//   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity", "e_minus_BremMultiplicity", -1, 3);
//   B_plus_DTFM_M_zero.setBins(50);
//
//   RooRealVar weight("weight", "weight", 0, -1000, 1000);
//   RooArgSet vars(B_plus_DTFM_M_zero, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weight);// NSPDWeight);
//
//   string trigCatCut;
//   trigCatCut = "(1)"; 
//   string bremCut;
//   if(nPhotons == 0) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) == 0)";
//   if(nPhotons == 1) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) == 1)";
//   if(nPhotons == 2) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) >= 2)";
//   if(nPhotons < 0) bremCut = "(1)";
//
//
//
//   if(trigCat == 0){ weight.SetName("PIDTrigDataCondSPDWeight_ETOSOnly"); weight.SetTitle("PIDTrigDataCondSPDWeight_ETOSOnly"); }
//   if(trigCat == 1){ weight.SetName("PIDTrigDataCondSPDWeight_HTOSOnly"); weight.SetTitle("PIDTrigDataCondSPDWeight_HTOSOnly"); } 
//   if(trigCat == 2){ weight.SetName("PIDTrigDataCondSPDWeight_TISOnly"); weight.SetTitle("PIDTrigDataCondSPDWeight_TISOnly"); } 
//   if(trigCat == 3){ weight.SetName("PIDTrigDataCondSPDWeight_TISINCL"); weight.SetTitle("PIDTrigDataCondSPDWeight_TISINCL"); } 
//   if(trigCat >= 4){ weight.SetName( ("PIDTrigDataCondSPDWeight_TrigCat"+i2s(trigCat)).c_str() ); weight.SetName( ("PIDTrigDataCondSPDWeight_TrigCat"+i2s(trigCat)).c_str() );}
//   if(forcedWeight != ""){ weight.SetName("forcedWeightNameForFit"); weight.SetTitle("forcedWeightNameForFit");}
//
//   RooDataSet data("data", "data", vars, Import(*t), Cut( ("B_plus_DTFM_M_zero>"+d2s(B_plus_M_min_MC-100)+" && B_plus_DTFM_M_zero <"+d2s(B_plus_M_max_MC+100)+" && "+trigCatCut+" && "+bremCut).c_str() ),
//         WeightVar(weight) );
//
//   RooWorkspace workspaceRareData(("workspaceRareDataTrig"+i2s(trigCat)+"Phot"+i2s(nPhotons)).c_str(), ("workspaceRareData"+i2s(trigCat)+"Phot"+i2s(nPhotons)).c_str());
//
//   workspaceRareData.import(data);
//
//   TFile fws(workspaceFileName.c_str(), "UPDATE");
//   workspaceRareData.Write("",TObject::kOverwrite );
//
//   cout<<"Workspace containing rare prc data ready: "<<endl;
//
//   workspaceRareData.Print();
//   cout<<"rare data has "<<data.sumEntries()<<" events"<<endl;
//
//
//   f.Close();
//   fws.Close();
//   fTrash.Close();
//}





void FitAndSplotKeeDataForTraining::prepareCharmPrcWorkspace(int trigCat, int nPhotons)
{
   string cutStringMC2(cutStringMC); 
   string cutStringMCSave(cutStringMC);
   boost::algorithm::replace_all(cutStringMC2, "yearLabbel == 2011", "yearLabbel > 0");

   cutStringMC = cutStringMC2;
   cout<<"PREPARE Charm Prc WORKSPACE"<<endl;
   prepareGenericMCWorkspace(trigCat, nPhotons, 2);

   cutStringMC = cutStringMCSave;
}


void FitAndSplotKeeDataForTraining::preparePieeWorkspace(int trigCat, int PhotCat)
{

   string cutStringMC2(cutStringMC); 
   string cutStringMCSave(cutStringMC);
   
   cutStringMC2 = cutStringMCPiee;

   boost::algorithm::replace_all(cutStringMC2, "yearLabbel == 2011", "yearLabbel > 0");

   string weightMCSave(weightMC);
   
   weightMC = weightMCPiee;
   cutStringMC = cutStringMC2;
   cout<<"PREPARE Piee WORKSPACE"<<endl;
   // prepareGenericMCWorkspace(trigCat, -1, 3);
   prepareGenericMCWorkspace(trigCat, PhotCat, 3);

   cutStringMC = cutStringMCSave;
   weightMC = weightMCSave;
   
}

void FitAndSplotKeeDataForTraining::prepareLbWorkspace(int trigCat)
{

   string cutStringMCSave(cutStringMC);
   string cutStringMC2 = "B_plus_DTFM_M_zero>"+d2s(B_plus_DTF_M_cut);
   boost::algorithm::replace_all(cutStringMC2, "yearLabbel == 2011", "yearLabbel > 0");

   string weightMCSave(weightMC);   
   weightMC = "K_Kst_PIDWeight_eff*passTrigCat"+i2s(trigCat);

   cutStringMC = cutStringMC2;
   cout<<"PREPARE Lb WORKSPACE"<<endl;
   prepareGenericMCWorkspace(trigCat, -1, 4);

   cutStringMC = cutStringMCSave;
   weightMC = weightMCSave;
   
}


void FitAndSplotKeeDataForTraining::prepareRarePrcWorkspace(int trigCat, int nPhotons)
{
   string cutStringMC2(cutStringMC); 
   string cutStringMCSave(cutStringMC);
   boost::algorithm::replace_all(cutStringMC2, "yearLabbel == 2011", "yearLabbel > 0");

   cutStringMC = cutStringMC2;
   cout<<"PREPARE Rare Prc WORKSPACE"<<endl;
   prepareGenericMCWorkspace(trigCat, nPhotons, 1);

   cutStringMC = cutStringMCSave;
}


void FitAndSplotKeeDataForTraining::prepareMCWorkspace(int trigCat, int nPhotons)
{
   cout<<"PREPARE MC WORKSPACE"<<endl;
   string weightMC2(""); 
   string weightMCSave(weightMC);
   if (turnOffWeightsSignal) weightMC = weightMC2;
   prepareGenericMCWorkspace(trigCat, nPhotons, 0);
   weightMC = weightMCSave;

}


void FitAndSplotKeeDataForTraining::prepareGenericMCWorkspace(int trigCat, int nPhotons, int whichMC)
{
   //first select right tree to get
   string workspaceName;
   string fileName;
   string treeName;


   if(whichMC == 0) 
   {
      fileName = tupleMCdir+"/"+tupleMCname;
      treeName = treeMCname;
      workspaceName = "workspaceMCDataTrig"+i2s(trigCat)+"Phot"+i2s(nPhotons);
   }

   if(whichMC == 1)
   {
      fileName =  tupleRarePrcDir+"/"+tupleRarePrcName;
      treeName = treeRarePrcName;
      workspaceName = "workspaceRarePrcDataTrig"+i2s(trigCat)+"Phot"+i2s(nPhotons);
   }

   if(whichMC == 2)
   {
      fileName = tupleCharmPrcDir+"/"+tupleCharmPrcName;
      treeName = treeCharmPrcName;
      workspaceName = "workspaceCharmPrcDataTrig"+i2s(trigCat)+"Phot"+i2s(nPhotons);
   }

   if(whichMC == 3)
   {
      fileName = tuplePieeDir+"/"+tuplePieeName;
      treeName = treePieeName;
      workspaceName = "workspacePieeDataTrig"+i2s(trigCat)+"Phot"+i2s(nPhotons);
   }

   if(whichMC == 4)
   {
      fileName = tupleLbDir+"/"+tupleLbName;
      treeName = treeLbName;
      workspaceName = "workspaceLbDataTrig"+i2s(trigCat);
   }

   //copy tree with new vars

   vector<string> varsToSwitchOn;
   varsToSwitchOn.push_back("B_plus_DTFM_M_zero");
   varsToSwitchOn.push_back("B_plus_M");
   varsToSwitchOn.push_back("e_plus_BremMultiplicity");
   varsToSwitchOn.push_back("e_minus_BremMultiplicity");
   varsToSwitchOn.push_back("J_psi_1S_M");
   varsToSwitchOn.push_back("FidelCastro_ETOS");
   varsToSwitchOn.push_back("FidelCastro_HTOS");
   varsToSwitchOn.push_back("FidelCastro_TIS");
   varsToSwitchOn.push_back("FidelCastro_KtoPi_trigCat0");
   varsToSwitchOn.push_back("FidelCastro_KtoPi_trigCat1");
   varsToSwitchOn.push_back("FidelCastro_KtoPi_trigCat2");
   varsToSwitchOn.push_back("BDT1ETOSRun1");


   if(cutStringMC == "") cutStringMC = "(1)";
   if(weightMC == "") weightMC = "(1)";

   system("rm /vols/lhcb/palvare1/toErase/toErase.root");

   cout<<"Copying tree with new vars. Cut applied: "<<endl<<cutStringMC<<endl<<endl;

   copyTreeWithNewVars("/vols/lhcb/palvare1/toErase/toErase.root", fileName, cutStringMC, weightMC, "weightForFit", varsToSwitchOn, "DecayTree");

   //recover tree and prepare workspace


   TFile fTrash("/vols/lhcb/palvare1/toErase/toErase.root");
   TTree* t = (TTree*)fTrash.Get("DecayTree");

   if(!t)
   {
      cout<<"Tree: "<<t<<endl;
      cerr<<"ERROR: in function prepareWorkspace, no tree found in "<<tupleMCdir<<"/"<<tupleMCname<<endl;
      return;
   }


   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M^{DTF}(B^{+})",B_plus_M_min_MC,B_plus_M_max_MC, "MeV");
   RooRealVar B_plus_M("B_plus_M", "M(B^{+})",minBMass_MC,maxBMass_MC, "MeV");
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity", "e_plus_BremMultiplicity", -1, 3);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity", "e_minus_BremMultiplicity", -1, 3);
   RooRealVar J_psi_1S_M("J_psi_1S_M", "J_psi_1S_M", 0, 10000);
   RooRealVar weight("weightForFit", "weightForFit", 0, 10);


   B_plus_DTFM_M_zero.setBins(50);
   B_plus_M.setBins(50);

   RooArgSet vars(B_plus_DTFM_M_zero, B_plus_M, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weight, J_psi_1S_M);

   string bremCut;
   if(nPhotons == 0) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) == 0)";
   if(nPhotons == 1) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) == 1)";
   if(nPhotons == 2) bremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) >= 2)";
   if(nPhotons < 0) bremCut = "(1)";

   RooDataSet data("data", "data", vars, Import(*t), Cut( ("B_plus_DTFM_M_zero>"+d2s(B_plus_DTF_M_cut)+" && B_plus_M >"+d2s(minBMass_MC)+" && B_plus_M <"+d2s(maxBMass_MC)+" && "+bremCut).c_str() ),
   WeightVar( weight.GetName()  ) );


   RooWorkspace workspaceMCData( workspaceName.c_str(), workspaceName.c_str());

   workspaceMCData.import(data);

   TFile fws(workspaceFileName.c_str(), "UPDATE");
   workspaceMCData.Write("",TObject::kOverwrite );

   cout<<"Workspace containing MC data ready: "<<endl;

   workspaceMCData.Print();
   cout<<"MC data has "<<data.sumEntries()<<" events"<<endl;

   fws.Close();
   fTrash.Close();
}




void FitAndSplotKeeDataForTraining::prepareDataWorkspace(int trigCat)
{
   string trigCatS("Trig"+i2s(trigCat));

   string totFileName(tupleDataDir+"/"+tupleDataName);

   TFile f( totFileName.c_str());
   TTree* t1 = (TTree*)f.Get(treeDataName.c_str());

   if(!t1)
   {
      cerr<<"ERROR: in function prepareDataWorkspace, no tree "<<treeDataName<<" found in "<<tupleDataDir<<"/"<<tupleDataName<<endl;
      return;
   }

   t1->SetBranchStatus("*", 0);
   t1->SetBranchStatus("B_plus_M", 1);
   t1->SetBranchStatus("B_plus_DTFM_M_zero", 1);
   t1->SetBranchStatus("e_plus_BremMultiplicity", 1);
   t1->SetBranchStatus("e_minus_BremMultiplicity", 1);
   t1->SetBranchStatus("K_Kst_L0HadronDecision_TOS", 1);
   t1->SetBranchStatus("J_psi_1S_M", 1);
   t1->SetBranchStatus( ("passTrigCat"+i2s(trigCat)).c_str() );

   setBranchStatusTTF(t1, cutStringData);

   system("rm /vols/lhcb/palvare1/toErase/toErase.root");
   TFile fTrash("/vols/lhcb/palvare1/toErase/toErase.root", "RECREATE");
   TTree* t;
   t =  t1->CopyTree(cutStringData.c_str());

   RooRealVar B_plus_M("B_plus_M", "M(B^{+})",minBMass_data,maxBMass_data, "MeV");
   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M^{DTF}(B^{+})",B_plus_M_min_data,B_plus_M_max_data, "MeV");
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity", "e_plus_BremMultiplicity", -1, 3);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity", "e_minus_BremMultiplicity", -1, 3);
   RooRealVar J_psi_1S_M("J_psi_1S_M", "J_psi_1S_M", 0, 10000);
   RooRealVar passTrigCat( ("passTrigCat"+i2s(trigCat)).c_str(), ("passTrigCat"+i2s(trigCat)).c_str(), -10, 10);

   RooArgSet vars(B_plus_M, B_plus_DTFM_M_zero, e_plus_BremMultiplicity, e_minus_BremMultiplicity);
   vars.add(J_psi_1S_M);
   vars.add(passTrigCat);

   RooDataSet data(("data"+trigCatS).c_str(), ("data"+trigCatS).c_str(), vars, Import(*t), Cut( ("B_plus_M>"+d2s(minBMass_data)+" && B_plus_M <"+d2s(maxBMass_data)+" && B_plus_DTFM_M_zero>"+d2s(B_plus_DTF_M_cut)).c_str() ) );

   RooWorkspace workspaceDataData(("workspaceDataDataTrig"+i2s(trigCat)).c_str(), ("workspaceDataDataData"+i2s(trigCat)).c_str());

   workspaceDataData.import(data);

   TFile fws(workspaceFileName.c_str(), "UPDATE");
   workspaceDataData.Write("",TObject::kOverwrite );

   cout<<"Workspace containing data data ready: "<<endl;

   workspaceDataData.Print();
   cout<<"data has "<<data.sumEntries()<<" events"<<endl;

   f.Close();
   fws.Close();
   fTrash.Close();
}


void FitAndSplotKeeDataForTraining::prepareDataWorkspaceWithPoissonOsc(int trigCat, int nPoisson)
{
   string trigCatS("Trig"+i2s(trigCat));

   string totFileName(tupleDataDir+"/"+tupleDataName);

   TFile f( totFileName.c_str());
   TTree* t1 = (TTree*)f.Get(treeDataName.c_str());

   if(!t1)
   {
      cerr<<"ERROR: in function prepareDataWorkspaceWithPoissonOsc, no tree "<<treeDataName<<" found in "<<tupleDataDir<<"/"<<tupleDataName<<endl;
      return;
   }


   //copy tree with new vars

   vector<string> varsToSwitchOn;
   varsToSwitchOn.push_back("B_plus_DTFM_M_zero");
   varsToSwitchOn.push_back("e_plus_BremMultiplicity");
   varsToSwitchOn.push_back("e_minus_BremMultiplicity");
   varsToSwitchOn.push_back("J_psi_1S_M");
   varsToSwitchOn.push_back("poissonWeight");
   varsToSwitchOn.push_back(("passTrigCat"+i2s(trigCat)).c_str());

   if(cutStringData == "") cutStringData = "(1)";

   system("rm /vols/lhcb/palvare1/toErase/toErase.root");

   cout<<"Copying tree with new vars for data. Cut applied: "<<endl<<cutStringData<<endl<<endl;

   copyTreeWithNewVars("/vols/lhcb/palvare1/toErase/toErase.root", totFileName, cutStringData, "Sum$(poissonWeight["+i2s(nPoisson)+"])", "weightForFit", varsToSwitchOn, "DecayTree");

   //recover tree and prepare workspace


   TFile fTrash("/vols/lhcb/palvare1/toErase/toErase.root");
   TTree* t = (TTree*)fTrash.Get("DecayTree");

   if(!t)
   {
      cout<<"Tree: "<<t<<endl;
      cerr<<"ERROR: in function prepareDataWorkspaceWithPoissonOsc, copied tree not found."<<endl;
      return;
   }


   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M^{DTF}(B^{+})",B_plus_M_min_data,B_plus_M_max_data, "MeV");
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity", "e_plus_BremMultiplicity", -1, 3);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity", "e_minus_BremMultiplicity", -1, 3);
   RooRealVar J_psi_1S_M("J_psi_1S_M", "J_psi_1S_M", 0, 10000);
   RooRealVar passTrigCat( ("passTrigCat"+i2s(trigCat)).c_str(), ("passTrigCat"+i2s(trigCat)).c_str(), -10, 10);
   RooRealVar weight("weightForFit", "weightForFit", 0, 10);


   B_plus_DTFM_M_zero.setBins(50);

   RooArgSet vars(B_plus_DTFM_M_zero, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weight, J_psi_1S_M, passTrigCat);


   RooDataSet data(("data"+trigCatS).c_str(), ("data"+trigCatS).c_str(), vars, Import(*t), Cut( ("B_plus_DTFM_M_zero>"+d2s(B_plus_M_min_data)+" && B_plus_DTFM_M_zero <"+d2s(B_plus_M_max_data)).c_str() ),
         WeightVar( weight.GetName()  ) );

   RooWorkspace workspaceDataData(("workspaceDataDataTrig"+i2s(trigCat)).c_str(), ("workspaceDataDataData"+i2s(trigCat)).c_str());

   workspaceDataData.import(data);

   TFile fws(workspaceFileName.c_str(), "UPDATE");
   workspaceDataData.Write("",TObject::kOverwrite );

   cout<<"Workspace containing poisson oscillated data data ready: "<<endl;

   workspaceDataData.Print();
   cout<<"data has "<<data.sumEntries()<<" events"<<endl;

   f.Close();
   fws.Close();
   fTrash.Close();
}


void FitAndSplotKeeDataForTraining::plotOnMC(RooAbsPdf *model, RooRealVar* var, int trigCat)
{



   vector<string> varsToSwitchOn;
   varsToSwitchOn.push_back("B_plus_DTFM_M_zero");
   varsToSwitchOn.push_back("B_plus_M");
   varsToSwitchOn.push_back("e_plus_BremMultiplicity");
   varsToSwitchOn.push_back("e_minus_BremMultiplicity");
   varsToSwitchOn.push_back("J_psi_1S_M");
   varsToSwitchOn.push_back("FidelCastro_ETOS");
   varsToSwitchOn.push_back("FidelCastro_HTOS");
   varsToSwitchOn.push_back("FidelCastro_TIS");
   varsToSwitchOn.push_back("FidelCastro_KtoPi_trigCat0");
   varsToSwitchOn.push_back("FidelCastro_KtoPi_trigCat1");
   varsToSwitchOn.push_back("FidelCastro_KtoPi_trigCat2");
   varsToSwitchOn.push_back("BDT1ETOSRun1");


   if(cutStringMC == "") cutStringMC = "(1)";
   if(weightMC == "") weightMC = "(1)";

   system("rm /vols/lhcb/palvare1/toErase/toErase.root");

   cout<<"Copying tree with new vars. Cut applied: "<<endl<<cutStringMC<<endl<<endl;

   if (turnOffWeightsSignal) copyTreeWithNewVars("/vols/lhcb/palvare1/toErase/toErase.root", (tupleMCdir+"/"+tupleMCname).c_str(), 
                                                 cutStringMC, "(1)", "weightForFit", varsToSwitchOn, "DecayTree");
   else copyTreeWithNewVars("/vols/lhcb/palvare1/toErase/toErase.root", (tupleMCdir+"/"+tupleMCname).c_str(), 
                                                 cutStringMC, weightMC, "weightForFit", varsToSwitchOn, "DecayTree");
   //recover tree and prepare workspace
   TFile fTrash("/vols/lhcb/palvare1/toErase/toErase.root");
   TTree* t = (TTree*)fTrash.Get("DecayTree");

   if(!t)
   {
      cout<<"Tree: "<<t<<endl;
      cerr<<"ERROR: in function prepareWorkspace, no tree found in "<<tupleMCdir<<"/"<<tupleMCname<<endl;
      return;
   }


   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M^{DTF}(B^{+})",B_plus_M_min_MC,B_plus_M_max_MC, "MeV");
   RooRealVar weight("weightForFit", "weightForFit", 0, 10);


   B_plus_DTFM_M_zero.setBins(50);

   RooArgSet vars(B_plus_DTFM_M_zero, *var, weight);


   RooDataSet data("data", "data", vars, 
                   Import(*t), 
                   Cut( ("B_plus_DTFM_M_zero>"+d2s(B_plus_DTF_M_cut)+" && B_plus_M >"+d2s(minBMass_MC)+" && B_plus_M <"
                         +d2s(maxBMass_MC)).c_str() ),
                   WeightVar( weight.GetName()) );   

  RooPlot *frame = var->frame();
  data.plotOn(frame, Binning(50));
  model->plotOn(frame);

  savePullPlot(*frame, plotdir+"pullPlot.root");
  TFile fpull((plotdir+"pullPlot.root").c_str());
  TCanvas* cpull = (TCanvas*)fpull.Get("pullplot");

  TCanvas canv("canv", "canv", 600, 600);
  frame->SetMinimum(0.5);
  frame->Draw();

  TCanvas canvTot("canvTot", "canvTot", 600, 600);

  canvTot.Divide(1,2);
  canvTot.cd(1);
  gPad->SetPad(0.005, 0.205, 0.995, 0.995);
  canv.DrawClonePad();
  canvTot.cd(2);
  gPad->SetPad(0.005, 0.005, 0.995, 0.2);
  cpull->DrawClonePad();

  canvTot.Print((plotdir+"totalMCPDF"+i2s(trigCat)+".pdf").c_str());

  canv.SetLogy();
  canvTot.cd(1);

  canv.DrawClonePad();
  canvTot.Print((plotdir+"totalMCPDFLogy"+i2s(trigCat)+".pdf").c_str());

  fpull.Close();
  fTrash.Close();


}



void FitAndSplotKeeDataForTraining::getFracPiee(int trigCat, RooRealVar& frac)
{
   // string cut("B_plus_DTFM_M_zero > "+d2s(B_plus_M_min_data)+" && B_plus_DTFM_M_zero < "+d2s(B_plus_M_max_data));
  string cut("B_plus_DTFM_M_zero > "+d2s(B_plus_DTF_M_cut)+" && B_plus_M > "+d2s(minBMass_data)+" && B_plus_M < "+d2s(maxBMass_MC));
  
  string weightBefore;
  string weightBeforePiee;
  string weightAfter;
  string weightAfterPiee;
  
  // using the B->JpsiK K->pi
  weightBefore = "e_plus_PIDWeight_eff*e_minus_PIDWeight_eff*DataCondKinWeight_MUTOS*trigWeightR_trigCat0";
  weightBeforePiee = "e_plus_PIDWeight_eff*e_minus_PIDWeight_eff*DataCondKinWeight_MUTOS*trigWeightR_trigCat0";
  weightAfter = weightMC;
  weightAfterPiee = weightMCPiee;

  // if (turnOffWeightsSignal) {
  //   // weightBefore = "( e_plus_PIDe > 3 && e_minus_PIDe > 3)";
  //   // weightAfter = "(K_Kst_ProbNNk>0.2 && e_plus_PIDe > 3 && e_minus_PIDe > 3)";
  //   weightBefore = "(1)";
  //   weightAfter = "(1)";
  // }

  // weightBefore = "dataCondWeight";
  // // weightAfter = weightMC;
  // weightAfter = "passTrigCat"+i2s(trigCat);

  double nPiBefore( getEntries(tuplePieeDir+tuplePieeName, cut, weightBeforePiee, treePieeName) );
  double nKBefore( getEntries(tupleMCdir+tupleMCname, cut, weightBefore, treeMCname) );
  
  if(cutStringMC != "") cut += (" && "+cutStringMC);

  double nPiAfter( getEntries(tuplePieeDir+tuplePieeName, cut, weightAfterPiee, treePieeName) );
  double nKAfter( getEntries(tupleMCdir+tupleMCname, cut, weightAfter, treeMCname) );

  double effPi(nPiAfter/nPiBefore);
  double effK (nKAfter/nKBefore);
  
  double sigmaK(1.027e-3);
  double sigmaPi(4.1e-5);

   frac.SetName( ("fracPieeTrig"+i2s(trigCat)).c_str() );
   frac.SetTitle( ("fracPieeTrig"+i2s(trigCat)).c_str() );
   frac.setVal( sigmaPi*effPi/ (sigmaK*effK) ); 
   frac.setError( 0.1*frac.getVal() );

   cout<<"PIEE PIEE PIEE FRAC: in trigcat "<<trigCat<<", "<<frac.getVal()<<" +- "<<frac.getError()<<endl;
}





//void FitAndSplotKeeDataForTraining::getFracRare(int trigCat, RooRealVar& frac)
//{
//   string cut("B_plus_DTFM_M_zero > "+d2s(B_plus_M_min_data)+" && B_plus_DTFM_M_zero < "+d2s(B_plus_M_max_data)); 
//
//   if(cutString != "") cut += (" && "+cutString);
//
//   double sigmaK(1.027e-3*0.05971);
//   double sigmaRare(4.5e-7);
//
//   double effQ2K(0.88);
//   double effQ2Rare(0.24);
//
//   string weightBefore;
//   string weightAfter;
//   
//
//   weightBefore = "dataCondWeight";
//   weightAfter = "PIDTrigDataCondWeight_TrigCat"+i2s(trigCat);
//
//   double nRareBefore( getEntries(tupleRareDir+tupleRareName, "", weightBefore, treeRareName) );
//   double nRareAfter( getEntries(tupleRareDir+tupleRareName, cut, weightAfter, treeRareName) );
//
//   double nKBefore( getEntries(tupleMCdir+tupleMCname, "", weightBefore, treeMCname) );
//   double nKAfter( getEntries(tupleMCdir+tupleMCname, cut, weightAfter, treeMCname) );
//   
//   double effRare(effQ2Rare*nRareAfter/nRareBefore);
//   double effK(effQ2K*nKAfter/nKBefore);
//
//   frac.SetName( ("fracRareTrig"+i2s(trigCat)).c_str() );
//   frac.SetTitle( ("fracRareTrig"+i2s(trigCat)).c_str() );
//   frac.setVal( sigmaRare*effRare/ (sigmaK*effK) );
//
//
//   frac.setError( 1.*frac.getVal() );
//
//   cout<<"RARE RARE RARE FRAC: in trigcat "<<trigCat<<", "<<frac.getVal()<<" +- "<<frac.getError()<<endl;
//
//}



void FitAndSplotKeeDataForTraining::performFullFitUseJPsiModeConstraints(int trigCat, RooRealVar& yield, bool fast, bool saveYield, string workspaceFileName)
{
   string trigCatS("Trig"+i2s(trigCat));

   TFile f(workspaceFileName.c_str());
   RooWorkspace* ws = (RooWorkspace*)f.Get(("workspaceDataFit"+trigCatS).c_str());

   if(!ws)
   {
      cerr<<"ERROR: in FitAndSplotKSeeDataForTraining::performFullFitUseJPsiModeConstraints, no workspace found."<<endl;
      f.Close();
      return;
   }

   RooRealVar* scaleRRV = ws->var(("sigmaScaleFactor"+trigCatS).c_str());
   RooRealVar* shiftRRV = ws->var(("meanShift"+trigCatS).c_str());

   if(!scaleRRV || !shiftRRV)
   {
      cerr<<"ERROR: in FitAndSplotKeeDataForTraining::performFullFitUseJPsiModeConstraints, problem downloading shift and scale. "<<scaleRRV<<" "<<shiftRRV<<endl;
      f.Close();
      return;
   }

   double scaleFix(scaleRRV->getVal());
   double shiftFix(shiftRRV->getVal());

   performFullFit(trigCat, yield, fast, saveYield, shiftFix, scaleFix);
}



void FitAndSplotKeeDataForTraining::performFullFit(int trigCat, RooRealVar& yield, bool fast, bool saveYield,  double shiftFix, double scaleFix)
{
   system( ("mkdir -p "+plotdir).c_str());

   prepareDataWorkspace(trigCat );

   prepareMCWorkspace(trigCat, -1);
   prepareMCWorkspace(trigCat, 0);
   prepareMCWorkspace(trigCat, 1);
   prepareMCWorkspace(trigCat, 2);

   prepareRarePrcWorkspace(trigCat, -1);
   prepareCharmPrcWorkspace(trigCat, -1);
   //prepareRareWorkspace(trigCat, -1);
   preparePieeWorkspace(trigCat,-1);

   fitPrc(trigCat, -1, "Rare"); //whichOne = "Rare" Or whichOne = "Charm"
   fitPrc(trigCat, -1, "Charm"); //whichOne = "Rare" Or whichOne = "Charm"
   fitPiee(trigCat, -1, true);

   fitMCAuto(trigCat, 0, fast);
   fitMCAuto(trigCat, 1, fast);
   fitMCAuto(trigCat, 2, fast);

   fitData(trigCat, yield, fast, saveYield, shiftFix, scaleFix);
}



void FitAndSplotKeeDataForTraining::fillHistoWithYields(string nameHistoFile, string _plotdir, int trigCat, int nBins, double binTab[], string varToCutOn, double shiftFix, double scaleFix)  
{
   TFile f(nameHistoFile.c_str(), "UPDATE" );
   TH1F h(("yieldsTrig"+i2s(trigCat)).c_str(), ("yieldsTrig"+i2s(trigCat)).c_str(), nBins, binTab); 
   h.Sumw2();

   string plotDir2;
   RooRealVar sig("sig", "sig", 0);

   system( ("mkdir -p "+_plotdir).c_str() );

   string saveCutMC(cutStringMC);
   string saveCutData(cutStringData);
   string saveDir(plotdir);

   for(int i(1); i<=nBins; ++i)
   {
      cutStringMC = varToCutOn +">"+d2s(binTab[i-1])+" && "+varToCutOn+"<"+d2s(binTab[i]);
      if(saveCutMC != "") cutStringMC += " && " + saveCutMC;

      cutStringData = varToCutOn +">"+d2s(binTab[i-1])+" && "+varToCutOn+"<"+d2s(binTab[i]);
      if(saveCutData != "") cutStringData += " && " + saveCutData;

      plotdir = _plotdir+"/bin"+i2s(i)+"/";

      performFullFit(trigCat, sig, true, false, shiftFix, scaleFix); 

      h.SetBinContent(i, sig.getVal());
      h.SetBinError(i, sig.getError());

      cout<<endl<<endl<<"fit "<<i<<"/"<<nBins<<" done"<<endl<<endl<<endl;
   }

   h.GetXaxis()->SetTitle(varToCutOn.c_str() );
   h.GetYaxis()->SetTitle("signal yield" );

   f.cd();
   h.Write("", TObject::kOverwrite );

   TCanvas canv("canv", "canv", 600, 600);
   h.Draw("E1");
   canv.Print((_plotdir+"/yieldsTrig"+i2s(trigCat)+".pdf").c_str());

   f.Close();

   plotdir = saveDir;
   cutStringMC = saveCutMC;
   cutStringData = saveCutData;
}



void FitAndSplotKeeDataForTraining::updatePoissonHisto(string nameHistoFile, string _plotdir, int trigCat, int binToFill, int nPoisson, string varToCutOn, bool refitMC, 
   double shiftFix, double scaleFix)  
{
   TFile f(nameHistoFile.c_str(), "UPDATE" );
   TH1F* h = (TH1F*)f.Get(("yieldsTrig"+i2s(trigCat)).c_str()); 

   if(!h)
   {
      cerr<<"ERROR: in FitAndSplotKeeDataForTraining::updatePoissonHisto, no histogram found!"<<endl;
      f.Close();
      return;
   }

   int nBins(h->GetNbinsX());

   double binTab[nBins+1];

   for(int i(0); i <= nBins; ++i) binTab[i] = h->GetXaxis()->GetBinUpEdge(i);

   if(binToFill < 1 || binToFill > nBins)
   {
      cerr<<"ERROR: in FitAndSplotKeeDataForTraining::updatePoissonHisto, invalid binToFill: "<<binToFill<<endl;
      f.Close();
      return;
   }

   string plotDir2;
   RooRealVar sig("sig", "sig", 0);



   string saveCutMC(cutStringMC);
   string saveCutData(cutStringData);
   string saveDir(plotdir);

   cutStringMC = varToCutOn +">"+d2s(binTab[binToFill-1])+" && "+varToCutOn+"<"+d2s(binTab[binToFill]);
   if(saveCutMC != "") cutStringMC += " && " + saveCutMC;

   cutStringData = varToCutOn +">"+d2s(binTab[binToFill-1])+" && "+varToCutOn+"<"+d2s(binTab[binToFill]);
   if(saveCutData != "") cutStringData += " && " + saveCutData;

   plotdir = _plotdir+"/bin"+i2s(binToFill)+"/";

   system( ("mkdir -p "+plotdir).c_str() );

   if(refitMC)
   {
      prepareMCWorkspace(trigCat, -1);
      prepareMCWorkspace(trigCat, 0);
      prepareMCWorkspace(trigCat, 1);
      prepareMCWorkspace(trigCat, 2);

      prepareRarePrcWorkspace(trigCat, -1);
      prepareCharmPrcWorkspace(trigCat, -1);
      preparePieeWorkspace(trigCat,-1);

      fitPrc(trigCat, -1, "Rare"); //whichOne = "Rare" Or whichOne = "Charm"
      fitPrc(trigCat, -1, "Charm"); //whichOne = "Rare" Or whichOne = "Charm"
      fitPiee(trigCat, -1, true);

      fitMCAuto(trigCat, 0, true);
      fitMCAuto(trigCat, 1, true);
      fitMCAuto(trigCat, 2, true);
   }

   prepareDataWorkspaceWithPoissonOsc(trigCat, nPoisson);
   fitData(trigCat, sig, true, false, shiftFix, scaleFix);

   h->SetBinContent(binToFill, sig.getVal());
   h->SetBinError(binToFill, sig.getError());

   h->GetXaxis()->SetTitle(varToCutOn.c_str() );
   h->GetYaxis()->SetTitle("signal yield" );

   f.cd();
   h->Write("", TObject::kOverwrite );

   TCanvas canv("canv", "canv", 600, 600);
   h->Draw("E1");
   canv.Print((_plotdir+"/yieldsTrig"+i2s(trigCat)+".pdf").c_str());

   f.Close();

   plotdir = saveDir;
   cutStringMC = saveCutMC;
   cutStringData = saveCutData;
}




void fitKeeDataPerformFitWithCut(string plotDir, int trigCat, RooRealVar& yield, string MCCut, string dataCut, string MCWeight)
{
   FitAndSplotKeeDataForTraining eefit;

   eefit.cutStringData = dataCut;
   eefit.cutStringMC = MCCut;
   eefit.plotdir = plotDir;
   eefit.weightMC = MCWeight;
   eefit.workspaceFileName = ("/vols/lhcb/palvare1/RK_analysis/FitWorkspaces/FitAndSplotKeeDataForTrainingWorkspaceCUT.root");
   system( ("mkdir -p "+plotDir).c_str());

   eefit.prepareDataWorkspace(trigCat );

   eefit.prepareRarePrcWorkspace(trigCat, -1);
   eefit.prepareCharmPrcWorkspace(trigCat, -1);
   //eefit.prepareRareWorkspace(trigCat, -1);
   eefit.preparePieeWorkspace(trigCat,-1);

   eefit.prepareMCWorkspace(trigCat, 0);
   eefit.prepareMCWorkspace(trigCat, 1);
   eefit.prepareMCWorkspace(trigCat, 2);

   eefit.fitPrc(trigCat, -1, "Rare"); //whichOne = "Rare" Or whichOne = "Charm"
   eefit.fitPrc(trigCat, -1, "Charm"); //whichOne = "Rare" Or whichOne = "Charm"
   eefit.fitPiee(trigCat, -1, true);
   //eefit.fitRare(trigCat, -1);

   eefit.fitMCAuto(trigCat, 0, true);
   eefit.fitMCAuto(trigCat, 1, true);
   eefit.fitMCAuto(trigCat, 2, true);

   eefit.fitData(trigCat, yield, true, false);
}


void fitKeeDataFillHistoWithYields(string nameHistoFile, string plotDir, int trigCat, int nBins, double binTab[], string varToCutOn, string MCCut, string dataCut, string MCWeight)  
{
   TFile f(nameHistoFile.c_str(), "UPDATE" );
   TH1F h(("yieldsTrig"+i2s(trigCat)).c_str(), ("yieldsTrig"+i2s(trigCat)).c_str(), nBins, binTab); 
   h.Sumw2();

   string totCutMC;
   string totCutData;
   string plotDir2;
   RooRealVar sig("sig", "sig", 0);

   system( ("mkdir -p "+plotDir).c_str() );

   for(int i(1); i<=nBins; ++i)
   {
      totCutMC = MCCut + " && " + varToCutOn +">"+d2s(binTab[i-1])+" && "+varToCutOn+"<"+d2s(binTab[i]);
      totCutData = dataCut + " && " + varToCutOn +">"+d2s(binTab[i-1])+" && "+varToCutOn+"<"+d2s(binTab[i]);
      plotDir2 = plotDir+"/bin"+i2s(i)+"/";

      fitKeeDataPerformFitWithCut(plotDir2, trigCat, sig, totCutMC, totCutData, MCWeight); 

      h.SetBinContent(i, sig.getVal());
      h.SetBinError(i, sig.getError());

      cout<<endl<<endl<<"fit "<<i<<"/"<<nBins<<" done"<<endl<<endl<<endl;
   }

   h.GetXaxis()->SetTitle(varToCutOn.c_str() );
   h.GetYaxis()->SetTitle("signal yield" );

   f.cd();
   h.Write("", TObject::kOverwrite );

   TCanvas canv("canv", "canv", 600, 600);
   h.Draw("E1");
   canv.Print((plotDir+"/yieldsTrig"+i2s(trigCat)+".pdf").c_str());

   f.Close();
}


void fitKeeDataFillHistoWithYieldsUniform(string nameHistoFile, string plotDir, int trigCat, int nBins, double binTab[], string varToCutOn, string MCCut, string dataCut,  string MCWeight)  
{
   TFile f(nameHistoFile.c_str(), "UPDATE" );
   TH1F h(("yieldsTrig"+i2s(trigCat)).c_str(), ("yieldsTrig"+i2s(trigCat)).c_str(), nBins, 0., 1.*nBins); 
   h.Sumw2();

   string totCutMC;
   string totCutData;

   string plotDir2;
   RooRealVar sig("sig", "sig", 0);

   system( ("mkdir "+plotDir).c_str() );

   for(int i(1); i<=nBins; ++i)
   {
      totCutMC = MCCut + " && " + varToCutOn +"=="+d2s(binTab[i-1]);
      totCutData = dataCut + " && " + varToCutOn +"=="+d2s(binTab[i-1]);
      plotDir2 = plotDir+"/bin"+i2s(i)+"/";

      fitKeeDataPerformFitWithCut(plotDir2, trigCat, sig, totCutMC, totCutData, MCWeight); 

      h.SetBinContent(i, sig.getVal());
      h.SetBinError(i, sig.getError());
      
      cout<<endl<<endl<<"fit "<<i<<"/"<<nBins<<" done"<<endl<<endl<<endl;
   }

   h.GetXaxis()->SetTitle(varToCutOn.c_str() );
   h.GetYaxis()->SetTitle("signal yield" );

   f.cd();
   h.Write("", TObject::kOverwrite );

   TCanvas canv("canv", "canv", 600, 600);
   h.Draw("E1");
   canv.Print((plotDir+"/yieldsTrig"+i2s(trigCat)+".pdf").c_str());

   f.Close();
}

//  LocalWords:  cutStringMC cutStringMCSave
