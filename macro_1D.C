#include "RooExpBinned.h"
#include "RooExponential.h"

void macro_1D()
{


  RooRealVar B_plus_M("B_plus_M", "M_{visible}", 4800, 6700, "MeV");
  RooRealVar misPT("misPT", "p_{#perp}", 0, 5000, "MeV");

  std::string combfile ="/vols/lhcb/palvare1/B2Kee_Kmumu/Marcin/B2XMuMu_Line_Strip21_Kemu_presel.root";
  std::string combtree ="DecayTree";
  std::string combcuts = "K_Kst_isMuon == 0";
  std::vector<std::string> varset = {"B_plus_M","K_Kst_isMuon"};
  

  TFile* fComb = new TFile(combfile.c_str());
  TTree* tComb = (TTree*)fComb->Get(combtree.c_str());
  
  RooRealVar *var;
  RooArgSet argset;
     
  tComb->SetBranchStatus("*",0);
  for(vector<string>::const_iterator i = varset.begin(); i != varset.end(); ++i) 
  {
    tComb->SetBranchStatus((*i).c_str(), 1);
    if (*i=="B_plus_M") continue;
    var = new RooRealVar((*i).c_str(),(*i).c_str(),-10,10);
    argset.add(*var);
  }
  
  argset.add(B_plus_M);

  RooBinning mass_binning(floor((-4800+6200)/40.),4800,6200,"mass_binning");
  
  B_plus_M.setBinning(mass_binning );

  RooRealVar l1Kee("l1Kee", "l1Kee", +2.3e-3, -1e-1, 5e-1);
  RooExpBinned* combPDF =  new RooExpBinned("combPDF", "combPDF", B_plus_M, 
                                                                    l1Kee, 
                                                                    mass_binning);

  // RooExponential* combPDF =  new RooExponential("combPDF", "combPDF", B_plus_M,  
  //                                                       l1Kee);
  
  RooDataSet *dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, combcuts.c_str());
  RooDataHist *dataHistComb = dataSetComb->binnedClone("dataHistComb","dataHistComb");
  
  
  combPDF->fitTo(*dataHistComb);
  

  TH1D *model_bin_x = (TH1D*) combPDF->createIntegral(RooArgSet(misPT))->createHistogram("model_bin",B_plus_M);
  
  
  // model_bin_x->Scale(dataSetComb->sumEntries()*1.*mass_binning.averageBinWidth()/model_bin_x->GetSumOfWeights());
  // model_bin_y->Scale(dataSetComb->sumEntries()*1.*misPT_binning.averageBinWidth()/model_bin_y->GetSumOfWeights());

  model_bin_x->Scale(dataSetComb->sumEntries()*1./model_bin_x->GetSumOfWeights());
  
  
  model_bin_x->SetLineColor(kRed);
  
  
  // B_plus_M.setVal(5500);
  // misPT.setVal(2000);
  // double val = combPDF->getVal();
  // cout<<"combPDF = "<<val<<endl;


  TCanvas *cv1 = new TCanvas("cv1","cv1");
  RooPlot *frame_m = B_plus_M.frame();
  dataHistComb->plotOn(frame_m);
  combPDF->plotOn(frame_m);
  frame_m->Draw();
  model_bin_x->Draw("same");
  


}


