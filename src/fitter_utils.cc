#include "fitter_utils.h"


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



void initiateParams(RooArgSet& parset)
{

  RooRealVar *var;
  TIterator *iter = parset.createIterator();

  while((var = (RooRealVar*) iter->Next()))
  {
    var->randomize();
  }
}


void prepare_PDFs(string workspacefile, string trigStr, string BDTcut, bool fit2D,
                  string signalfile, string partrecofile, string combfile,
                  string signaltree, string partrecotree, string combtree)
{


  //***********Get the datasets
  TFile* fSignal = new TFile(signalfile.c_str());
  TTree* tSignal = (TTree*)fSignal->Get(signaltree.c_str());
  TFile* fPartReco = new TFile(partrecofile.c_str());
  TTree* tPartReco = (TTree*)fPartReco->Get(partrecotree.c_str());
  TFile* fComb = new TFile(combfile.c_str());
  TTree* tComb = (TTree*)fComb->Get(combtree.c_str()); 


  //**********Define variables
  RooRealVar trigVar(trigStr.c_str(), trigStr.c_str(), -10, 10);
  RooRealVar BDTKeeBig("BDTKeeBig3", "BDTKeeBig3", -1,1);
  RooRealVar B_plus_M("B_plus_M", "M_{visible}", 4880, 5700, "MeV/c^{2}");
  RooRealVar B_plus_M_corr("B_plus_M_corr", "M_{cor}", 5000, 7000, "MeV/c^{2}");
  RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV/c^{2}"); 
  RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
  RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);


  //***********Set Binning
  B_plus_M.setBins(20);
  B_plus_M_corr.setBins(15);
  B_plus_DTFM_M_zero.setBins(100);

  RooArgSet argset(BDTKeeBig, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity);

  cout<<"getting the datasets:"<<endl;
  
  RooDataSet* dataSetSignalZeroGamma;
  RooDataSet* dataSetSignalOneGamma;
  RooDataSet* dataSetSignalTwoGamma;
  RooDataSet* dataSetPartReco;
  RooDataSet* dataSetComb;

  TFile* fw;

  dataSetSignalZeroGamma = new RooDataSet("dataSetSignalZeroGamma", "dataSetSignalZeroGamma", tSignal, argset,( " ("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5)").c_str());
  dataSetSignalOneGamma = new RooDataSet("dataSetSignalOneGamma", "dataSetSignalOneGamma", tSignal, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5)").c_str());
  dataSetSignalTwoGamma = new RooDataSet("dataSetSignalTwoGamma", "dataSetSignalTwoGamma", tSignal, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+") && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5)").c_str());
  dataSetPartReco = new RooDataSet("dataSetPartReco", "dataSetPartReco", tPartReco, argset, ("("+trigStr+"  > 0.9) && (BDTKeeBig3> "+BDTcut+")").c_str());
  dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, argset, ("("+trigStr+"  > 0.9)").c_str());


   cout<<"binning the datasets:"<<endl;

   RooArgSet argset2(B_plus_M);
   if (fit2D) argset2.add(B_plus_M_corr);

   RooDataHist dataHistSignalZeroGamma("dataHistSignalZeroGamma", "dataHistSignalZeroGamma", argset2, *dataSetSignalZeroGamma); 
   RooDataHist dataHistSignalOneGamma("dataHistSignalOneGamma", "dataHistSignalOneGamma", argset2, *dataSetSignalOneGamma); 
   RooDataHist dataHistSignalTwoGamma("dataHistSignalTwoGamma", "dataHistSignalTwoGamma", argset2, *dataSetSignalTwoGamma); 
   RooDataHist dataHistComb("dataHistComb", "dataHistComb", argset2, *dataSetComb); 
   RooDataHist dataHistPartReco("dataHistPartReco", "dataHistPartReco", argset2, *dataSetPartReco); 

   //***************Create 2D histogram estimates from data


   cout<<"Preparing the 3 2D histPdf: 1"<<endl;
   //   RooArgSet argset2(B_plus_M);
   RooHistPdf histPdfSignalZeroGamma("histPdfSignalZeroGamma", "histPdfSignalZeroGamma", argset2, dataHistSignalZeroGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfSignalOneGamma("histPdfSignalOneGamma", "histPdfSignalOneGamma", argset2, dataHistSignalOneGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfSignalTwoGamma("histPdfSignalTwoGamma", "histPdfSignalTwoGamma", argset2, dataHistSignalTwoGamma,2); cout<<" 2"<<endl;
   RooHistPdf histPdfPartReco("histPdfPartReco", "histPdfPartReco", argset2, dataHistPartReco,2); cout<<" 3"<<endl;


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
   workspace.import(*dataSetSignalZeroGamma);
   workspace.import(*dataSetSignalOneGamma);
   workspace.import(*dataSetSignalTwoGamma);
   workspace.import(*dataSetPartReco);
   workspace.import(*dataSetComb);
   workspace.import(histPdfSignalZeroGamma);
   workspace.import(histPdfSignalOneGamma);
   workspace.import(histPdfSignalTwoGamma);
   workspace.import(histPdfPartReco);
   workspace.import(*combPDF);
   workspace.writeToFile(workspacefile.c_str());
   

}
