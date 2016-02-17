#include<iostream>
#include<boost/filesystem/operations.hpp>
#include<boost/filesystem/path.hpp>
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
#include"RooMcorMvisTsallis.h"
#include"usefulFunctions.h"
#include"fitter_utils.h"

namespace fs = boost::filesystem;


using namespace std;
using namespace RooFit;

int main(int argc, char* argv[])
{
   //***********Number of events to generata



   if(argc != 9)
   {
      cout<<"maintables take the following arguments:"<<endl;  
      cout<<"name of the file"<<endl;
      cout<<"name of the tree in the file"<<endl;
      cout<<"nGen Signal"<<endl;
      cout<<"nGen PartReco"<<endl;
      cout<<"nGen combinatorial"<<endl;
      cout<<"nGen J/psi leak"<<endl;
      cout<<"name output"<<endl;
      cout<<"want output in latex? (1=yes, 0=no)"<<endl;
   }

   string filename;
   string treename;
   int nGenSignal(203);
   int nGenPartReco(79);
   int nGenComb(31);
   int nGenJpsiLeak(10);
   string outputname;
   bool wantLatex(false);

   if(argc == 9)
   {
      filename = argv[1];
      treename = argv[2];
      nGenSignal = atoi(argv[3]);
      nGenPartReco = atoi(argv[4]);
      nGenComb = atoi(argv[5]);
      nGenJpsiLeak = atoi(argv[6]);
      outputname = argv[7];
      if(*argv[8] == '1') wantLatex = true;
   }

   ofstream out(outputname.c_str());

   cout<<"HELLO1"<<endl;
   makeTableResults(filename, treename, nGenSignal, nGenPartReco, nGenComb, nGenJpsiLeak, out, wantLatex);
   cout<<"HELLO2 "<<outputname<<endl;

   out.close();

   return 0;
}


