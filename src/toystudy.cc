#include<iostream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
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
#include "fitter_utils.h"

namespace fs = boost::filesystem;


using namespace std;
using namespace RooFit;

int main(int argc, char* argv[])
{

   //***********Number of events to generate
  
  
  int nGenSignal(203);
  int nGenPartReco(79);
  int nGenComb(31);

  double nGenFracZeroGamma(0.368);
  double nGenFracOneGamma(0.484);

  
  //*********** Get arguments and set stuff
  string trigStr("L0ETOSOnly_d");
  bool wantOldDataSet(false);
  string BDTcut("0.257");
  int ntoys(1000);
  bool fit2D(0);
  bool constPartReco(0);
  string outputfolder="fit_result";
  

  if(argc != 11)
  {
    cout<<"toystudy:  Launchs a ToyMC study of the fitter specified by the options"<<endl;
    cout<<"Syntax: "<<argv[0]<<" <no_reload> <nsignal> <npartreco> <ncomb> <trigger_cat> <ntoys> <bdt_cut> <constPartReco> <ndims> <output_folder>"<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"no_reload: 0 (build PDFs from control channel), 1 (use PDFs in workspace)"<<endl;
    cout<<"trig_cat: 0 (ETOS), 1 (HTOS), 2 (TIS)"<<endl;
    cout<<"constPartReco: Constrain part. reco. fraction? 0 (No), 1 (Yes)"<<endl;
    cout<<"ndims: 1 (Mvis), 2 (Mvis x Mcorr)"<<endl;
  
    return 0;
  }
  

  if(argc == 11)
  {
    if(*argv[1] == '1') wantOldDataSet = true;
    
    nGenSignal = atoi(argv[2]);
    nGenPartReco = atoi(argv[3]);
    nGenComb = atoi(argv[4]);
    
    int trigInt(atoi(argv[5]));
    if(trigInt == 1) 
    {
      trigStr = "L0HTOSOnly_d";
      nGenFracZeroGamma = 0.322;
      nGenFracOneGamma = 0.484;
    }      
    if(trigInt == 2) 
    {
      trigStr = "L0TISOnly_d";
      nGenFracZeroGamma = 0.330;
      nGenFracOneGamma = 0.495;
    }
    
    BDTcut = argv[7];
    ntoys = atoi(argv[6]);

    if(*argv[8] == '1') constPartReco = true;
    if(*argv[9] == '2') fit2D = true;
    outputfolder = argv[10];


    fs::path data_dir("./"+outputfolder);
    if (!fs::is_directory(data_dir)){
      cout<<"Output directory "<<outputfolder<<" does not exist."<<endl;
      return 0;
    }
    
  }
  
  //*************** Output files

  string plotsfile = outputfolder+"/plots2DHistBremCatTsallisBkgfit.root";
  string resultsfile = outputfolder+"/toystudy2DHistBremCatTsallisBkg_results.root";
  string outfile = outputfolder+"/fitResult2D.dat";
  string tablefile = outputfolder+"/resultToy.dat";
  

   //***********Get the datasets
   string fSignal = "/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldtrees/B2Kee_Strip21_BDT_ctrl_trigged.root";
   string fPartReco = "/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldtrees/B2Kee_Strip21_BDT_prc_trigged.root";
   string fComb = "/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldtrees/B2Kee_Strip21_piee_trigged.root";
   string workspacename = outputfolder+"/mypdfs.root";
   
   if (!wantOldDataSet) prepare_PDFs(workspacename, trigStr, BDTcut, fit2D, fSignal, fPartReco, fComb);


   //***************Prepare the stuff to generate events


   TFile f(resultsfile.c_str(),"recreate");
   TTree t("paramsFloatingExpFloatingFracPartReco_L0ETOS", "paramsFloatingExpFloatingFracPartReco_L0ETOS");

   bool wantPlots(true);
   bool update(false);
   
   ofstream out(outfile);

   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Generation and fit number "<<i<<endl;

      generate_and_fit( workspacename, fit2D, wantPlots, constPartReco, 
                nGenSignal,   nGenPartReco,   nGenComb,
                nGenFracZeroGamma,   nGenFracOneGamma,
                out, &t,  update, plotsfile);


      wantPlots = false;
      update = true;
   }

   out.close();

   f.cd();
   t.Write();

   ofstream out2(tablefile, std::ios::app);
   if (fit2D) out2<<"2D Mvis X MCorr with Tsallis bkg:"<<endl;
   else out2<<"1D Mvis:"<<endl;
   makeTableResults(&t, nGenSignal, nGenPartReco, nGenComb, out2 );
   out2.close();


   f.Close();

   return 0;
}


