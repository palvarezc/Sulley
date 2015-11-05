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
  
  
  int nGenSignal(203);
  int nGenPartReco(79);
  int nGenComb(31);
  int nGenJpsiLeak(10);

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
  bool wantHOPCut(false);
  double minBMass(4880);
  double maxBMass(5700);

  if(argc != 15)
  {
    cout<<"toystudy:  Launchs a ToyMC study of the fitter specified by the options"<<endl;
    cout<<"Syntax: "<<argv[0]<<" <no_reload> <nsignal> <npartreco> <ncomb> <nJpsiLeak> <trigger_cat> <ntoys> <bdt_cut> <constPartReco> <ndims> <output_folder> <HOP cut> <minBMass> <maxBMass>"<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"no_reload: 0 (build PDFs from control channel), 1 (use PDFs in workspace)"<<endl;
    cout<<"trig_cat: 0 (ETOS), 1 (HTOS), 2 (TIS)"<<endl;
    cout<<"constPartReco: Constrain part. reco. fraction? 0 (No), 1 (Yes)"<<endl;
    cout<<"ndims: 1 (Mvis), 2 (Mvis x Mcorr)"<<endl;
    cout<<"Want HOP cut: 1 (HOP cut applied), 0 (no HOP cut)"<<endl;
  
    return 1;
  }
  

  if(argc == 15)
  {
    if(*argv[1] == '1') wantOldDataSet = true;
    
    nGenSignal = atoi(argv[2]);
    nGenPartReco = atoi(argv[3]);
    nGenComb = atoi(argv[4]);
    nGenJpsiLeak = atoi(argv[5]);
    
    int trigInt(atoi(argv[6]));
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
    
    BDTcut = argv[8];
    ntoys = atoi(argv[7]);

    if(*argv[9] == '1') constPartReco = true;
    if(*argv[10] == '2') fit2D = true;
    outputfolder = argv[11];

    if(*argv[12] == '1') wantHOPCut = true;
    if(*argv[13] == '0') wantHOPCut = false;

    minBMass = atof(argv[13]);
    maxBMass = atof(argv[14]);

    fs::path data_dir("./"+outputfolder);
    if (!fs::is_directory(data_dir)){
      cout<<"Output directory "<<outputfolder<<" does not exist."<<endl;
      return 1;
    }
    
  }
  
  if(minBMass < 4280 || maxBMass > 6280) 
  {
      cout<<"ERROR: fitting range out of tuple range. Stop."<<endl;
      return 1;
  }

  //*************** Output files

  string plotsfile, resultsfile, outfile, tablefile, workspacename;
  tablefile = outputfolder+"/resultToy.dat";

  string extraString;
  if(fit2D && wantHOPCut) extraString = "2D_HOPCut";
  if(fit2D && !wantHOPCut) extraString = "2D";
  if(!fit2D && wantHOPCut) extraString = "1D_HOPCut";
  if(!fit2D && !wantHOPCut) extraString = "1D";

  plotsfile = outputfolder+"/plotsHistBremCatTsallisBkgfit"+extraString+".root";
  resultsfile = outputfolder+"/toystudyHistBremCatTsallisBkg_results"+extraString+trigStr+".root";
  outfile = outputfolder+"/fitResult"+extraString+".dat";
  workspacename = outputfolder+"/workspace"+extraString+".root";


   //***********Get the datasets

   if(wantHOPCut) extraString = "_MH";
   if(!wantHOPCut) extraString = "";
   string fSignal("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_BDT_ctrl_trigged"+extraString+".root");
   string fPartReco("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/BJpsiX_Strip21_MC2012_ctrlNoDTF_trigged_rarebkgs"+extraString+".root");
   string fComb("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_piee_trigged"+extraString+".root");
   string fJpsiLeak("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/BJpsiX_Strip21_MC2012_signal_leakage_trigged"+extraString+".root");

   
   //RooMsgService::instance().addStream(DEBUG, Topic(Integration));
   if (!wantOldDataSet) prepare_PDFs(workspacename, trigStr, BDTcut, fit2D, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);


   //***************Prepare the stuff to generate events


   TFile f(resultsfile.c_str(),"recreate");
   TTree t(("paramsFloatingExpFloatingFracPartReco_"+trigStr).c_str(), ("paramsFloatingExpFloatingFracPartReco_"+trigStr).c_str());

   bool wantPlots(true);
   bool update(false);
   
   ofstream out(outfile);

   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Generation and fit number "<<i<<endl;

      generate_and_fit( workspacename, fit2D, wantPlots, constPartReco, 
                nGenSignal,   nGenPartReco,   nGenComb, nGenJpsiLeak,
                nGenFracZeroGamma,   nGenFracOneGamma, 0.1,
                out, &t,  update, plotsfile);


      wantPlots = false;
      update = true;
   }

   out.close();

   f.cd();
   t.Write();

   int w(10);


   //save toy characteristics and results

   ofstream out2(tablefile, std::ios::app);
   out2<<endl; for(int i(0); i<50; ++i) out2<<"="; out2<<endl;
   if (fit2D) out2<<"2D Mvis X MCorr with Tsallis bkg";
   else out2<<"1D Mvis with expo bkg";
   if(wantHOPCut) out2<<" HOP CUT:"<<endl;
   else out2<<":"<<endl;
   makeTableResults(&t, nGenSignal, nGenPartReco, nGenComb, out2 );
   out2<<endl;
   out2<<setw(w)<<"HOP CUT"<<setw(w)<<"min M"<<setw(w)<<"max M"<<setw(w)<<"Trigger"<<setw(w)<<"constrain"<<endl;
   out2<<setw(w)<<(wantHOPCut==true ? "yes":"no")<<setw(w)<<minBMass<<setw(w)<<maxBMass<<setw(w)<<trigStr.substr(0, trigStr.size()-6)<<setw(w)<<(constPartReco==true? "yes":"no")<<endl;
   out2<<endl;
   out2<<endl; for(int i(0); i<50; ++i) out2<<"="; out2<<endl;
   out2.close();


   f.Close();

   return 0;
}


