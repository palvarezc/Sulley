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
   double BDTCutVal(0.257);
   string BDTVar;
   int ntoys(1000);
   bool fit2D(true);
   bool constPartReco(0);
   string outputfolder="fit_result";
   bool wantHOPCut(false);
   double minBMass(4880);
   double maxBMass(5700);
   string paramToChange("");
   double nSigmas(0);

   if(argc != 17)
   {
      cout<<"toystudy:  Launchs a ToyMC systematic study of the fitter specified by the options"<<endl;
      cout<<"Syntax: "<<argv[0]<<" <no_reload> <nsignal> <npartreco> <ncomb> <nJpsiLeak> <trigger_cat> <ntoys> <bdt_var_name> <bdt_cut> <constPartReco> <paramToChange> <nSigmasVariation> <output_folder> <HOP cut> <minBMass> <maxBMass>"<<endl;
      cout<<endl;
      cout<<endl;
      cout<<"no_reload: 0 (build PDFs from control channel), 1 (use PDFs in workspace)"<<endl;
      cout<<"trig_cat: 0 (ETOS), 1 (HTOS), 2 (TIS)"<<endl;
      cout<<"constPartReco: Constrain part. reco. fraction? 0 (No), 1 (Yes)"<<endl;
      cout<<"ndims: 1 (Mvis), 2 (Mvis x Mcorr)"<<endl;
      cout<<"Want HOP cut: 1 (HOP cut applied), 0 (no HOP cut)"<<endl;
      cout<<"paramToChange: either T or n"<<endl;
      cout<<"nSigmasVariation: i.e. if paramToChange=T and nSigmasVariations = -1.5, T -> T-1.5*err(T)"<<endl;

      return 1;
   }

   //1 <no_reload>
   //2 <nsignal>
   //3 <npartreco>
   //4 <ncomb>
   //5 <nJpsiLeak>
   //6 <trigger_cat>
   //7 <ntoys>
   //8 <bdt_var_name>
   //9 <bdt_cut>
   //10 <constPartReco>
   //11 <paramToChange>
   //12 <nSigmasVariation>
   //13 <output_folder>
   //14 <HOP cut>
   //15 <minBMass>
   //16 <maxBMass>

   if(argc == 17)
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

      ntoys = atoi(argv[7]);

      BDTVar = argv[8];
      BDTCutVal = atof(argv[9]);


      if(*argv[10] == '1') constPartReco = true;
      paramToChange = argv[11];
      if(paramToChange != "T" && paramToChange != "n")
      {
         cout<<"ERROR: parameter to change not known. exit"<<endl;
         return 1;
    }
   nSigmas = atof(argv[12]);

    outputfolder = argv[13];

    if(*argv[14] == '1') wantHOPCut = true;
    if(*argv[14] == '0') wantHOPCut = false;

    minBMass = atof(argv[15]);
    maxBMass = atof(argv[16]);

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


  string fSignal("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/B2Kee_Strip21_MC_trigged.root");
  string fComb("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/B2Kemu_Strip21_trigged.root");
  string fJpsiLeak("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/total_signal_leakage_trigged.root");
  string fPartReco("/home/hep/th1011/B2KeeData/tuples/strip21/tupleThibaud/newTuples/trigged/BJpsiX_Strip21_MC2012_ctrl_noDTFMAllmass_trigged_rarebkgs.root");

   
   //RooMsgService::instance().addStream(DEBUG, Topic(Integration));
   

   FitterUtils fu(nGenSignal,nGenPartReco, nGenComb, nGenJpsiLeak, nGenFracZeroGamma, nGenFracOneGamma, fit2D, workspacename);

   if (!wantOldDataSet) fu.prepare_PDFs(trigStr, BDTVar, BDTCutVal, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);


   //***************Prepare the stuff to generate events


   TFile f(resultsfile.c_str(),"update");
   TTree t(("params_"+trigStr).c_str(), ("params_"+trigStr).c_str());

   bool wantPlots(true);
   bool update(false);
   
   ofstream out(outfile);

   fu.generate();

   TFile fw(workspacename.c_str(), "update");
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
   cout<<"workspace before change:"<<endl;
   workspace->Print();
   RooRealVar* param = (RooRealVar*)workspace->var(paramToChange.c_str());
   param->setVal(param->getVal()+nSigmas*param->getError());
   cout<<"workspace after change:"<<endl;
   workspace->Print();
   workspace->Write("", TObject::kOverwrite);
   cout<<"HELLO1"<<endl;
   fw.Close();
   cout<<"HELLO2"<<endl;
   f.cd();
   cout<<"HELLO3"<<endl;

   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Fit number "<<i<<endl;

      fu.fit( wantPlots, constPartReco, 0.1, out, &t,  update, plotsfile);

      wantPlots = false;
      update = true;
   }

   out.close();

   f.cd();
   t.Write();

   int w(10);


   //save toy characteristics and results

   ofstream out2(tablefile, std::ios::app);
   out2<<endl; for(int i(0); i<70; ++i) out2<<"="; out2<<endl;
   if (fit2D) out2<<"2D Mvis X MCorr with Tsallis bkg";
   else out2<<"1D Mvis with expo bkg";
   if(wantHOPCut) out2<<" HOP CUT:"<<endl;
   else out2<<":"<<endl;
   makeTableResults(&t, nGenSignal, nGenPartReco, nGenComb, nGenJpsiLeak, out2 );
   out2<<endl;
   out2<<setw(w)<<"HOP CUT"<<setw(w)<<"min M"<<setw(w)<<"max M"<<setw(w)<<"Trigger"<<setw(w)<<"constrain"<<setw(w)<<"JPsi leak"<<endl;
   out2<<setw(w)<<(wantHOPCut==true ? "yes":"no")<<setw(w)<<minBMass<<setw(w)<<maxBMass<<setw(w)<<trigStr.substr(0, trigStr.size()-6)<<setw(w)<<(constPartReco==true? "yes":"no")<<setw(w)<<nGenJpsiLeak<<endl;
   out2<<endl;
   out2<<endl; for(int i(0); i<70; ++i) out2<<"="; out2<<endl;
   out2.close();


   f.Close();

   return 0;
}


