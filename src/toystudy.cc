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
#include"fitter_utils_simultaneous.h"
#include"fitter_utils_1Dsimultaneous.h"
#include "fitter_utils_ExpOfPolyTimesX.h"
#include"fitter_utils_simultaneous_ExpOfPolyTimesX.h"

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
   int nGenKemu(200);

   double nGenFracZeroGamma(0.368);
   double nGenFracOneGamma(0.484);


   //*********** Get arguments and set stuff
   string trigStr("L0ETOSOnly_d");
   bool wantOldDataSet(false);
   double BDTCutVal(0.257);
   string BDTVar;
   int ntoys(1000);
   bool fit2D(0);
   bool constPartReco(0);
   string outputfolder="fit_result";
   bool wantHOPCut(false);
   double minBMass(4880);
   double maxBMass(5700);
   int fitMode(1);
   double PPerpCut(600);

   if(argc != 16 && argc != 17 && argc!=18)
   {
      cout<<"toystudy:  Launchs a ToyMC study of the fitter specified by the options"<<endl;
      cout<<"Syntax: "<<argv[0]<<" <no_reload> <nsignal> <npartreco> <ncomb> <nJpsiLeak> <trigger_cat> <ntoys> <bdt_var_name> <bdt_cut> <constPartReco> <mode> <output_folder> <HOP cut> <minBMass> <maxBMass> <nKemu (optional)> <misPTcut (optional)"<<endl;
      cout<<endl;
      cout<<endl;
      cout<<"no_reload: 0 (build PDFs from control channel), 1 (use PDFs in workspace)"<<endl;
      cout<<"trig_cat: 0 (ETOS), 1 (HTOS), 2 (TIS)"<<endl;
      cout<<"constPartReco: Constrain part. reco. fraction? 0 (No), 1 (Yes)"<<endl;
      cout<<"mode: 1 (1D fit Mvis), 2 (2D fit Mvis x misPT), 3 (2D fit simultaneous with Kemu), 4 (binned), 5 (2D fit Mvis x misPT, RooPolyTimesX fit function),..."<<endl;
      cout<<"... 6 (2D fit Mvis x misPT simultaneous with Kemu, RooPolyTimesX fit function)"<<endl;
      cout<<"Want HOP cut: 1 (HOP cut applied), 0 (no HOP cut)"<<endl;
      cout<<"nKemu: OPTIONAL, if mode = 3, 4 or 6, must specify number of Kemu events"<<endl;
      cout<<"misPT: OPTIONAL, if mode = 4, must put PT cut"<<endl;

      return 1;
   }


   if(argc == 16 || argc == 17 || argc == 18)
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

      fitMode = atoi(argv[11]);
      if(fitMode == 2 or fitMode == 3 or fitMode == 6 or fitMode == 5) fit2D = true;
      outputfolder = argv[12];

      if(*argv[13] == '1') wantHOPCut = true;
      if(*argv[13] == '0') wantHOPCut = false;

      minBMass = atof(argv[14]);
      maxBMass = atof(argv[15]);

      fs::path data_dir("./"+outputfolder);
      if (!fs::is_directory(data_dir)){
         cout<<"Output directory "<<outputfolder<<" does not exist."<<endl;
         return 1;
      }

   }

   if(argc == 16 && (fitMode == 3 or fitMode == 4 or fitMode == 6) )
   {
      cerr<<"ERROR : fit mode set to simultaneous but no Kemu yield specified. Exit"<<endl;
      return 1;
   }

   if(argc < 18 && fitMode == 4)
   {
      cerr<<"ERROR : fit mode set to binned simultaneous but no mis PT cut set. Exit."<<endl;
      return 1;
   }

   if( (argc == 17 || argc == 18)  && (fitMode == 3 or fitMode == 4 or fitMode == 6) ) nGenKemu = atoi(argv[16]);
   if(argc == 18 && fitMode ==4) PPerpCut = atof(argv[17]);

   if(minBMass < 4280 || maxBMass > 6280) 
   {
      cout<<"ERROR: fitting range out of tuple range. Stop."<<endl;
      return 1;
   }

   //*************** Output files

   string plotsfile, resultsfile, outfile, tablefile, workspacename;
   tablefile = outputfolder+"/resultToy.dat";

   string extraString("Mode"+i2s(fitMode));
   //  if(fit2D && wantHOPCut) extraString = "2D_HOPCut";
   //  if(fit2D && !wantHOPCut) extraString = "2D";
   //  if(!fit2D && wantHOPCut) extraString = "1D_HOPCut";
   //  if(!fit2D && !wantHOPCut) extraString = "1D";

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



   //   string fSignal("/vols/lhcbdesk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_BDT_ctrl_trigged"+extraString+".root");
   //   string fPartReco("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/BJpsiX_Strip21_MC2012_ctrlNoDTF_trigged_rarebkgs"+extraString+".root");
   // //  string fPartReco("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_BDT_prc_trigged.root");
   //  // string fComb("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kee_Strip21_piee_trigged"+extraString+".root");
   //   string fComb("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/B2Kemu_Strip21_trigged.root");
   //   string fJpsiLeak("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/total_signal_leakage_trigged"+extraString+".root");
   //
   ////   string fSignal("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldStrippingTrees/B2Kee_Strip21_BDT_ctrl_trigged.root");
   ////   string fPartReco("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldStrippingTrees/B2Kee_Strip21_BDT_prc_trigged.root");
   ////   string fComb("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/oldStrippingTrees/B2Kee_Strip21_piee_trigged.root");
   ////   string fJpsiLeak("/vols/lhcbdisk04/thibaud/tuples/B2Kee/tuples/strip21/tupleThibaud/total_signal_leakage_trigged.root");


   //RooMsgService::instance().addStream(DEBUG, Topic(Integration));



   FitterUtils fu(nGenSignal,nGenPartReco, nGenComb, nGenJpsiLeak, nGenFracZeroGamma, nGenFracOneGamma, fit2D, workspacename);
   FitterUtilsSimultaneous fuS(nGenKemu, nGenSignal,nGenPartReco, nGenComb, nGenJpsiLeak, nGenFracZeroGamma, nGenFracOneGamma, workspacename);
   FitterUtils1DSimultaneous fuS1D(nGenKemu, nGenSignal,nGenPartReco, nGenComb, nGenJpsiLeak, nGenFracZeroGamma, nGenFracOneGamma, PPerpCut, workspacename);
   FitterUtilsExpOfPolyTimesX fuExpOfPoly(nGenSignal,nGenPartReco, nGenComb, nGenJpsiLeak, nGenFracZeroGamma, nGenFracOneGamma, workspacename);
   FitterUtilsSimultaneousExpOfPolyTimesX fuSExpOfPoly(nGenKemu, nGenSignal,nGenPartReco, nGenComb, nGenJpsiLeak, nGenFracZeroGamma, nGenFracOneGamma, workspacename);



   if (!wantOldDataSet)
   { 
      if(fitMode <= 2) fu.prepare_PDFs(trigStr, BDTVar, BDTCutVal, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);
      if(fitMode == 3) fuS.prepare_PDFs(trigStr, BDTVar, BDTCutVal, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);
      if(fitMode == 4) fuS1D.prepare_PDFs(trigStr, BDTVar, BDTCutVal, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);
      if(fitMode == 5) fuExpOfPoly.prepare_PDFs(trigStr, BDTVar, BDTCutVal, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);
      if(fitMode == 6) fuSExpOfPoly.prepare_PDFs(trigStr, BDTVar, BDTCutVal, fSignal, fPartReco, fComb, fJpsiLeak, minBMass, maxBMass);
   }



   //***************Prepare the stuff to generate events


   TFile f(resultsfile.c_str(),"update");
   TTree t(("params_"+trigStr).c_str(), ("params_"+trigStr).c_str());
   TTree tKemu("tKemu", "tKemu");

   bool wantPlots(true);
   bool update(false);

   ofstream out(outfile);


   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Generation and fit number "<<i<<endl;

      if(fitMode <=2)
      {
         fu.generate();
         fu.fit( wantPlots, constPartReco, 0.1, out, &t,  update, plotsfile);
      }
      if(fitMode == 3)
      {
         fuS.generate();
         fuS.fit( wantPlots, constPartReco, 0.1, out, &t,  update, plotsfile);
      }
      if(fitMode == 4)
      {
         fuS1D.generate();
         fuS1D.fit( wantPlots, constPartReco, 0.1, out, &t, &tKemu,  update, plotsfile);
      }
      if(fitMode == 5)
      {
         fuExpOfPoly.generate(wantPlots, plotsfile);
         fuExpOfPoly.fit( wantPlots, constPartReco, 0.1, out, &t,  update, plotsfile);
      }
      if(fitMode == 6)
      {
         fuSExpOfPoly.generate(wantPlots, plotsfile);
         fuSExpOfPoly.fit( wantPlots, constPartReco, 0.1, out, &t,  update, plotsfile);
      }

      wantPlots = false;
      update = true;
   }

   out.close();

   f.cd();
   t.Write();
   if(fitMode == 4) tKemu.Write();

   int w(10);


   //save toy characteristics and results

   ofstream out2(tablefile, std::ios::app);
   out2<<endl; for(int i(0); i<70; ++i) out2<<"="; out2<<endl;
   if (fit2D) out2<<"2D Mvis X MisPT";
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


