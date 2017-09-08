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
#include"fitter_utils_ana.h"
// #include"fitter_utils_simultaneous.h"
// #include"fitter_utils_1Dsimultaneous.h"
// #include "fitter_utils_ExpOfPolyTimesX.h"
// #include"fitter_utils_simultaneous_ExpOfPolyTimesX.h"
#include "fitter_utils_HistFact.h"
// #include "fitter_utils_HistFact_1D.h"
#include "TObjectTable.h"
#include "Utils.h"

namespace fs = boost::filesystem;


using namespace std;
using namespace RooFit;
using namespace Sulley;

int main(int argc, char* argv[])
{

   //***********Number of events to generata

  Options opts = Options("mytoystudy");

  opts.nGenSignal = 203;
  opts.nGenPartReco = 79;
  opts.nGenComb = 31;
  opts.nGenJpsiLeak = 10;
  opts.nGenKemu = 200;
  opts.nGenSignal_mumu = 353;
  opts.nGenBkg_mumu = 95;

  opts.nGenFracZeroGamma = 0.368;
  opts.nGenFracOneGamma = 0.484;
  opts.trigCat = -1;
  opts.trigStr = "B_plus_ETA"; // no requirement (trigStr>0.9)
  opts.weightStr = "PIDTrigDataCondWeight_ETOSOnly";
  opts.BDTcut = 0.257;
  opts.minBMass = 4880;
  opts.maxBMass = 5700;
  opts.minBMass_mumu = 5180;
  opts.maxBMass_mumu = 5700;
  opts.fitSimultaneous = false;
  opts.fittype = RooFit1D; 
  opts.interpolate = 2;
  opts.massshift = 0;
  
  // opts.defaultbinningname = "default";

  int fitsample = 0; // ee = 0,  mumu = 1,  both = 2

  bool wantOldDataSet(false);
  int ntoys(1000);
  string outputfolder = "fit_result";
  int fitMode(1);
  double PPerpCut(600);
  bool wantHOPCut(false);


   if(argc != 17 && argc != 18 && argc!=19)
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


   if(argc == 17 || argc == 18 || argc == 19)
   {
      if(*argv[1] == '1') wantOldDataSet = true;

      opts.nGenSignal = atoi(argv[2]);
      opts.nGenPartReco = atoi(argv[3]);
      opts.nGenComb = atoi(argv[4]);
      opts.nGenJpsiLeak = atoi(argv[5]);

      int trigInt(atoi(argv[6]));
      opts.trigCat = trigInt;
      
      if(trigInt == 0) 
      {
        opts.weightStr = "FidelCastro_ETOS";
        opts.trigStr = "passTrigCat0";
        opts.nGenFracZeroGamma = 0.368;
        opts.nGenFracOneGamma = 0.484;
      }      
      if(trigInt == 1) 
      {
        opts.weightStr = "FidelCastro_HTOS";
        opts.trigStr = "passTrigCat1";
        opts.nGenFracZeroGamma = 0.322;
        opts.nGenFracOneGamma = 0.484;
      }      
      if(trigInt == 2) 
      {
        opts.weightStr = "FidelCastro_TIS";
        opts.trigStr = "passTrigCat2";
        opts.nGenFracZeroGamma = 0.330;
        opts.nGenFracOneGamma = 0.495;
      }

      ntoys = atoi(argv[7]);

      opts.BDTVar = argv[8];
      opts.BDTcut = atof(argv[9]);


      if(*argv[10] == '1') opts.constPartReco = true;

      fitMode = atoi(argv[11]);

      if(fitMode == 2) opts.fittype = RooFit2D_RooPTMVis;
      if(fitMode == 3) {
	opts.fittype = RooFit2D_RooPTMVis;
	opts.fitSimultaneous = true;
      }
      if(fitMode == 5) opts.fittype = RooFit2D;
      if(fitMode == 6) {
	opts.fittype = RooFit2D;
	opts.fitSimultaneous = true;
      }
      if(fitMode == 7) opts.fittype = HistFact2D;
      if(fitMode == 8) opts.fittype = HistFact2D;// WRONG!
      

      if(fitMode == 2 or fitMode == 3 or fitMode == 6 or fitMode == 5 or fitMode==7) opts.fit2D = true;
      outputfolder = argv[12];

      if(*argv[13] == '1') wantHOPCut = true;
      if(*argv[13] == '0') wantHOPCut = false;

      opts.minBMass = atof(argv[14]);
      opts.maxBMass = atof(argv[15]);

      fs::path data_dir("./"+outputfolder);
      if (!fs::is_directory(data_dir)){
         cout<<"Output directory "<<outputfolder<<" does not exist."<<endl;
         return 1;
      }

      opts.Run = atoi(argv[16]);

   }

   if(argc == 17 && (fitMode == 3 or fitMode == 4 or fitMode == 6) )
   {
      cerr<<"ERROR : fit mode set to simultaneous but no Kemu yield specified. Exit"<<endl;
      return 1;
   }

   if(argc < 19 && fitMode == 4)
   {
      cerr<<"ERROR : fit mode set to binned simultaneous but no mis PT cut set. Exit."<<endl;
      return 1;
   }

   if( (argc == 18 || argc == 19)  && (fitMode == 3 or fitMode == 4 or fitMode == 6) ) 
     {
       opts.nGenKemu = atoi(argv[17]);
       opts.fitSimultaneous = true;
     }

   if(argc == 19 && fitMode ==4) PPerpCut = atof(argv[18]);

   if(opts.minBMass < 4280 || opts.maxBMass > 6280) 
   {
      cout<<"ERROR: fitting range out of tuple range. Stop."<<endl;
      return 1;
   }

   //*************** Output files

   string plotsfile, resultsfile, outfile, tablefile, workspacename;
   tablefile = outputfolder+"/resultToy.dat";

   string extraString("Mode"+i2s(fitMode));
   //  if(opts.fit2D && wantHOPCut) extraString = "2D_HOPCut";
   //  if(opts.fit2D && !wantHOPCut) extraString = "2D";
   //  if(!opts.fit2D && wantHOPCut) extraString = "1D_HOPCut";
   //  if(!opts.fit2D && !wantHOPCut) extraString = "1D";

   plotsfile = outputfolder+"/plotsHistBremCatTsallisBkgfit"+extraString+".root";
   resultsfile = outputfolder+"/toystudyHistBremCatTsallisBkg_results"+extraString+opts.trigStr+".root";
   outfile = outputfolder+"/fitResult"+extraString+".dat";
   workspacename = outputfolder+"/workspace"+extraString+".root";

   opts.plotsfile = plotsfile;
   opts.plotsdir = outputfolder;
   opts.workspacefile = workspacename;

   //***********Get the datasets
   if(wantHOPCut) extraString = "_MH";
   if(!wantHOPCut) extraString = "";

   if(opts.Run==1)
   {  
     // opts.signalfile_mumu = "/vols/lhcb/th1011/KmumuTuples/MC/oldTuples/HltTOS/KJpsi_mumu_HltTOS.root";
     // opts.signalfile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/Kee_RedoCalo_signalPresel_HltTrigged.root";
     // // opts.combfile = "/vols/lhcb/palvare1/B2Kee_Kmumu/Marcin/B2XMuMu_Line_Strip21_Kemu_presel.root";
     // // opts.combfile = "/vols/lhcb/th1011/extraData/KemuForToys/B2XMuMu_Line_Strip21_Kemu_presel_prepared.root";
     // opts.combfile = "/vols/lhcb/palvare1/RK_analysis/TuplesFit/B2XMuMu_Line_Strip21_Kemu_presel_prepared.root";
     // opts.JpsiLeakfile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_charmBkg.root";
     // opts.partrecofile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root";
     

     opts.signalfile_mumu = "/vols/lhcb/th1011/KmumuTuples/MC/oldTuples/HltTOS/KJpsi_mumu_HltTOS.root";
     opts.signalfile = "/vols/lhcb/th1011/KeeTuples/MC/presel/Kee_RedoCalo_signalPresel.root";
     opts.combfile = "/vols/lhcb/th1011/KeeTuples/data/trigged/B2Kemu_Strip21_trainingPresel_trigged.root";
     // opts.JpsiLeakfile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_charmBkg.root";
     // opts.partrecofile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root";
     opts.JpsiLeakfile = "/vols/lhcb/th1011/KeeTuples/MC/presel/KJpsi_RedoCalo_signalPresel.root";
     // opts.partrecofile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root";
     opts.partrecofile = "/vols/lhcb/palvare1/KeeTuples/MC/presel/PartReco_signalPresel_Kmixture.root";



   }
   else
   {
     opts.signalfile_mumu = "/vols/lhcb/th1011/KmumuTuples/MC/oldTuples/HltTOS/KJpsi_mumu_HltTOS.root";
     opts.signalfile = "/vols/lhcb/th1011/KeeTuples/MC/presel/Kee_run2_signalPresel.root";
     opts.combfile = "/vols/lhcb/th1011/KeeTuples/data/trigged/B2Kemu_run2_trainingPresel_trigged.root";
     // opts.JpsiLeakfile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_charmBkg.root";
     // opts.partrecofile = "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root";

     // //uBDTs_Run2 = BDTu1alpha4Run1R       
     // opts.JpsiLeakfile = "/vols/lhcb/palvare1/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_charmBkg.root";
     // opts.partrecofile = "/vols/lhcb/palvare1/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root";
     // opts.partrecofile = "/vols/lhcb/palvare1/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root";

     opts.partrecofile = "/vols/lhcb/palvare1/KeeTuples/MC/presel/PartReco_signalPresel_Kmixture.root";
     opts.JpsiLeakfile = "/vols/lhcb/th1011/KeeTuples/MC/presel/KJpsi_run2_signalPresel.root";
   }

   if (wantHOPCut)
   {
     opts.signalfile = "/vols/lhcb/palvare1/RK_analysis/TuplesFit/Kee_RedoCalo_signalPresel_HltTrigged_MH.root";  
     opts.combfile = "/vols/lhcb/palvare1/RK_analysis/TuplesFit/B2XMuMu_Line_Strip21_Kemu_presel_prepared_MH.root";
     opts.JpsiLeakfile = "/vols/lhcb/palvare1/RK_analysis/TuplesFit/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_charmBkg_MH.root";
     opts.partrecofile = "/vols/lhcb/palvare1/RK_analysis/TuplesFit/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg_MH.root";
     
     opts.nGenComb = floor(70171./100740.*opts.nGenComb);
   }
   

   FitterUtilsAna fu("fitterRoofit",opts);
   FitterUtilsHistFact fuHistFact("fitterHistFact", opts);
   // FitterUtilsHistFact1D fuHistFact1D(workspacename, opts);

   fu.enablePlotting();
   fuHistFact.enablePlotting();


   if (!wantOldDataSet)
   { 
     if(fitMode <= 6) fu.prepare_PDFs(fitsample);
     if(fitMode == 7) fuHistFact.prepare_PDFs(fitsample);
     // if(fitMode == 8) fuHistFact1D.prepare_PDFs();
   }

   cout<<"PDFs prepared!"<<endl;

   //***************Prepare the stuff to generate events
   TFile f(resultsfile.c_str(),"update");
   TTree t(("params_"+opts.trigStr).c_str(), ("params_"+opts.trigStr).c_str());
   TTree tKemu("tKemu", "tKemu");

   bool update(false);
   ofstream out(outfile);

   // fuHistFact1D.enablePlotting();
   
   for(int i(0); i<ntoys; ++i)
   {
      cout<<endl<<endl<<endl<<endl<<endl;
      cout<<"************************************"<<endl<<"*********** NEW FIT *****************"<<endl<<"*************************"<<endl;
      cout<<"Generation and fit number "<<i<<endl;

      if(fitMode <=6)
      {
        fu.run_toy( 0.1, out, &t,  update, fitsample);
      }
      if(fitMode == 7)
      {
        fuHistFact.run_toy( 0.1, out, &t,  update);
      }
      // if(fitMode == 8)
      // {
      //   fuHistFact1D.generate();
      //   fuHistFact1D.fit( 0.1, out, &t,  update);
      // }

      fu.disablePlotting();
      fuHistFact.disablePlotting();
      // fuHistFact1D.disablePlotting();
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
   if (opts.fit2D) out2<<"2D Mvis X MisPT";
   else out2<<"1D Mvis with expo bkg";
   if(wantHOPCut) out2<<" HOP CUT:"<<endl;
   else out2<<":"<<endl;
   // makeTableResults(&t, opts.nGenSignal, opts.nGenPartReco, opts.nGenComb, opts.nGenJpsiLeak, out2 );
   out2<<endl;
   out2<<setw(w)<<"HOP CUT"<<setw(w)<<"min M"<<setw(w)<<"max M"<<setw(w)<<"Trigger"<<setw(w)<<"constrain"<<setw(w)<<"JPsi leak"<<endl;
   out2<<setw(w)<<(wantHOPCut==true ? "yes":"no")<<setw(w)<<opts.minBMass<<setw(w)<<opts.maxBMass<<setw(w)<<opts.trigStr.substr(0, opts.trigStr.size()-6)<<setw(w)<<(opts.constPartReco==true? "yes":"no")<<setw(w)<<opts.nGenJpsiLeak<<endl;
   out2<<endl;
   out2<<endl; for(int i(0); i<70; ++i) out2<<"="; out2<<endl;
   out2.close();


   f.Close();

   return 0;
}


