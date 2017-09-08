#ifndef FITANDSPLOTKEEDATAFORTRAINING 
#define FITANDSPLOTKEEDATAFORTRAINING 

#include "usefulFunctions.h"
//#include "RooCBShape.h"
//#include "doubleCrystalBall.hpp"
//#include "RooTripleGaussian.hpp"
#include<string>

using namespace RooFit;
using namespace std;


class FitAndSplotKeeDataForTraining
{


   public:

      FitAndSplotKeeDataForTraining()
         :workspaceFileName("/vols/lhcb/palvare1/RK_analysis/FitWorkspaces/WorkspaceKee.root"),
         //:workspaceFileName("/home/hep/th1011/Documents/Kee/fits/workspaces/FitAndSplotKeeDataTruthInvestigationWorkspace.root"),
         fHistoYieldName("/vols/lhcb/palvare1/RK_analysis/FitWorkspaces/HistoYieldKee.root"),
         plotdir("/home/hep/palvare1/public_html/Jpsiee_fit/"),
         tupleMCdir("/vols/lhcb/palvare1/KeeTuples/MC/presel/"),
         tupleMCname("Kee_RedoCalo_signalPresel.root"),
         treeMCname("DecayTree"), 
         tupleDataDir("/vols/lhcb/th1011/KeeTuples/data/trigged/"),
         // tupleDataName("B2Kee_Strip21_data_JpsiKPresel_trigged.root"),
         tupleDataName("B2Kee_Strip21_data_JpsiKPresel_trigged.root"),
         treeDataName("DecayTree"),
         tupleRarePrcDir("/vols/lhcb/th1011/KeeTuples/MC/HltTOS/"),
          // tupleRarePrcDir("/vols/lhcb/palvare1/KeeTuples/MC/HltTOS/"), // uBDTs_Run2 = BDTu1alpha4Run1R  
         tupleRarePrcName("B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root"),
         treeRarePrcName("DecayTree"),
         tupleCharmPrcDir("/vols/lhcb/th1011/KeeTuples/MC/HltTOS/"),
          // tupleCharmPrcDir("/vols/lhcb/palvare1/KeeTuples/MC/HltTOS/"), // uBDTs_Run2 = BDTu1alpha4Run1R
         tupleCharmPrcName("B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_charmBkg.root"),
         treeCharmPrcName("DecayTree"),
         // tuplePieeDir("/vols/lhcb/th1011/KeeTuples/MC/HltTOS/"),
         // tuplePieeName("Bu_JpsiX_2012RedoCalo_presel_piee_HltTrigged_rareBkg.root"),
         tuplePieeDir("/vols/lhcb/palvare1/KeeTuples/MC/presel/"),
         tuplePieeName("KJpsi_RedoCalo_JpsiKPresel_KtoPi.root"),
         treePieeName("DecayTree"),
         tupleLbDir("/vols/lhcb/th1011/KeeTuples/MC/HltTOS/"),
         tupleLbName("Lb_JpsiX_2012RedoCalo_presel_KisKorPiorP_HltTrigged.root"),
         // tupleLbDir("/vols/lhcb/palvare1/KeeTuples/MC/HltTOS/"),
         // tupleLbName("Lb_JpsiX_2012RedoCalo_presel_KisKorPiorP_HltTrigged_rareBkg.root"),
         treeLbName("DecayTree"),
         tupleRareDir("/vols/lhcb/th1011/KeeTuples/MC/HltTOS/"),
         tupleRareName("Kee_RedoCalo_JpsiKPresel_HltTrigged.root"),
         treeRareName("DecayTree"),
         cutStringData(""),
         cutStringMC(""),
         cutStringMCPiee(""),
         weightMC(""),
         weightMCPiee(""),
         B_plus_M_min_data(5080), 
         B_plus_M_max_data(5680),
         B_plus_M_min_MC(4980), 
          B_plus_M_max_MC(5680),
          minBMass_MC(4880),
          maxBMass_MC(5700),
          minBMass_data(4880),
          maxBMass_data(5700),
          B_plus_DTF_M_cut(5150),
          turnOffWeightsSignal(false)
   {}


      string workspaceFileName;
      string fHistoYieldName;
      string plotdir;
      string tupleMCdir;
      string tupleMCname;
      string treeMCname;
      string tupleDataDir;
      string tupleDataName;
      string treeDataName;
      string tupleRarePrcDir;
      string tupleRarePrcName;
      string treeRarePrcName;
      string tupleCharmPrcDir;
      string tupleCharmPrcName;
      string treeCharmPrcName;
      string tuplePieeDir;
      string tuplePieeName;
      string treePieeName;
      string tupleLbDir;
      string tupleLbName;
      string treeLbName;
      string tupleRareDir;
      string tupleRareName;
      string treeRareName;
      string cutStringData;
      string cutStringMC;
      string cutStringMCPiee;
      string weightMC;
      string weightMCPiee;

      double B_plus_M_min_data;
      double B_plus_M_max_data;
      double B_plus_M_min_MC;
      double B_plus_M_max_MC;

  double minBMass_MC;
  double maxBMass_MC;
  double minBMass_data;
  double maxBMass_data;
  double B_plus_DTF_M_cut;

  bool turnOffWeightsSignal;
  

      void initiateHistoYield();

      void fitPrc(int trigCat, int nPhotons, string whichOne); //whichOne = "Rare" Or whichOne = "Charm"
      void fitAllPrc();

      void prepareGenericMCWorkspace(int trigCat, int nPhotons, int whichMC);

      void prepareDataWorkspace(int trigCat );
      void prepareDataWorkspaceWithPoissonOsc(int trigCat, int nPoisson);

      void prepareAllPrcMCWorkspaces();

      void prepareMCWorkspace(int trigCat, int nPhotons);
      void prepareRarePrcWorkspace(int trigCat, int nPhotons);
      void prepareCharmPrcWorkspace(int trigCat, int nPhotons);
  void preparePieeWorkspace(int trigCat, int nPhotons);
      void prepareLbWorkspace(int trigCat);

      //void prepareRareWorkspace(int trigCat, int nPhotons);

  void addMassShiftScale(int trigCat, string extracut, string label);
  


  void plotOnMC(RooAbsPdf *model, RooRealVar *var, int trigCat);
  

      void fitMC(int trigCat, int nPhotons, double mean0 = 5280, double sigma0 = 7.3, double fracSigma0 = 1.2, double al0 = 0.5, double ar0 = 0.5, double nl0 = 3.5, double nr0 = 3.5, double fracCB20 = 0.66, bool fast = false );
      void fitMCAuto(int trigCat, int nPhotons, bool fast);
  void fitPiee(int trigCat, int PhotCat, bool fast = false);
      void fitLb(int trigCat, bool fast = false);
      void fitData(int trigCat, RooRealVar& yield, bool fast = false, bool saveYield = true,  double shiftFix = -2000, double scaleFix = -2000);
      void fitData1PhotCat(int trigCat, int PhotCat, RooRealVar& yield);
      void fitRare(int trigCat, int nPhotons);
      void mergeAllTriggers();
      void computeSWeight(int trigCat);
      void makeSWeightedTree(int trigCat, string extracut, string label);

      void getFracPiee(int trigCat, RooRealVar& frac);
      //void getFracRare(int trigCat, RooRealVar& frac);

      void performFullFit(int trigCat, RooRealVar& yield, bool fast = true, bool saveYield = false,  double shiftFix = -2000, double scaleFix = -2000);
      void performFullFitUseJPsiModeConstraints(int trigCat, RooRealVar& yield, bool fast, bool saveYield, string workspaceFileName);
      void fillHistoWithYields(string nameHistoFile, string _plotdir, int trigCat, int nBins, double binTab[], string varToCutOn,  double shiftFix = -2000, double scaleFix = -2000);
      void updatePoissonHisto(string nameHistoFile, string _plotdir, int trigCat, int binToFill, int nPoisson, string varToCutOn, bool refitMC,                                
         double shiftFix = -2000, double scaleFix = -2000); 
};

void fitKeeDataPerformFitWithCut(string plotDir, int trigCat, RooRealVar& yield, string MCCut, string dataCut, string MCWeight);
void fitKeeDataFillHistoWithYields(string nameHistoFile, string plotDir, int trigCat, int nBins, double binTab[], string varToCutOn, string MCCut, string dataCut, string MCWeight);
void fitKeeDataFillHistoWithYieldsUniform(string nameHistoFile, string plotDir, int trigCat, int nBins, double binTab[], string varToCutOn, string MCCut, string dataCut, string MCWeight);

#endif
