#ifndef FITTER_UTILS_1DSIMULTANEOUS_H
#define FITTER_UTILS_1DSIMULTANEOUS_H

#include "RooPTMVis.h"
#include "fitter_utils.h"
#include "RooBinning.h"
#include "RooRandom.h"
#include "RooFitResult.h"

class FitterUtils1DSimultaneous: public FitterUtils
{
   public:

      FitterUtils1DSimultaneous(int nGenKemu_, int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_, double nGenFracOneGamma_, double PPerpCut_, string workspacename_);

      void prepare_PDFs(string trigStr, string BDTVar, double BDTcut,
            string signalfile, string partrecofile, string combinatorialfile, string JpsiLeakfile,
            double minBMass = 4880, double maxBMass = 5700,
            string signaltree = "DecayTree", string partrecotree = "DecayTree", string combinatorialtree = "DecayTree", string JpsiLeaktree = "DecayTree");

      void initiateParamsKee(RooRealVar& nSignal, RooRealVar& nPartReco, RooRealVar& nComb, RooRealVar& nJpsiLeak,
            RooRealVar& fracZeroGammaLowPPerp, RooRealVar const& fracZeroGammaLowPPerpSigma, RooRealVar& fracOneGammaLowPPerp, RooRealVar const& fracOneGammaLowPPerpSigma, 
            RooRealVar& fracZeroGammaHighPPerp, RooRealVar const& fracZeroGammaHighPPerpSigma, RooRealVar& fracOneGammaHighPPerp, RooRealVar const& fracOneGammaHighPPerpSigma, 
            RooRealVar& expoConstCombLowPPerp, RooRealVar const& trueExpoConstCombLowPPerp, RooRealVar& expoConstCombHighPPerp, RooRealVar const& trueExpoConstCombHighPPerp,
            RooRealVar& fracCombLowOverTot, RooRealVar const& fracCombLowOverTotFromKemuFit);

      void initiateParamsKemu(RooRealVar& expoConstKemu, RooRealVar const& trueExpoConstKemu, RooRealVar& T, RooRealVar const& trueT, RooRealVar& n, RooRealVar const& trueN);

      void generate();

      void fit(bool wantplot, bool constPartReco,
            double fracPartReco_const,
            ofstream& out, TTree* tKee, TTree* tKemu, bool update, string plotsfile);

   void plot_fit_result_Kee(RooRealVar B_plus_M, RooCategory& category, string plotsfile, RooAbsPdf &totPdf, RooDataSet& dataGenTot);
   void plot_fit_result_Kemu(string plotsfile, RooAbsPdf &totPdf, RooDataSet& dataGenTot);

   void display();

   protected:

      double nGenFracZeroGammaLowPPerp;
      double nGenFracOneGammaLowPPerp;
      double nGenFracZeroGammaHighPPerp;
      double nGenFracOneGammaHighPPerp;


      int nGenSignalLowPPerp;
      int nGenSignalZeroGammaLowPPerp;
      int nGenSignalOneGammaLowPPerp;
      int nGenSignalTwoGammaLowPPerp;
      int nGenPartRecoLowPPerp;
      int nGenCombLowPPerp;
      int nGenJpsiLeakLowPPerp;

      int nGenSignalHighPPerp;
      int nGenSignalZeroGammaHighPPerp;
      int nGenSignalOneGammaHighPPerp;
      int nGenSignalTwoGammaHighPPerp;
      int nGenPartRecoHighPPerp;
      int nGenCombHighPPerp;
      int nGenJpsiLeakHighPPerp;

      double nGenFracSigLowOverTot;
      double nGenFracSigZeroGammaLowOverTot;
      double nGenFracSigOneGammaLowOverTot;
      double nGenFracSigTwoGammaLowOverTot;
      double nGenFracCombLowOverTot;
      double nGenFracPartRecoLowOverTot;
      double nGenFracJpsiLeakLowOverTot;

      int nGenSignalZeroGamma;
      int nGenSignalOneGamma;
      int nGenSignalTwoGamma;


      RooRealVar FracCombRealVarLowOverTot;
      double FracCombFromFitLowOverTot;
      int nCombFromFitLowPPerp;
      int nCombFromFitHighPerp;


      int nGenKemu;
      double PPerpCut;



};



#endif
