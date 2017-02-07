#include "fitter_utils_simultaneous_ExpOfPolyTimesX.h"
#include "RooExpOfPolyTimesX.h"
#include "RooBinning.h"
#include "RooRandom.h"


FitterUtilsSimultaneousExpOfPolyTimesX::FitterUtilsSimultaneousExpOfPolyTimesX(int nGenKemu_, int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_,
    double nGenFracZeroGamma_, double nGenFracOneGamma_, string workspacename_)
   :FitterUtilsExpOfPolyTimesX(nGenSignal_, nGenPartReco_, nGenComb_, nGenJpsiLeak_, nGenFracZeroGamma_, nGenFracOneGamma_, workspacename_), nGenKemu(nGenKemu_)
   {}


void FitterUtilsSimultaneousExpOfPolyTimesX::initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma,
      RooRealVar& nKemu, RooRealVar& nSignal, RooRealVar& nPartReco,
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne,
      RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma,
      RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee,
      RooRealVar& l1Kemu, RooRealVar& l2Kemu, RooRealVar& l3Kemu, RooRealVar& l4Kemu, RooRealVar& l5Kemu,
      RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen 

      )
{
  FitterUtilsExpOfPolyTimesX::initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma,
           nSignal, nPartReco, nComb, fracZero, fracOne, nJpsiLeak, constPartReco, fracPartRecoSigma,
           l1Kee, l2Kee, l3Kee, l4Kee, l5Kee,
                                             l1KeeGen, l2KeeGen, l3KeeGen, l4KeeGen, l5KeeGen ,0,0); 



  TRandom rand;
  rand.SetSeed();

  nKemu.setVal(rand.Uniform(nGenKemu-5*sqrt(nGenKemu), nGenKemu+5*sqrt(nGenKemu)));
  nKemu.setRange(nGenKemu-10*sqrt(nGenKemu), nGenKemu+10*sqrt(nGenKemu));

  l1Kemu.setVal(rand.Uniform( l1KeeGen.getVal() - 5*l1KeeGen.getError(), l1KeeGen.getVal() + 5*l1KeeGen.getError() ) );
  l1Kemu.setRange( l1KeeGen.getVal() - 10*l1KeeGen.getError(), l1KeeGen.getVal() + 10*l1KeeGen.getError() );

  l2Kemu.setVal(rand.Uniform( l2KeeGen.getVal() - 5*l2KeeGen.getError(), l2KeeGen.getVal() + 5*l2KeeGen.getError() ) );
  l2Kemu.setRange( l2KeeGen.getVal() - 10*l2KeeGen.getError(), l2KeeGen.getVal() + 10*l2KeeGen.getError() );

  l3Kemu.setVal(rand.Uniform( l3KeeGen.getVal() - 5*l3KeeGen.getError(), l3KeeGen.getVal() + 5*l3KeeGen.getError() ) );
  l3Kemu.setRange( l3KeeGen.getVal() - 10*l3KeeGen.getError(), l3KeeGen.getVal() + 10*l3KeeGen.getError() );

  l4Kemu.setVal(rand.Uniform( l4KeeGen.getVal() - 5*l4KeeGen.getError(), l4KeeGen.getVal() + 5*l4KeeGen.getError() ) );
  l4Kemu.setRange( l4KeeGen.getVal() - 10*l4KeeGen.getError(), l4KeeGen.getVal() + 10*l4KeeGen.getError() );

  l5Kemu.setVal(rand.Uniform( l5KeeGen.getVal() - 5*l5KeeGen.getError(), l5KeeGen.getVal() + 5*l5KeeGen.getError() ) );
  l5Kemu.setRange( l5KeeGen.getVal() - 10*l5KeeGen.getError(), l5KeeGen.getVal() + 10*l5KeeGen.getError() );

}

void FitterUtilsSimultaneousExpOfPolyTimesX::generate(bool wantPlots, string plotsfile)
{
   FitterUtilsExpOfPolyTimesX::generate(wantPlots, plotsfile);
   TFile fw(workspacename.c_str(), "UPDATE");
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");

   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooDataSet* dataSetCombExt = (RooDataSet*)workspace->data("dataSetCombExt");
   RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
//   RooRealVar *l1KeeGen = workspace->var("l1KeeGen");  
//   RooRealVar *l2KeeGen = workspace->var("l2KeeGen");  
//   RooRealVar *l3KeeGen = workspace->var("l3KeeGen");  
//   RooRealVar *l4KeeGen = workspace->var("l4KeeGen");  
//   RooRealVar *l5KeeGen = workspace->var("l5KeeGen");  
//
//
//   RooExpOfPolyTimesX kemuPDF("kemuPDF", "kemuPDF",  *B_plus_M, *misPT,  *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);
//
//   RooAbsPdf::GenSpec* GenSpecKemu = kemuPDF.prepareMultiGen(RooArgSet(*B_plus_M, *misPT), RooFit::Extended(1), NumEvents(nGenKemu));
//
//   cout<<"Generating Kemu"<<endl;
//   RooDataSet* dataGenKemu = kemuPDF.generate(*GenSpecKemu);//(argset, 100, false, true, "", false, true);
//   dataGenKemu->SetName("dataGenKemu"); dataGenKemu->SetTitle("dataGenKemu");
//
//
//   RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
//   workspaceGen->import(*dataGenKemu);
//
//   workspaceGen->Write("", TObject::kOverwrite);
//   fw.Close();
//   delete dataGenKemu;
//   delete GenSpecKemu;


   TVectorD rho(2);
   rho[0] = 2.5;
   rho[1] = 1.5;
   misPT->setRange(-2000, 5000);
   RooNDKeysPdf kemuPDF("kemuPDF", "kemuPDF", RooArgList(*B_plus_M, *misPT), *dataSetCombExt, rho, "ma",3, true);
   misPT->setRange(0, 5000);


   RooAbsPdf::GenSpec* GenSpecKemu = kemuPDF.prepareMultiGen(RooArgSet(*B_plus_M, *misPT), RooFit::Extended(1), NumEvents(nGenKemu));

   cout<<"Generating Kemu"<<endl;
   RooDataSet* dataGenKemu = kemuPDF.generate(*GenSpecKemu);//(argset, 100, false, true, "", false, true);
   dataGenKemu->SetName("dataGenKemu"); dataGenKemu->SetTitle("dataGenKemu");


   RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
   workspaceGen->import(*dataGenKemu);

   if(wantPlots) PlotShape(*dataSetComb, *dataGenKemu, kemuPDF, plotsfile, "cKemuKeys", *B_plus_M, *misPT);

   fw.cd();
   workspaceGen->Write("", TObject::kOverwrite);
   fw.Close();
   delete dataGenKemu;
   delete GenSpecKemu;
}


void FitterUtilsSimultaneousExpOfPolyTimesX::fit(bool wantplot, bool constPartReco,
      double fracPartReco_const,
      ofstream& out, TTree* t, bool update, string plotsfile)
{

   //***************Get the PDFs from the workspace

   TFile fw(workspacename.c_str());   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooRealVar *l1Kee = workspace->var("l1Kee");
   RooRealVar *l2Kee = workspace->var("l2Kee");
   RooRealVar *l3Kee = workspace->var("l3Kee");
   RooRealVar *l4Kee = workspace->var("l4Kee");
   RooRealVar *l5Kee = workspace->var("l5Kee");
   RooRealVar *l1KeeGen = workspace->var("l1KeeGen");
   RooRealVar *l2KeeGen = workspace->var("l2KeeGen");
   RooRealVar *l3KeeGen = workspace->var("l3KeeGen");
   RooRealVar *l4KeeGen = workspace->var("l4KeeGen");
   RooRealVar *l5KeeGen = workspace->var("l5KeeGen");
   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");


   RooRealVar l1Kemu(*l1Kee);
   l1Kemu.SetName("l1Kemu"); l1Kemu.SetTitle("l1Kemu");    
   RooRealVar l2Kemu(*l2Kee);
   l2Kemu.SetName("l2Kemu"); l2Kemu.SetTitle("l2Kemu");    
   RooRealVar l3Kemu(*l3Kee);
   l3Kemu.SetName("l3Kemu"); l3Kemu.SetTitle("l3Kemu");    
   RooRealVar l4Kemu(*l4Kee);
   l4Kemu.SetName("l4Kemu"); l4Kemu.SetTitle("l4Kemu");    
   RooRealVar l5Kemu(*l5Kee);
   l5Kemu.SetName("l5Kemu"); l5Kemu.SetTitle("l5Kemu");    


   RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
   RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
   RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
   RooHistPdf *histPdfJpsiLeak(0);
   if(nGenJpsiLeak>0) histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");

   //Here set in the Kemu PDF the parameters that have to be shared

   RooExpOfPolyTimesX* combPDF = new RooExpOfPolyTimesX("combPDF", "combPDF",  *B_plus_M, *misPT,  *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee);
   RooExpOfPolyTimesX* KemuPDF = new RooExpOfPolyTimesX("kemuPDF", "kemuPDF",  *B_plus_M, *misPT,  l1Kemu, *l2Kee, *l3Kee, *l4Kee, *l5Kee);



   RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
   RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGamma");
   RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen->data("dataGenSignalOneGamma");
   RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGamma");
   RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen->data("dataGenPartReco");
   RooDataSet* dataGenComb = (RooDataSet*)workspaceGen->data("dataGenComb");
   RooDataSet* dataGenKemu = (RooDataSet*)workspaceGen->data("dataGenKemu");
   RooDataSet* dataGenJpsiLeak(0);
   if(nGenJpsiLeak>0) dataGenJpsiLeak = (RooDataSet*)workspaceGen->data("dataGenJpsiLeak");


   if(wantplot)
   {
      //**************Must get the datasets

      RooDataSet* dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
      RooDataSet* dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
      RooDataSet* dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
      RooDataSet* dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
      RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
      RooDataSet* dataSetJpsiLeak = (RooDataSet*)workspace->data("dataSetJpsiLeak");

      //**************Plot all the different components

      cout<<"dataGenSignalZeroGamma: "<<dataGenSignalZeroGamma<<endl;
      PlotShape(*dataSetSignalZeroGamma, *dataGenSignalZeroGamma, *histPdfSignalZeroGamma, plotsfile, "cSignalZeroGamma", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalOneGamma, *dataGenSignalOneGamma, *histPdfSignalOneGamma, plotsfile, "cSignalOneGamma", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalTwoGamma, *dataGenSignalTwoGamma, *histPdfSignalTwoGamma, plotsfile, "cSignalTwoGamma", *B_plus_M, *misPT);
      PlotShape(*dataSetPartReco, *dataGenPartReco, *histPdfPartReco, plotsfile, "cPartReco", *B_plus_M, *misPT);
      PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cComb", *B_plus_M, *misPT);
      if(nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeak, *dataGenJpsiLeak, *histPdfJpsiLeak, plotsfile, "cJpsiLeak", *B_plus_M, *misPT);
   }

   //***************Merge datasets

   RooDataSet* dataGenTot(dataGenPartReco);
   dataGenTot->append(*dataGenSignalZeroGamma);
   dataGenTot->append(*dataGenSignalOneGamma);
   dataGenTot->append(*dataGenSignalTwoGamma);
   dataGenTot->append(*dataGenComb);
   if(nGenJpsiLeak>0) dataGenTot->append(*dataGenJpsiLeak);

   //**************Create index category and join samples

   RooCategory category("category", "category");
   category.defineType("Kee");
   category.defineType("Kemu");

   RooDataSet dataGenSimultaneous("dataGenSimultaneous", "dataGenSimultaneous", RooArgSet(*B_plus_M, *misPT), Index(category), Import("Kee", *dataGenTot), Import("Kemu", *dataGenKemu));

   //**************Prepare fitting function

   RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   RooRealVar nKemu("nKemu", "#nKemu", 1.*nGenKemu, nGenKemu-7*sqrt(nGenKemu), nGenKemu+7*sqrt(nGenKemu));
   RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
   RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);
   RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(nPartReco,nSignal));
   RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));

   RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", RooArgList(*histPdfSignalZeroGamma, *histPdfSignalOneGamma, *histPdfSignalTwoGamma), RooArgList(fracZero, fracOneRec));

   RooArgList pdfList(histPdfSignal, *histPdfPartReco, *combPDF);
   RooArgList yieldList(nSignal, nPartReco, nComb);

   if(nGenJpsiLeak>0)
   {
      pdfList.add(*histPdfJpsiLeak);
      yieldList.add(nJpsiLeak); 
   }
   RooAddPdf totPdf("totPdf", "totPdf", pdfList, yieldList);
   RooExtendPdf totKemuPdf("totKemuPdf", "totKemuPdf", *KemuPDF, nKemu);

   //**************** Prepare simultaneous PDF

   RooSimultaneous simPdf("simPdf", "simPdf", category);
   simPdf.addPdf(totPdf, "Kee");
   simPdf.addPdf(totKemuPdf, "Kemu");

   //**************** Constrain the fraction of zero and one photon

   int nGenSignalZeroGamma(floor(nGenFracZeroGamma*nGenSignal));
   int nGenSignalOneGamma(floor(nGenFracOneGamma*nGenSignal));
   int nGenSignalTwoGamma(floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));

   RooRealVar fracZeroConstMean("fracZeroConstMean", "fracZeroConstMean", nGenSignalZeroGamma*1./nGenSignal);
   RooRealVar fracZeroConstSigma("fracZeroConstSigma", "fracZeroConstSigma", sqrt(nGenSignalZeroGamma)/nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroConstMean, fracZeroConstSigma); 

   RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", nGenSignalOneGamma*1./nGenSignal/(1-fracZeroConstMean.getVal()));
   RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", sqrt(nGenSignalOneGamma)/nGenSignal/(1-fracZeroConstMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneConstMean, fracOneConstSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", nGenPartReco/(1.*nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());

   RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

   RooRealVar JpsiLeakMean("JpsiLeakMean", "JpsiLeakMean", nGenJpsiLeak);
   RooRealVar JpsiLeakSigma("JpsiLeakSigma", "JpsiLeakSigma", nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
   RooGaussian JpsiLeakConst("JpsiLeakConst", "JpsiLeakConst", nJpsiLeak, JpsiLeakMean, JpsiLeakSigma); 


   //**************** fit
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;


   initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
         nKemu, nSignal, nPartReco, nComb, fracZero, fracOne,
         nJpsiLeak, constPartReco, fracPartRecoSigma, 
         *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, l1Kemu, l2Kemu, l3Kemu, l4Kemu, l5Kemu, 
         *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);

   RooArgSet constraints(fracZeroConst, fracOneConst);
   if (constPartReco) constraints.add(fracPartRecoConst);
   if(nGenJpsiLeak>0) constraints.add(JpsiLeakConst);

   RooAbsReal* nll = simPdf.createNLL(dataGenSimultaneous, Extended(), ExternalConstraints(constraints));
   RooMinuit minuit(*nll);
   minuit.setStrategy(2);


   int migradRes(1);
   int hesseRes(4);

   vector<int> migradResVec;
   vector<int> hesseResVec;

   double edm(10);
   int nrefit(0);

   RooFitResult* fitRes(0);
   vector<RooFitResult*> fitResVec;

   bool hasConverged(false);

   for(int i(0); (i<15) && !hasConverged ; ++i)
   {
      initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
            nKemu, nSignal, nPartReco, nComb, fracZero, fracOne,
            nJpsiLeak, constPartReco, fracPartRecoSigma, 
            *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, l1Kemu, l2Kemu, l3Kemu, l4Kemu, l5Kemu, 
            *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen);

      cout<<"FITTING: starting with nsignal = "<<nSignal.getValV()<<" refit nbr. "<<i<<endl;
      //if(fitRes != NULL && fitRes != 0) delete fitRes;

      migradRes = minuit.migrad();
      hesseRes = minuit.hesse();

      fitRes = minuit.save();
      edm = fitRes->edm();

      fitResVec.push_back(fitRes); 
      migradResVec.push_back(migradRes);
      hesseResVec.push_back(hesseRes);

      if( migradRes == 0 && hesseRes == 0 && edm < 1e-3 ) hasConverged = true;

      ++nrefit;

      cout<<"Fitting nbr "<<i<<" done. Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   if(!hasConverged)
   {
      double minNll(1e20);
      int minIndex(-1);
      for(unsigned int i(0); i<fitResVec.size(); ++i)
      {
         if( fitResVec.at(i)->minNll() < minNll)
         {
            minIndex = i;
            minNll = fitResVec[i]->minNll();
         }
      }
      
      migradRes = migradResVec.at(minIndex);
      hesseRes = hesseResVec.at(minIndex);
      cout<<"Fit not converged, choose fit "<<minIndex<<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(t, fitRes,  update, migradRes, hesseRes, hasConverged);

   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);
   //totPdf.fitTo(*dataGenTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(wantplot) plot_fit_result(plotsfile, totPdf, *dataGenTot);
   if(wantplot) plot_kemu_fit_result(plotsfile, totKemuPdf, *dataGenKemu);

   fw.Close();
   //delete and return
   delete nll;
   delete workspace;
   delete workspaceGen;
   delete combPDF;
   delete KemuPDF;
}


void FitterUtilsSimultaneousExpOfPolyTimesX::plot_kemu_fit_result(string plotsfile, RooAbsPdf &totKemuPdf, RooDataSet const& dataGenKemu)
{

   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totKemuPdf.getObservables(dataGenKemu);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   while((var = (RooRealVar*) iter->Next()))
   {
      frame = var->frame();
      dataGenKemu.plotOn(frame);
      totKemuPdf.plotOn(frame, LineColor(kRed));

      plots.push_back(frame);
   }

   if (!(plots.size())) return;

   TCanvas cFit("cKemuFit", "cKemuFit", 600, 800);
   cFit.Divide(1,2);
   cFit.cd(1); plots[0]->Draw();
   if (plots.size()>1){ 
      cFit.cd(2); plots[1]->Draw();
   }

   cFit.Write();
   f2.Close();
}

