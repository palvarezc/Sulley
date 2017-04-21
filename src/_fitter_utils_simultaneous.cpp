#include "fitter_utils_simultaneous.h"
#include "RooBinning.h"
#include "RooRandom.h"


FitterUtilsSimultaneous::FitterUtilsSimultaneous(string workspacename_, const Options opts_)
  :FitterUtils(workspacename_,opts_)
   {

     Kemu = Combinatorial(workspacename+"Kemu",opts.fittype);
   }




void FitterUtilsSimultaneous::prepare_PDFs()
{

  FitterUtils::prepare_PDFs();
  
   TFile fw(workspacename.c_str(), "UPDATE");
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");

   RooCategory category("category", "category");
   category.defineType("Kee");
   category.defineType("Kemu");

   workspace->import(category);
   workspace->Write("", TObject::kOverwrite);
   fw.Close();
   
}

void FitterUtilsSimultaneous::generate()
{
   FitterUtils::generate();
   TFile fw(workspacename.c_str(), "UPDATE");
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");

   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooArgSet *parsettrue = (RooArgSet*) workspace->set((combinatorial.get_name()+"_parsettrue").c_str());
   RooArgSet *combparset = (RooArgSet*) workspace->set((combinatorial.get_name()+"_parset").c_str());

   RooArgSet argset(*B_plus_M);
   if (opts.fit2D) argset.add(*misPT);


   RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
   Kemu.generate("Kemu", &argset, parsettrue, opts.nGenKemu, workspaceGen,m_wantplot);
   workspaceGen->Write("", TObject::kOverwrite);


   RooArgSet Kemuset((Kemu.get_name()+"_parset").c_str());
   TIterator *iter = combparset->createIterator();
   RooRealVar *combpar;
   RooRealVar *expoConstKemu;   

   while((combpar = (RooRealVar*) iter->Next()))
   {

     combpar->setConstant(false);

     if (combpar->GetName()=="expoConst")
     {
       expoConstKemu = new RooRealVar(*combpar);
       expoConstKemu->SetName("expoConstKemu");
       expoConstKemu->SetTitle("expoConstKemu");
       Kemuset.add(*expoConstKemu);
       continue;
     }
     
     Kemuset.add(*combpar);

   }

   workspace->import(Kemuset);
   workspace->Write("", TObject::kOverwrite);
   fw.Close();

   delete expoConstKemu;
   
}





void FitterUtilsSimultaneous::run_toy(bool constPartReco,
                          double fracPartReco_const,
                          ofstream& out, TTree* t, bool update)
{

  //***************Generate
  generate();

  TFile fw(workspacename.c_str());

  //***************Retrieve generated samples
  RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
  RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGamma");
  RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen->data("dataGenSignalOneGamma");
  RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGamma");
  RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen->data("dataGenPartReco");
  RooDataSet* dataGenComb = (RooDataSet*)workspaceGen->data("dataGenComb");
  RooDataSet* dataGenKemu = (RooDataSet*)workspaceGen->data("dataGenKemu");
  RooDataSet* dataGenJpsiLeak(0);
  if(opts.nGenJpsiLeak>0) dataGenJpsiLeak = (RooDataSet*)workspaceGen->data("dataGenJpsiLeak");
  


  //***************Merge datasets
  RooDataSet* dataGenTot(dataGenPartReco);
  dataGenTot->append(*dataGenSignalZeroGamma);
  dataGenTot->append(*dataGenSignalOneGamma);
  dataGenTot->append(*dataGenSignalTwoGamma);
  dataGenTot->append(*dataGenComb);
  if(opts.nGenJpsiLeak>0) dataGenTot->append(*dataGenJpsiLeak);


  //***************Fit
  RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
  fit(dataGenTot, dataGenKemu, workspace, constPartReco, fracPartReco_const, out, t, update);

  fw.Close();


}


void FitterUtilsSimultaneous::fit(RooDataSet* datasetTot, RooDataSet* datasetKemu, 
                                  RooWorkspace* workspace,
                                  bool constPartReco,
                                  double fracPartReco_const,
                                  ofstream& out, TTree* t, bool update)
{

  //***************Get the PDFs from the workspace
  RooRealVar *B_plus_M = workspace->var("B_plus_M");
  RooRealVar *misPT = workspace->var("misPT");
  RooCategory *category = workspace->cat("category");
  RooArgSet *combparset = (RooArgSet*) workspace->set((combinatorial.get_name()+"_parset").c_str());
  RooArgSet *Kemuparset = (RooArgSet*) workspace->set((Kemu.get_name()+"_parset").c_str());
  RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");
  RooArgSet observables(*B_plus_M);
  if (opts.fit2D) observables.add(*misPT);


  //***************Simultaneous dataset
  RooDataSet dataGenSimultaneous("dataGenSimultaneous", "dataGenSimultaneous", observables, Index(*category), 
                                 Import("Kee", *datasetTot), Import("Kemu", *datasetKemu));

  //**************Prepare fitting function    
  RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
  RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
  RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
  RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
  RooHistPdf *histPdfJpsiLeak(0);
  if(opts.nGenJpsiLeak>0) histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");

  RooAbsPdf* combPDF =  combinatorial.build_model(&observables, 0, combparset);
  RooAbsPdf* KemuPDF =  Kemu.build_model(&observables, 0, Kemuparset);

  // expoConst->setVal(trueExp->getVal());
 

  //**************Prepare fitting function
  
  RooRealVar nSignal("nSignal", "#signal events", 1.*opts.nGenSignal, opts.nGenSignal-7*sqrt(opts.nGenSignal), 
                     opts.nGenSignal+7*sqrt(opts.nGenSignal));
  RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*opts.nGenPartReco, opts.nGenPartReco-7*sqrt(opts.nGenPartReco), 
                       opts.nGenPartReco+7*sqrt(opts.nGenPartReco));
  RooRealVar nComb("nComb", "#nComb", 1.*opts.nGenComb, opts.nGenComb-7*sqrt(opts.nGenComb), opts.nGenComb+7*sqrt(opts.nGenComb));
  RooRealVar nKemu("nKemu", "#nKemu", 1.*opts.nGenKemu, opts.nGenKemu-7*sqrt(opts.nGenKemu), opts.nGenKemu+7*sqrt(opts.nGenKemu));
  RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*opts.nGenJpsiLeak, opts.nGenJpsiLeak-7*sqrt(opts.nGenJpsiLeak), 
                       opts.nGenJpsiLeak+7*sqrt(opts.nGenJpsiLeak));
  RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
  RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);
  RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(nPartReco,nSignal));
  RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));

  RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", 
                          RooArgList(*histPdfSignalZeroGamma, *histPdfSignalOneGamma, *histPdfSignalTwoGamma), 
                          RooArgList(fracZero, fracOneRec));

  RooArgList pdfList(histPdfSignal, *histPdfPartReco, *combPDF);
  RooArgList yieldList(nSignal, nPartReco, nComb);
  
  if(opts.nGenJpsiLeak>0)
  {
    pdfList.add(*histPdfJpsiLeak);
    yieldList.add(nJpsiLeak); 
  }
  RooAddPdf totPdf("totPdf", "totPdf", pdfList, yieldList);
  RooExtendPdf totKemuPdf("totKemuPdf", "totKemuPdf", *KemuPDF, nKemu);
  
  //**************** Prepare simultaneous PDF
  
  RooSimultaneous simPdf("simPdf", "simPdf", *category);
  simPdf.addPdf(totPdf, "Kee");
  simPdf.addPdf(totKemuPdf, "Kemu");

  //**************** Constrain the fraction of zero and one photon

  int nGenSignalZeroGamma(floor(opts.nGenFracZeroGamma*opts.nGenSignal));
  int nGenSignalOneGamma(floor(opts.nGenFracOneGamma*opts.nGenSignal));
  int nGenSignalTwoGamma(floor(opts.nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));
  
  RooRealVar fracZeroConstMean("fracZeroConstMean", "fracZeroConstMean", nGenSignalZeroGamma*1./opts.nGenSignal);
  RooRealVar fracZeroConstSigma("fracZeroConstSigma", "fracZeroConstSigma", sqrt(nGenSignalZeroGamma)/opts.nGenSignal);
  RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroConstMean, fracZeroConstSigma); 
  
  RooRealVar fracOneConstMean("fracOneConstMean", "fracOneConstMean", 
                              nGenSignalOneGamma*1./opts.nGenSignal/(1-fracZeroConstMean.getVal()));
  RooRealVar fracOneConstSigma("fracOneConstSigma", "fracOneConstSigma", 
                               sqrt(nGenSignalOneGamma)/opts.nGenSignal/(1-fracZeroConstMean.getVal()));
  RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneConstMean, fracOneConstSigma); 
  
  RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", opts.nGenPartReco/(1.*opts.nGenSignal));
  RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());
  
  RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);
  
  RooRealVar JpsiLeakMean("JpsiLeakMean", "JpsiLeakMean", opts.nGenJpsiLeak);
  RooRealVar JpsiLeakSigma("JpsiLeakSigma", "JpsiLeakSigma", opts.nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
  RooGaussian JpsiLeakConst("JpsiLeakConst", "JpsiLeakConst", nJpsiLeak, JpsiLeakMean, JpsiLeakSigma); 
  

  //**************** fit
  
  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;

  RooArgSet *par_set = simPdf.getParameters(dataGenSimultaneous);
  initiateParams(par_set);
  
  // initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, *trueExp,
  //                *trueT, *trueN,
  //                nKemu, nSignal, nPartReco,
  //                nComb, fracZero, fracOne,
  //                *expoConst, expoConstKemu, *T, *n,
  //                nJpsiLeak, constPartReco, fracPartRecoSigma);

   RooArgSet constraints(fracZeroConst, fracOneConst);
   if (constPartReco) constraints.add(fracPartRecoConst);
   if(opts.nGenJpsiLeak>0) constraints.add(JpsiLeakConst);

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

   for(int i(0); (i<10) && !hasConverged ; ++i)
   {
      // initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, *trueExp,
      //       *trueT, *trueN,
      //       nKemu, nSignal, nPartReco,
      //       nComb, fracZero, fracOne,
      //       *expoConst, expoConstKemu, *T, *n,
      //       nJpsiLeak, constPartReco, fracPartRecoSigma);
     initiateParams(par_set);
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

   if(m_wantplot) plot_fit_result(totPdf, *datasetTot);
   if(m_wantplot) plot_kemu_fit_result(totKemuPdf, *datasetKemu);

   //delete and return
   delete nll;
   delete combPDF;
   delete KemuPDF;
}


void FitterUtilsSimultaneous::plot_kemu_fit_result(RooAbsPdf &totKemuPdf, RooDataSet const& datasetKemu)
{

   //**************Prepare TFile to save the plots

   TFile f2(opts.plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totKemuPdf.getObservables(datasetKemu);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   while((var = (RooRealVar*) iter->Next()))
   {
      frame = var->frame();
      datasetKemu.plotOn(frame);
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

void FitterUtilsSimultaneous::initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, 
                                             RooRealVar const& expoConstGen,
                                             RooRealVar const& TGen, RooRealVar const& nGen,
                                             RooRealVar& nKemu, RooRealVar& nSignal, RooRealVar& nPartReco,
                                             RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne,
                                             RooRealVar& expoConst, RooRealVar& expoConstKemu,
                                             RooRealVar& T, RooRealVar& n,
                                             RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma)
{
  FitterUtils::initiateParams(expoConstGen,
                              nSignal, nPartReco, nComb, fracZero, fracOne, 
                              expoConst, nJpsiLeak, constPartReco, fracPartRecoSigma); 
  
  TRandom rand;
  rand.SetSeed();

  nKemu.setVal(rand.Uniform(opts.nGenKemu-5*sqrt(opts.nGenKemu), opts.nGenKemu+5*sqrt(opts.nGenKemu)));
  nKemu.setRange(opts.nGenKemu-10*sqrt(opts.nGenKemu), opts.nGenKemu+10*sqrt(opts.nGenKemu));

  T.setVal(rand.Uniform( TMath::Max(0.,TGen.getVal() - 5*TGen.getError()), TGen.getVal() + 5*TGen.getError()));
  T.setRange(TMath::Max(TGen.getVal() - 10*TGen.getError(),0.), TGen.getVal() + 10*TGen.getError());

  n.setVal(rand.Uniform(TMath::Max( nGen.getVal() - 5*nGen.getError(),2.), nGen.getVal() + 5*nGen.getError()));
  n.setRange(TMath::Max(nGen.getVal() - 10*nGen.getError(),2.), nGen.getVal() + 10*nGen.getError());

  expoConstKemu.setVal(rand.Uniform( expoConstGen.getVal() - 5*expoConstGen.getError(), 
                                     expoConstGen.getVal() + 5*expoConstGen.getError() ) );
  expoConstKemu.setRange( expoConstGen.getVal() - 10*expoConstGen.getError(), 
                          expoConstGen.getVal() + 10*expoConstGen.getError() );
}
