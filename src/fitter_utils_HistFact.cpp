#include "fitter_utils_HistFact.h"
#include "RooBinning.h"
#include "RooRandom.h"
#include "RooExpOfPolyTimesXBinned.h"
#include "RooExpOfPolyTimesX.h"
#include "RooExpOfPoly.h"

#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooRealSumPdf.h"
#include "RooProduct.h"
#include "TIterator.h"

string measurementName = "my_measurement";
string channelName = "B2Kee";


FitterUtilsHistFact::FitterUtilsHistFact(string workspacename_, const Options opts_)
  :FitterUtils(workspacename_,opts_)   
{}


void FitterUtilsHistFact::prepare_PDFs(int fitmode)
{

  bool templateStat = opts.templatestat;
  
  FitterUtils::prepare_PDFs(fitmode);

  TFile fw(workspacename.c_str(),"UPDATE");
  RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
  TH2D* h_SignalZero = (TH2D*) workspace->obj("SignalZeroGamma");
  TH2D* h_SignalOne = (TH2D*) workspace->obj("SignalOneGamma");
  TH2D* h_SignalTwo = (TH2D*) workspace->obj("SignalTwoGamma");
  TH2D* h_PartReco = (TH2D*) workspace->obj("PartReco");
  TH2D* h_JpsiLeak = (TH2D*) workspace->obj("JpsiLeak");
  

  //**************Create PDFs
  cout<<"Creating pdfs..."<<endl;

  double fzero = 0.368;
  double fone = 0.484;
  double ftwo = 1- fzero - fone;
  
  RooStats::HistFactory::Measurement meas(measurementName.c_str(),measurementName.c_str());
  string measPath = workspacename;
  measPath.resize(measPath.find("/"));
  meas.SetOutputFilePrefix(measPath+"/");
  meas.SetExportOnly(kTRUE);
  
  meas.SetPOI("nSignal");

  

  meas.SetLumi(1.0);
  meas.SetLumiRelErr(0.05);  
  
  RooStats::HistFactory::Channel chan(channelName.c_str());
  chan.SetStatErrorConfig(1e-5,"Poisson");

  //Setup the zero brem 
  RooStats::HistFactory::Sample sigzero("SignalZeroSample");
  sigzero.SetHisto(h_SignalZero);
  if(templateStat) sigzero.ActivateStatError();
  sigzero.SetNormalizeByTheory(kFALSE);
  sigzero.AddNormFactor("nSignal", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
  sigzero.AddNormFactor("fracZero", fzero,0.,1.);
  sigzero.AddNormFactor("mcNorm_sigzero", 1./h_SignalZero->Integral(), 1e-9, 1.);
  // sigzero.AddOverallSys("fracZeroConstraint",0.9,1.1); // Relative to fracZero
  chan.AddSample(sigzero);


  //Setup the one brem 
  RooStats::HistFactory::Sample sigone("SignalOneSample");
  if(templateStat) sigone.ActivateStatError();
  sigone.SetHisto(h_SignalOne);
  sigone.SetNormalizeByTheory(kFALSE);
  sigone.AddNormFactor("nSignal",  1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
  sigone.AddNormFactor("fracOneRec", fone,0.,1.);
  sigone.AddNormFactor("mcNorm_sigone", 1./h_SignalOne->Integral(), 1e-9, 1.);
  // sigone.AddOverallSys("fracOneConstraint",0.9,1.1); // Relative to fracOne
  chan.AddSample(sigone);



  //Setup the two brem 
  RooStats::HistFactory::Sample sigtwo("SignalTwoSample");
  if(templateStat) sigtwo.ActivateStatError();
  sigtwo.SetHisto(h_SignalTwo);
  sigtwo.SetNormalizeByTheory(kFALSE);
  sigtwo.AddNormFactor("nSignal",  1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
  sigtwo.AddNormFactor("fracTwo", ftwo,0.,1.);
  sigtwo.AddNormFactor("mcNorm_sigtwo", 1./h_SignalTwo->Integral(), 1e-9, 1.);
  // sigtwo.AddOverallSys("fracTwoConstraint",0.9,1.1); // Relative to fracTwo (need to figure out how to lock it to 1-zero-one)
  chan.AddSample(sigtwo);
  

  //Setup the part reco
  RooStats::HistFactory::Sample partreco("PartRecoSample");
  if(templateStat) partreco.ActivateStatError();
  partreco.SetHisto(h_PartReco);
  partreco.SetNormalizeByTheory(kFALSE);
  partreco.AddNormFactor("nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
  partreco.AddNormFactor("mcNorm_partreco", 1./h_PartReco->Integral(), 1e-9, 1.);
  // if (opts.constPartReco) partreco.AddOverallSys("PartRecoConstraint",0.9,1.1); // Relative to nPartReco
  chan.AddSample(partreco);
  

  //Setup the jpsi leak
    
  RooStats::HistFactory::Sample jpsileak("JpsiLeakSample");
  // if(templateStat) jpsileak.ActivateStatError();
  jpsileak.SetHisto(h_JpsiLeak);
  jpsileak.SetNormalizeByTheory(kFALSE);
  jpsileak.AddNormFactor("nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
  jpsileak.AddNormFactor("mcNorm_jpsileak", 1./h_JpsiLeak->Integral(), 1e-9, 1.);
  // if (constJpsiLeak) jpsileak.AddOverallSys("JpsiLeakConstraint",0.9,1.1); // Relative to nJpsiLeak
  if (nGenJpsiLeak>1)    chan.AddSample(jpsileak); 
  
  

  meas.AddChannel(chan);
  meas.AddPreprocessFunction("fracOneRec","(1-fracZero)*fracOne","fracOne[0,1],fracZero[0,1]");
  meas.AddPreprocessFunction("fracTwo","1-fracZero-(1-fracZero)*fracOne","fracOne[0,1],fracZero[0,1]");

  // meas.CollectHistograms();


  //Create the workspace that is going to save the config
  RooWorkspace *workspace_meas;
  workspace_meas = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);

  // workspace_meas = new RooWorkspace_Meas("combined");
  
  
  // ModelConfig *mc = (ModelConfig*) workspace_meas->obj("ModelConfig"); // Get model (noComb) manually
   
  // // Lets tell roofit the right names for our histogram variables
  // RooArgSet *obs = (RooArgSet*) mc->GetObservables();
  // RooRealVar *B_plus_M_mod = (RooRealVar*) obs->find(("obs_x_"+channelName).c_str());
  // RooRealVar *misPT_mod = (RooRealVar*) obs->find(("obs_y_"+channelName).c_str());
  // // B_plus_M_mod->SetName("B_plus_M");
  // B_plus_M_mod->SetTitle("B_plus_M");
  // B_plus_M_mod->setUnit("MeV/c^{2}");
  // misPT_mod->SetTitle("misPT");
  // // misPT_mod->SetName("misPT");
  // misPT_mod->setUnit("MeV/c");



  // //***************Create combinatorial from fit to data
  // // fit_combinatorial(combfile, combtree, CombCuts, argsetComb,
  // //                   *B_plus_M_mod, *misPT_mod, workspace_meas);


  // fit_combinatorial(combfile, combtree, CombCuts, argsetComb,
  //                   B_plus_M, misPT, workspace_meas);
  
  RooRealVar nComb("nComb", "#nComb", 1.*opts.nGenComb, opts.nGenComb-7*sqrt(opts.nGenComb), 
                   opts.nGenComb+7*sqrt(opts.nGenComb));
  RooRealVar nComb_fit("nComb_fit", "#nComb_fit", 1.*opts.nGenComb, opts.nGenComb-7*sqrt(opts.nGenComb), 
                       opts.nGenComb+7*sqrt(opts.nGenComb));
  workspace_meas->import(nComb);  
  workspace_meas->import(nComb_fit);

  //***************Save some extra factors
  // workspace_meas->factory(Form("nMC_SignalZero[%f]",h_SignalZero->Integral()));
  // workspace_meas->factory(Form("nMC_SignalOne[%f]",h_SignalOne->Integral()));
  // workspace_meas->factory(Form("nMC_SignalTwo[%f]",h_SignalTwo->Integral()));
  // workspace_meas->factory(Form("nMC_PartReco[%f]",h_PartReco->Integral()));
  // if (nGenJpsiLeak>1) workspace_meas->factory(Form("nMC_JpsiLeak[%f]",h_JpsiLeak->Integral()));
  // workspace_meas->factory(Form("nMC_Comb[%f]",h_Comb->Integral()));
  // workspace_meas->factory(Form("nMC_CombExt[%f]",h_CombExt->Integral()));


  workspace_meas->Print();
  workspace->merge(*workspace_meas);
  workspace->Write("", TObject::kOverwrite);
   
  cout<<"PDFs prepared!"<<endl;
  delete workspace_meas;  
  fw.Close();

}


void FitterUtilsHistFact::run_toy(double fracPartReco_const,
                          ofstream& out, TTree* t, bool update)
{

  //***************Generate
  FitterUtils::generate(0);

  TFile fw(workspacename.c_str());

  //***************Retrieve generated samples
  RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");
  RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGamma");
  RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen->data("dataGenSignalOneGamma");
  RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGamma");
  RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen->data("dataGenPartReco");
  RooDataSet* dataGenComb = (RooDataSet*)workspaceGen->data("dataGenComb");
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
  fit(dataGenTot, workspace, fracPartReco_const, out, t, update);

  fw.Close();


}

void FitterUtilsHistFact::fit(RooDataSet* dataGenTot_chan, RooWorkspace* workspace,
                      double fracPartReco_const,
                      ofstream& out, TTree* t, bool update)

{


  bool constFracs = 0;  
  bool constComb = 0;  

  // Get observables
  ModelConfig *mc = (ModelConfig*) workspace->obj("ModelConfig");
  RooArgSet *obs = (RooArgSet*) mc->GetObservables();
  RooRealVar *B_plus_M = (RooRealVar*) obs->find(("obs_x_"+channelName).c_str());
  RooRealVar *misPT = (RooRealVar*) obs->find(("obs_y_"+channelName).c_str());
  B_plus_M->SetTitle("B_plus_M");
  B_plus_M->setUnit("MeV/c^{2}");
  misPT->SetTitle("misPT");
  misPT->setUnit("MeV/c");
  cout<<"B_plus_M nbins: "<<B_plus_M->getBins("default")<<endl;
  cout<<"misPT nbins: "<<misPT->getBins("default")<<endl;
  RooArgSet *combparset = (RooArgSet*) workspace->set((combinatorial.get_name()+"_parset").c_str());


  RooArgSet observables(*B_plus_M);
  if (opts.fit2D) observables.add(*misPT);

  

  RooRealVar* nSignal = workspace->var("nSignal");
  RooRealVar* nPartReco = workspace->var("nPartReco");
  RooRealVar* nComb = workspace->var("nComb");
  RooRealVar* nComb_fit = workspace->var("nComb_fit");
  RooRealVar* fracZero = workspace->var("fracZero");
  RooRealVar* fracOne = workspace->var("fracOne");
  RooFormulaVar* fracTwo = (RooFormulaVar*) workspace->obj("fracTwo");
  RooRealVar* nJpsiLeak;
  if (opts.nGenJpsiLeak>0)  nJpsiLeak = workspace->var("nJpsiLeak");
  else nJpsiLeak = new RooRealVar("nJpsiLeak", "#nJpsiLeak", 1.*opts.nGenJpsiLeak, opts.nGenJpsiLeak-7*sqrt(opts.nGenJpsiLeak), opts.nGenJpsiLeak+7*sqrt(opts.nGenJpsiLeak));
  
  RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");
  RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(*nPartReco,*nSignal));


  //***************** Build model
   
  // Get HistFactory pdf (no combinatorial)
  RooSimultaneous *model_nocomb = (RooSimultaneous*)mc->GetPdf();
  RooArgList obsTerms;
  RooArgList pois_constraints;
  FactorizeHistFactoryPdf(*obs, *model_nocomb, obsTerms, pois_constraints);

  RooRealSumPdf* model_nocomb_pdf = ( RooRealSumPdf*) obsTerms.first();
   
  RooArgList *model_nocomb_pdf_coefList = (RooArgList*) model_nocomb_pdf->coefList().Clone();
  RooArgList *model_nocomb_pdf_funcList = (RooArgList*) model_nocomb_pdf->funcList().Clone();


  // Define combinatorial model
  RooAbsPdf *combPDF_unnorm = combinatorial.build_model(&observables, combparset);
  RooAbsReal *combIntegral = combPDF_unnorm->createIntegral(RooArgSet(*B_plus_M, *misPT));
  combIntegral->SetNameTitle("combIntegral","combIntegral");
  // RooFormulaVar *combNorm = new RooFormulaVar("combNorm","1./combIntegral",RooArgList(*combIntegral));
  RooRealVar *combNorm = new RooRealVar("combNorm","combNorm",1./combIntegral->getVal());
  // RooRealVar *combNorm = new RooRealVar("combNorm","combNorm",1.);
  RooProduct *combPDF = new RooProduct("combPDF","combPDF",RooArgList(*combPDF_unnorm,*combNorm));  
   

  model_nocomb_pdf_coefList->add(*nComb);
  model_nocomb_pdf_funcList->add(*combPDF); //normalization

   // Combined model
   RooRealSumPdf *model_pdf = new RooRealSumPdf((channelName+"_comb").c_str(),(channelName+"_comb").c_str(),
                                                *model_nocomb_pdf_funcList,*model_nocomb_pdf_coefList,1);  
   model_pdf->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator") ;
   model_pdf->specialIntegratorConfig(kTRUE)->method2D().setLabel("RooBinIntegrator") ;
   model_pdf->specialIntegratorConfig(kTRUE)->methodND().setLabel("RooBinIntegrator") ;
   model_pdf->forceNumInt();
   model_pdf->setAttribute("GenerateBinned"); // for use with RooSimultaneous::generate in mixed mode

   // Add constraints
   RooArgList *model_and_constraints = (RooArgList*)pois_constraints.Clone();
   model_and_constraints->add(*model_pdf);
   RooProdPdf *model_pdf_constrain = new RooProdPdf("RK_kinematic_comb_constrain","RK_kinematic_comb_constrain",*model_and_constraints);

   // Total model
   RooCategory *idx = (RooCategory*) obs->find("channelCat");
   RooSimultaneous model("modified_pdf","modified_pdf",*idx);
   model.addPdf(*model_pdf_constrain,idx->getLabel());  
   // model.addPdf(*model_nocomb_pdf,idx->getLabel());  

   HistFactorySimultaneous totPdf( model );
   // HistFactorySimultaneous totPdf( *model_nocomb );
   totPdf.SetNameTitle("totPdf", "totPdf");
   

   //***************** Get data


   cout<<"CACA3"<<endl;


   // if(m_wantplot)
   // {
   //    //**************Must get the datasets

   //    RooDataSet* dataSetSignalZeroGamma = (RooDataSet*)workspace->data("dataSetSignalZeroGamma");
   //    RooDataSet* dataSetSignalOneGamma = (RooDataSet*)workspace->data("dataSetSignalOneGamma");
   //    RooDataSet* dataSetSignalTwoGamma = (RooDataSet*)workspace->data("dataSetSignalTwoGamma");
   //    RooDataSet* dataSetPartReco = (RooDataSet*)workspace->data("dataSetPartReco");
   //    RooDataSet* dataSetComb = (RooDataSet*)workspace->data("dataSetComb");
   //    RooDataSet* dataSetJpsiLeak = (RooDataSet*)workspace->data("dataSetJpsiLeak");

   //    //**************Plot all the different components
   // cout<<"CACA4"<<endl;

   //    cout<<"dataGenSignalZeroGamma: "<<dataGenSignalZeroGamma<<endl;
   //    PlotShape(*dataSetSignalZeroGamma, *dataGenSignalZeroGamma, *histPdfSignalZeroGamma, plotsfile, "cSignalZeroGamma", *B_plus_M, *misPT);
   //    PlotShape(*dataSetSignalOneGamma, *dataGenSignalOneGamma, *histPdfSignalOneGamma, plotsfile, "cSignalOneGamma", *B_plus_M, *misPT);
   //    PlotShape(*dataSetSignalTwoGamma, *dataGenSignalTwoGamma, *histPdfSignalTwoGamma, plotsfile, "cSignalTwoGamma", *B_plus_M, *misPT);
   //    PlotShape(*dataSetPartReco, *dataGenPartReco, *histPdfPartReco, plotsfile, "cPartReco", *B_plus_M, *misPT);
   //    PlotShape(*dataSetComb, *dataGenComb, *combPDF, plotsfile, "cComb", *B_plus_M, *misPT);
   //    PlotShape(*dataSetComb, *dataGenComb, *histPdfComb, plotsfile, "cCombTemplate", *B_plus_M, *misPT);

   // cout<<"CACA5"<<endl;
   //    if(opts.nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeak, *dataGenJpsiLeak, *histPdfJpsiLeak, plotsfile, "cJpsiLeak", *B_plus_M, *misPT);
   // }

   //***************Merge datasets

   RooDataHist *dataGenTot_chan_bin = dataGenTot_chan->binnedClone();
   RooDataHist *dataGenTot = new RooDataHist("dataGenTot","dataGenTot",RooArgSet(*B_plus_M, *misPT),
                                             Index(*idx),Import(idx->getLabel(),*dataGenTot_chan_bin));

   cout<<"Binned data number of bins: "<<dataGenTot->numEntries()<<endl;
   

   //**************Setup constant parameters
  // ((RooRealVar*)(mc->GetNuisanceParameters()->find("fracZero")))->setConstant(kTRUE);
  // ((RooRealVar*)(mc->GetNuisanceParameters()->find("fracOne")))->setConstant(kTRUE);
  // ((RooRealVar*)(mc->GetNuisanceParameters()->find("fracTwo")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_sigzero")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_sigone")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_sigtwo")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_partreco")))->setConstant(kTRUE);
  if(opts.nGenJpsiLeak>1) ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_jpsileak")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);
  // ((RooRealVar*)(mc->GetNuisanceParameters()->find("l1Kee")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("l2Kee")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("l3Kee")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("l4Kee")))->setConstant(kTRUE);
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("l5Kee")))->setConstant(kTRUE);
  
  // nJpsiLeak->setVal(opts.nGenJpsiLeak);
  // nJpsiLeak->setConstant(kTRUE);
  
   cout<<"CACA7"<<endl;
  

   //**************** Constrain the fraction of zero and one photon


   int nGenSignalZeroGamma(floor(opts.nGenFracZeroGamma*opts.nGenSignal));
   int nGenSignalOneGamma(floor(opts.nGenFracOneGamma*opts.nGenSignal));
   int nGenSignalTwoGamma(floor(opts.nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));

   RooRealVar fracZeroMean("fracZeroMean", "fracZeroMean", nGenSignalZeroGamma*1./opts.nGenSignal);
   RooRealVar fracZeroSigma("fracZeroSigma", "fracZeroSigma", sqrt(nGenSignalZeroGamma)/opts.nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", *fracZero, fracZeroMean, fracZeroSigma); 

   RooRealVar fracOneMean("fracOneMean", "fracOneMean", 
                               nGenSignalOneGamma*1./opts.nGenSignal/(1-fracZeroMean.getVal()));
   RooRealVar fracOneSigma("fracOneSigma", "fracOneSigma", 
                                sqrt(nGenSignalOneGamma)/opts.nGenSignal/(1-fracZeroMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", *fracOne, fracOneMean, fracOneSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", opts.nGenPartReco/(1.*opts.nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());
   RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

   RooRealVar nJpsiLeakMean("nJpsiLeakMean", "nJpsiLeakMean", opts.nGenJpsiLeak);
   RooRealVar nJpsiLeakSigma("nJpsiLeakSigma", "nJpsiLeakSigma", opts.nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
   RooGaussian nJpsiLeakConst("nJpsiLeakConst", "nJpsiLeakConst", *nJpsiLeak, nJpsiLeakMean, nJpsiLeakSigma); 
   

   fracZero->setVal(fracZeroMean.getVal());
   fracOne->setVal(fracOneMean.getVal());
   
   
   cout<<"Evaluating fracTwo... "<<endl;
   cout<<"    fracZero = "<<fracZero->getVal()<<endl;
   cout<<"    fracOne = "<<fracOne->getVal()<<endl;
   cout<<"    fracTwo = "<<fracTwo->getVal()<<endl;
   
   cout<<"CACA10"<<endl;
   //Extra TEST CONSTRAINT


   //RooRealVar combConstMean("combConstMean", "combConstMean", opts.nGenComb);
   //RooRealVar combConstSigma("combConstSigma", "combConstSigma", 7.7);
   //RooGaussian combConst("combConst", "combConst", nComb, combConstMean, combConstSigma);

   //**************** fit
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;

   RooArgSet *par_set = totPdf.getParameters(*dataGenTot);
   RooArgSet par_set_const;

   RooArgSet constraints(fracZeroConst, fracOneConst);
   par_set_const.add(fracZeroMean);
   par_set_const.add(fracZeroSigma);
   par_set_const.add(fracOneMean);
   par_set_const.add(fracOneSigma);
   if (opts.constPartReco) 
     {
       constraints.add(fracPartRecoConst);
       par_set_const.add(fracPartRecoMean);
       par_set_const.add(fracPartRecoSigma);
     }
   if(opts.nGenJpsiLeak>0) 
     {
       constraints.add(nJpsiLeakConst);
       par_set_const.add(nJpsiLeakMean);
       par_set_const.add(nJpsiLeakSigma);
     }

   initiateParams(*par_set, constraints, par_set_const);

   RooArgList* gammas = new RooArgList();
   ParamHistFunc* param_func=NULL;
   // bool hasStatUncert = getStatUncertaintyFromChannel(model_nocomb_pdf , param_func, gammas );
   cout<<"CACA11b"<<endl;
   
   // bool hasStatUncert = getStatUncertaintyFromChannel(model_pdf_constrain , param_func, gammas );

   // cout<< "Param_func nbins: "<<param_func->numBins()<<endl;
   

   RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Verbose(kTRUE), ExternalConstraints(constraints));
   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, ExternalConstraints(constraints),Offset(kTRUE));
   // RooAbsReal* nll = totPdf.createNLL(*dataGenTot, Extended(), ExternalConstraints(constraints));
   RooMinuit minuit(*nll);
   minuit.setStrategy(2);
   minuit.setErrorLevel(0.5);
   minuit.setPrintLevel(1);   

   // std::cout << "Minimizing the Minuit (Migrad)" << std::endl;
   // w->saveSnapshot("TMCPARS",*allpars,kTRUE);
   // minuit_hf->setStrategy(2);
   // minuit_hf->fit("smh");
   // RooFitResult *tempResult=minuit_hf->save("TempResult","TempResult");
   // cout << tempResult->edm() << endl;
   // if (useMinos) minuit_hf->minos(RooArgSet(*poi));
   // result = minuit_hf->save("Result","Result");
 
   int migradRes(1);
   int hesseRes(4);

   vector<int> migradResVec;
   vector<int> hesseResVec;

   double edm(10);
   int nrefit(0);

   RooFitResult* fitRes(0);
   vector<RooFitResult*> fitResVec;

   bool hasConverged(false);

   // double int_comb = 1;
   

   for(int i(0); (i<15) && !hasConverged ; ++i)
   {

   initiateParams(*par_set, constraints, par_set_const);

      // initiateParams(nGenSignalZeroGamma, nGenSignalOneGamma, nGenSignalTwoGamma, 
      //     *nSignal, *nPartReco, *nComb, *fracZero, *fracOne, *nJpsiLeak,  opts.constPartReco, fracPartRecoSigma,
      //                *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee, *l1KeeGen, *l2KeeGen, *l3KeeGen, *l4KeeGen, *l5KeeGen, constFracs, constComb);

      cout<<"FITTING: starting with nsignal = "<<nSignal->getValV()<<" refit nbr. "<<i<<endl;
      //if(fitRes != NULL && fitRes != 0) delete fitRes;

      // cout<<"Eval nll: "<<nll->getVal()<<endl;
      
      

      migradRes = minuit.migrad();
      hesseRes = minuit.hesse();

      fitRes = minuit.save();
      edm = fitRes->edm();


      // int_comb = combIntegral->getVal();
      // nComb->setVal(nComb_fit->getVal()*int_comb);
      // nComb->setError(nComb_fit->getError()*int_comb);
      
      

      fitResVec.push_back(fitRes); 
      migradResVec.push_back(migradRes);
      hesseResVec.push_back(hesseRes);

      if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;

      ++nrefit;


      cout<<"Fitting nbr "<<i<<" done. Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   if(!hasConverged)
   {
     cout<<"Fit did not converged in the end"<<endl;
     

      // double minNll(1e20);
      // int minIndex(-1);
      // for(unsigned int j(0); j<fitResVec.size(); ++j)
      // {
      //    if( fitResVec.at(j)->minNll() < minNll)
      //    {
      //       minIndex = j;
      //       minNll = fitResVec.at(j)->minNll();
      //    }
      // }
      

      // migradRes = migradResVec.at(minIndex);
      // hesseRes = hesseResVec.at(minIndex);
      // cout<<"Fit not converged, choose fit "<<minIndex<<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   if (hasConverged) 
   {
     cout<<"Filling tree"<<endl;

     // fillTreeResult(t, fitRes,  update, migradRes, hesseRes, hasConverged);


     TIterator *pariter = par_set->createIterator();
     RooRealVar *pari = (RooRealVar*) pariter->Next();
     while (pari)
     {
       cout<<pari->GetName()<<" value = "<<pari->getVal()<<endl;
       pari = (RooRealVar*) pariter->Next();
       
     }
     

     fillTreeResultSimple(t, par_set,  fitRes,  update, migradRes, hesseRes, hasConverged);
     
   }
   
   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);
   //totPdf.fitTo(*dataGenTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   if (hasConverged) out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(m_wantplot) plot_fit_result(opts.plotsfile, totPdf, dataGenTot,
                                nSignal->getVal(), nPartReco->getVal(), 
                                nComb->getVal()*combIntegral->getVal()*combNorm->getVal(), nJpsiLeak->getVal() );


   // this should remain commented out
   // {
     
   //   //**************Prepare TFile to save the plots
     
   //   TFile f2(plotsfile.c_str(), "UPDATE");
   //   //**************Plot the results of the fit
     
   //   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   //   TIterator *iter = var_set->createIterator();
   //   RooRealVar *var;
     
   //   std::vector<RooPlot*> plots;
   //   RooPlot* frame;
     
     
     
   //   while((var = (RooRealVar*) iter->Next()))
   //   {
       
   //     frame = var->frame();
   //     dataGenTot->plotOn(frame,Name("Data"),DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.6),DrawOption("ZP"));
       
   //     cout<<"here bins"<<endl;
       
   //     totPdf.plotOn(frame, Components("histPdfPartReco"), LineColor(kBlue));
   //     totPdf.plotOn(frame, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
   //     totPdf.plotOn(frame, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
   //     totPdf.plotOn(frame, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
   //     totPdf.plotOn(frame, Components("histPdfJpsiLeak"), LineColor(14));
   //     totPdf.plotOn(frame, Components("combPDF"), LineColor(kBlack));
   //     totPdf.plotOn(frame, LineColor(kRed));
       
   //     plots.push_back(frame);
       
   //   }  
     
   //   if (!(plots.size())) return;
     
   //   TCanvas cFit("cFit", "cFit", 600, 800);
   //   cFit.Divide(1,2);
   //   cFit.cd(1); plots[0]->Draw();
   //   if (plots.size()>1){ 
   //     cFit.cd(2); plots[1]->Draw();
   //   }
     
   //   cFit.Write();
   //   f2.Close();
     
     
   // }
   
   // fw->Close();
   // delete fw;
   //delete and return
   delete nll;
   delete par_set;
   // delete workspace;
   // delete workspaceGen;
   cout<<"deleting binned combinatorial"<<endl;
   delete combPDF_unnorm;
   cout<<"deleted!"<<endl;
   delete combPDF;
   delete gammas;
   delete dataGenTot;
   delete model_pdf;
   delete model_pdf_constrain;
   delete combNorm;
   
   

}


void FitterUtilsHistFact::plot_fit_result(string plotsfile, RooSimultaneous &totPdf, RooDataHist *dataGenTot,
                                          double nsig, double npartreco, double ncomb, double njpsi)
{

   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   RooCategory* idx = (RooCategory*) (&totPdf.indexCat());
   
   double ntot = nsig+ncomb+npartreco+njpsi;
   

   while((var = (RooRealVar*) iter->Next()))
   {

     if (strcmp(var->ClassName(),"RooCategory")==0) continue;
     
      frame = var->frame();
      dataGenTot->plotOn(frame,Name("Data"),DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.6),DrawOption("ZP"));


      totPdf.plotOn(frame,Name("comb"), Slice(*idx),ProjWData(*idx,*dataGenTot),DrawOption("F"),FillColor(kGray),
                             LineWidth(0));


      totPdf.plotOn(frame,Name("comb_line"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                    LineColor(kRed),
                    Components("*combPDF*"),LineWidth(2),
                    Normalization(ntot,RooAbsReal::NumEvent));


      // resids[i]=frame->pullHist();

      totPdf.plotOn(frame,Name("jpsileak"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                    DrawOption("F"),FillColor(kCyan),Components("*Zero*,*One*,*Two*,*PartReco*,*Leak*"),
                    LineWidth(0), Normalization(ntot,RooAbsReal::NumEvent));


      totPdf.plotOn(frame,Name("partreco"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                    DrawOption("F"),FillColor(kMagenta),Components("*Zero*,*One*,*Two*,*PartReco*"),
                    LineWidth(0),Normalization(ntot,RooAbsReal::NumEvent));

      totPdf.plotOn(frame,Name("signal"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                    DrawOption("F"),FillColor(kBlue),Components("*Zero*,*One*,*Two*"),
                    LineWidth(0), Normalization(ntot,RooAbsReal::NumEvent));
      
      totPdf.plotOn(frame,Name("sigzero"), Slice(*idx),ProjWData(*idx,*dataGenTot),Components("*Zero*"),
                    LineWidth(1),LineColor(kYellow),
                    Normalization(ntot,RooAbsReal::NumEvent));

      totPdf.plotOn(frame,Name("sigone"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                    LineColor(kViolet),Components("*One*"),LineWidth(1),
                    Normalization(ntot,RooAbsReal::NumEvent));

      totPdf.plotOn(frame,Name("sigtwo"), Slice(*idx),ProjWData(*idx,*dataGenTot),
                    LineColor(kGreen),Components("*Two*"),LineWidth(1),
                    Normalization(ntot,RooAbsReal::NumEvent));
      
      
      dataGenTot->plotOn(frame,Name("Data"),DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.8),DrawOption("ZP"));
      cout<<"here bins"<<endl;



      // totPdf.plotOn(frame, Components("histPdfPartReco"), LineColor(kBlue));
      // totPdf.plotOn(frame, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
      // totPdf.plotOn(frame, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
      // totPdf.plotOn(frame, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
      // totPdf.plotOn(frame, Components("histPdfJpsiLeak"), LineColor(14));
      // totPdf.plotOn(frame, Components("combPDF"), LineColor(kBlack));
      // totPdf.plotOn(frame, LineColor(kRed));

      plots.push_back(frame);

   }  

   if (!(plots.size())) return;

   TCanvas cFit("cFit", "cFit", 600, 800);
   cFit.Divide(1,2);
   cFit.cd(1); plots[0]->Draw();
   if (plots.size()>1){ 
      cFit.cd(2); plots[1]->Draw();
   }

   cFit.Write();
   f2.Close();


}

void FitterUtilsHistFact::PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
{
   PlotShape2D(originDataSet, genDataSet, shape, plotsfile, canvName, B_plus_M, misPT);
}

void FitterUtilsHistFact::PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
{
   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");

   //**************Plot Signal Zero Gamma

   TH2F* th2fKey = (TH2F*)shape.createHistogram("th2Shape", B_plus_M, Binning(20), YVar(misPT, Binning(20)));
   cout<<genDataSet.sumEntries()<<endl;
   TH2F* th2fGen = (TH2F*)genDataSet.createHistogram("th2fGen", B_plus_M, Binning(20), YVar(misPT, Binning(20)));

   RooPlot* plotM = B_plus_M.frame();
   originDataSet.plotOn(plotM);
   shape.plotOn(plotM);

   RooPlot* plotMisPT = misPT.frame();
   originDataSet.plotOn(plotMisPT);
   shape.plotOn(plotMisPT);

   TCanvas canv(canvName.c_str(), canvName.c_str(), 800, 800);
   canv.Divide(2,2);
   canv.cd(1); th2fGen->Draw("lego");
   canv.cd(2); th2fKey->Draw("surf");
   canv.cd(3); plotM->Draw();
   canv.cd(4); plotMisPT->Draw();

   canv.Write();

   f2.Close();
}


void FitterUtilsHistFact::initiateParams(int nGenSignalZeroGamma, int nGenSignalOneGamma, int nGenSignalTwoGamma, RooRealVar& nSignal, RooRealVar& nPartReco, 
      RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, RooRealVar&  nJpsiLeak, bool constPartReco, RooRealVar const& fracPartRecoSigma, 
      RooRealVar& l1Kee, RooRealVar& l2Kee, RooRealVar& l3Kee, RooRealVar& l4Kee, RooRealVar& l5Kee, 
                                         RooRealVar const& l1KeeGen, RooRealVar const& l2KeeGen, RooRealVar const& l3KeeGen, RooRealVar const& l4KeeGen, RooRealVar const& l5KeeGen , bool constFracs, bool constComb)
{
   TRandom rand;
   rand.SetSeed();

   int nGenSignal = nGenSignalZeroGamma + nGenSignalOneGamma + nGenSignalTwoGamma;

   double nGenSignal2;
   double nGenPartReco2;
   if(!constPartReco)
   {
      nGenSignal2 = rand.Uniform(nGenSignal-5*sqrt(nGenSignal), nGenSignal+5*sqrt(nGenSignal));
      nGenPartReco2 = rand.Uniform(nGenPartReco-5*sqrt(nGenPartReco), nGenPartReco+5*sqrt(nGenPartReco));
   }
   if(constPartReco)
   { 
      double nGenSigPartReco( nGenSignal+nGenPartReco );
      double nGenSigPartReco2( rand.Uniform( nGenSigPartReco-5*sqrt(nGenSigPartReco), nGenSigPartReco+5*sqrt(nGenSigPartReco) ) );
      double fracPartReco1( nGenPartReco/(1.*nGenSignal));
      double fracPartReco2( rand.Uniform(fracPartReco1-5*fracPartRecoSigma.getVal(), fracPartReco1+5*fracPartRecoSigma.getVal()) ); 

      nGenPartReco2 = fracPartReco2*nGenSigPartReco2 / (1+fracPartReco2); 
      nGenSignal2 = nGenSigPartReco2 / (1+fracPartReco2); 
   }
   double nGenComb2 = rand.Uniform(nGenComb-5*sqrt(nGenComb), nGenComb+5*sqrt(nGenComb));
   double nGenJpsiLeak2 = rand.Uniform(nGenJpsiLeak-5*sqrt(nGenJpsiLeak), nGenJpsiLeak+5*sqrt(nGenJpsiLeak));


   nSignal.setVal(nGenSignal2);
   nSignal.setRange(TMath::Max(0.,nGenSignal2-10.*sqrt(nGenSignal)) , nGenSignal2+10*sqrt(nGenSignal));

   nPartReco.setVal(nGenPartReco2);
   nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(nGenPartReco)), nGenPartReco2+10*sqrt(nGenPartReco));


   nComb.setVal(nGenComb2);
   nComb.setRange(TMath::Max(0.,nGenComb2-20.*sqrt(nGenComb)), nGenComb2+20*sqrt(nGenComb));

   nJpsiLeak.setVal(nGenJpsiLeak2);
   nJpsiLeak.setRange(TMath::Max(0., nGenJpsiLeak2-10*sqrt(nGenJpsiLeak)), nGenJpsiLeak2+10*sqrt(nGenJpsiLeak));

   double fracGenZero(nGenSignalZeroGamma/(1.*nGenSignal));
   double fracGenOne(nGenSignalOneGamma/(1.*nGenSignal));

   if (!constFracs)
   {
     
   fracZero.setVal(rand.Gaus(fracGenZero, sqrt(nGenSignalZeroGamma)/(1.*nGenSignal))) ;
   fracZero.setRange(0., 1.);
   fracOne.setVal(rand.Gaus(fracGenOne, sqrt(nGenSignalOneGamma)/(1.*nGenSignal))) ;
   fracOne.setRange(0., 1.);
   }
   

   if (!constComb)
   {
     
   l1Kee.setVal(rand.Uniform( l1KeeGen.getVal() - 5*l1KeeGen.getError(), l1KeeGen.getVal() + 5*l1KeeGen.getError() ) );
   l1Kee.setRange( l1KeeGen.getVal() - 10*l1KeeGen.getError(), l1KeeGen.getVal() + 10*l1KeeGen.getError() );

   l2Kee.setVal(rand.Uniform( l2KeeGen.getVal() - 5*l2KeeGen.getError(), l2KeeGen.getVal() + 5*l2KeeGen.getError() ) );
   l2Kee.setRange( l2KeeGen.getVal() - 10*l2KeeGen.getError(), l2KeeGen.getVal() + 10*l2KeeGen.getError() );

   l3Kee.setVal(rand.Uniform( l3KeeGen.getVal() - 5*l3KeeGen.getError(), l3KeeGen.getVal() + 5*l3KeeGen.getError() ) );
   l3Kee.setRange( l3KeeGen.getVal() - 10*l3KeeGen.getError(), l3KeeGen.getVal() + 10*l3KeeGen.getError() );

   l4Kee.setVal(rand.Uniform( l4KeeGen.getVal() - 5*l4KeeGen.getError(), l4KeeGen.getVal() + 5*l4KeeGen.getError() ) );
   l4Kee.setRange( l4KeeGen.getVal() - 10*l4KeeGen.getError(), l4KeeGen.getVal() + 10*l4KeeGen.getError() );

   l5Kee.setVal(rand.Uniform( l5KeeGen.getVal() - 5*l5KeeGen.getError(), l5KeeGen.getVal() + 5*l5KeeGen.getError() ) );
   l5Kee.setRange( l5KeeGen.getVal() - 10*l5KeeGen.getError(), l5KeeGen.getVal() + 10*l5KeeGen.getError() );
   }

}
