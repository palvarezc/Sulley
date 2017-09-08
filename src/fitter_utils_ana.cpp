#include "fitter_utils_ana.h"
#include "RooBinning.h"
#include "RooRandom.h"
#include "TROOT.h"



FitterUtilsAna::FitterUtilsAna(string name_, const Options &opts_)
  :name(name_),
   m_wantplot(false),
   defaultbinningname("default")
{

  opts = Options(opts_);
  opts.name = name+"_opts";
  combinatorial = Combinatorial(name+"Comb",opts.fittype,opts.plotsfile);
  if (opts.fitSimultaneous) Kemu = Combinatorial(name+"Kemu",opts.fittype,opts.plotsfile);
  TFile fw(opts.workspacefile.c_str(),"RECREATE");
  fw.Close();
}



void FitterUtilsAna::initiateParams(RooArgSet parset, RooArgSet constraints, RooArgSet parset_const)
{

   cout<<"Here!"<<endl;

   TRandom rand;
   rand.SetSeed();

  RooAbsPdf* pdf;
  TIterator *iterc = constraints.createIterator();
  TIterator *iterpc;

  RooArgSet parset_constrained;
  RooRealVar* varc;
  double meanc, sigmac, newvarcval;
  
  // RooRealVar *meancv = NULL;
  // RooRealVar *sigmacv = NULL;


  string namevarc;
  string pdfname;
  string vetoname="fracPartRecoConst";

  while((pdf = (RooAbsPdf*) iterc->Next())) 
    {
      cout<<pdf->GetName()<<endl;
      pdfname=pdf->GetName();
      if (pdfname==vetoname){
        cout<<"Yes"<<endl;
        continue;
      }
      
      iterpc = (pdf->getDependents(parset))->createIterator();
      while(( varc = (RooRealVar*) iterpc->Next()))
	{
	  namevarc = varc->GetName();
	  
	  // meancv = (RooRealVar*) parset_const.find((namevarc+"Mean").c_str());
	  // sigmacv = (RooRealVar*) parset_const.find((namevarc+"Sigma").c_str());

	  // cout<<namevarc<<endl;
	  // if (!meancv || !sigmacv) cout<<"Didnt find constraint parameters"<<endl;

	  meanc = ((RooRealVar*) parset_const.find((namevarc+"Mean").c_str()))->getVal();
	  sigmac = ((RooRealVar*) parset_const.find((namevarc+"Sigma").c_str()))->getVal();
	  

	  // newvarcval = rand.Uniform(TMath::Max(0., meanc-5*sqrt(meanc)), meanc + 5*sqrt(meanc));

	  // varc->setVal(newvarcval);
	  // varc->setRange(TMath::Max(0., meanc-10*sqrt(meanc)), meanc+10*sqrt(meanc));

	  newvarcval = rand.Uniform(meanc-10*sigmac, meanc + 10*sigmac);

	  varc->setVal(newvarcval);
	  // varc->setRange(TMath::Max(0., meanc-10*sqrt(meanc)), meanc+10*sqrt(meanc));
	 
	}
      parset_constrained.add(*(pdf->getDependents(parset)));
    }


  RooRealVar *var;
  TIterator *iter = parset.createIterator();

   while((var = (RooRealVar*) iter->Next()))
   {

     if ( var->isConstant() || parset_constrained.contains(*var) ) continue;
     var->randomize();
   }

   cout<<"Randomized!"<<endl;
}





double FitterUtilsAna::prepare_component_PDF_ee(string filename, string treename, 
					     string compname, string cuts,  
					     RooArgSet varset, RooWorkspace *workspace,
					     string binningName)
{
// This saves both the histograms and the RooHistPDFs
// so it can be used with HistFactory or RooFit

  // gROOT->cd();

   TFile* fComp = new TFile(filename.c_str());
   TTree* tComp = (TTree*)fComp->Get(treename.c_str());
   
   tComp->SetBranchStatus("*",0);

   TIterator *variter = varset.createIterator();
   RooRealVar *var = (RooRealVar*) variter->Next();
   RooRealVar *B_plus_M = (RooRealVar*) varset.find("B_plus_M");
   RooRealVar *misPT = (RooRealVar*) varset.find("misPT");
   
   while (var)
   {
     tComp->SetBranchStatus(var->GetName(), 1);
     var = (RooRealVar*) variter->Next();
   }

   // TH1 *hComp;
   // string histcontent = "";


   // if (opts.fit2D)
   // {
   //   hComp = new TH2D(("h"+compname).c_str(),("h"+compname).c_str(),
   //                         B_plus_M->getBins(binningName.c_str()), B_plus_M->getMin(binningName.c_str()), 
   //                         B_plus_M->getMax(binningName.c_str()),
   //                         misPT->getBins(binningName.c_str()), misPT->getMin(binningName.c_str()), 
   //                         misPT->getMax(binningName.c_str()));

     
   //   hComp->Sumw2();
   //   histcontent = misPT->GetName();
   //   histcontent += ":";
   //   histcontent += B_plus_M->GetName();
   //   histcontent += ">>";
   //   histcontent += hComp->GetName();
   //   tComp->Draw(histcontent.c_str(),("("+cuts+")*("+opts.weightStr+")").c_str(),"qoff");
   // }
   // else
   // {
   //   hComp = new TH1D(("h"+compname).c_str(),("h"+compname).c_str(),
   //                         B_plus_M->getBins(binningName.c_str()), B_plus_M->getMin(binningName.c_str()), 
   //                         B_plus_M->getMax(binningName.c_str()));
   
   //   hComp->Sumw2();
     

   //   histcontent = B_plus_M->GetName();
   //   histcontent += ">>";
   //   histcontent += hComp->GetName();
   //   tComp->Draw(histcontent.c_str(),("("+cuts+")*("+opts.weightStr+")").c_str(),"qoff");
   // }
   
   


   RooDataSet* dataSetComp = new RooDataSet(("dataSet"+compname).c_str(), ("dataSet"+compname).c_str(), varset, 
                                            Import(*tComp), 
                                            Cut(cuts.c_str()), 
                                            WeightVar(opts.weightStr.c_str()));


   cout<<"binning the datasets:"<<endl;

   RooArgSet argset(*B_plus_M);
   if (opts.fit2D) argset.add(*misPT);   

   
   B_plus_M->setBinning(B_plus_M->getBinning(binningName.c_str()));
   if (opts.fit2D) misPT->setBinning(misPT->getBinning(binningName.c_str()));
   

   // RooDataHist dataHistComp(("dataHist"+compname).c_str(), ("dataHist"+compname).c_str(), argset, hComp);
   // RooDataHist dataHistComp(("dataHist"+compname).c_str(), ("dataHist"+compname).c_str(), argset, binningName.c_str());
   // dataHistComp.add(*dataSetComp);   

   RooDataHist dataHistComp(("dataHist"+compname).c_str(), ("dataHist"+compname).c_str(), argset, *dataSetComp);   

   //***************Create 2D histogram estimates from data
   RooHistPdf histPdfComp(("histPdf"+compname).c_str(), ("histPdf"+compname).c_str(), argset, 
                          dataHistComp,opts.interpolate);


   //*************** Compute Error on Component
   cout<<"Number of "<<compname<<": "<< dataSetComp->sumEntries()<<endl;
   double ErrorComp(0);
   if(dataSetComp->sumEntries(cuts.c_str()) > 0) ErrorComp = 1./sqrt(dataSetComp->sumEntries(cuts.c_str()));

   // cout<<"Number of entries in histo: "<<hComp->GetSumOfWeights()<<endl;

   // hComp->SetDirectory(0);

   // workspace->import(*hComp);
   workspace->import(dataHistComp);
   workspace->import(*dataSetComp);
   workspace->import(histPdfComp);

   // delete hComp;
   fComp->Close();
   // delete dataSetComp;
   // delete fComp;

   return ErrorComp;
}


void FitterUtilsAna::prepare_PDFs(int fitmode)
{

  if (fitmode==0) prepare_PDFs_ee();
  if (fitmode==1) prepare_PDFs_mumu();
  if (fitmode==2) 
    {
      prepare_PDFs_ee();
      prepare_PDFs_mumu();
    }
}


void FitterUtilsAna::prepare_PDFs_ee()
{

   //**********Define variables
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", opts.minBMass, opts.maxBMass, "MeV");
   RooRealVar misPT("misPT", "p_{#perp}", opts.minmisPT, opts.maxmisPT, "MeV");
   RooRealVar trigVar(opts.trigStr.c_str(), opts.trigStr.c_str(), -10, 100);
   RooRealVar weightVar(opts.weightStr.c_str(), opts.weightStr.c_str(), -10, 100);
   RooRealVar BDTRooRealVar(opts.BDTVar.c_str(), opts.BDTVar.c_str(), -1,1);
   string BDTVar_forPrc = opts.BDTVar;
   // if (opts.Run==2) BDTVar_forPrc.replace(BDTVar_forPrc.find("Run2"), sizeof("Run2")-1, "Run1"); 
   RooRealVar BDT_forPrcRooRealVar(BDTVar_forPrc.c_str(), BDTVar_forPrc.c_str(), -1,1);
   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV"); 
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);

   //***********Set Binning HARDCODED
   RooBinning defaultMBins(floor((B_plus_M.getMax() - B_plus_M.getMin())/(40.)), B_plus_M.getMin(), B_plus_M.getMax() ); 
   RooBinning defaultMisPTBins(floor(40), misPT.getMin(), misPT.getMax()); 
   RooBinning broaderMBins(floor((B_plus_M.getMax() - B_plus_M.getMin())/(80.)), B_plus_M.getMin(), B_plus_M.getMax()); 
   RooBinning broaderMisPTBins(floor(40), misPT.getMin(), misPT.getMax());

   string broaderbinningname = "broaderBins";
   
   B_plus_M.setBinning( defaultMBins, defaultbinningname.c_str());
   misPT.setBinning( defaultMisPTBins, defaultbinningname.c_str() );
   B_plus_M.setBinning( broaderMBins, broaderbinningname.c_str());
   misPT.setBinning( broaderMisPTBins, broaderbinningname.c_str());
   B_plus_DTFM_M_zero.setBins(100);


   //***********Cuts of each component
   RooArgSet observables(B_plus_M);
   if (opts.fit2D) observables.add(misPT);
   
   RooArgSet argset(observables); 
   argset.add(B_plus_DTFM_M_zero);
   argset.add(trigVar);
   argset.add(e_plus_BremMultiplicity);
   argset.add(BDTRooRealVar);
   RooArgSet argset_comb(argset);
   argset.add(e_minus_BremMultiplicity);

   RooArgSet argset_prc(observables);
   argset_prc.add(B_plus_DTFM_M_zero);
   argset_prc.add(trigVar);
   argset_prc.add(e_plus_BremMultiplicity);
   argset_prc.add(BDT_forPrcRooRealVar);
   argset_prc.add(e_minus_BremMultiplicity);
   argset_prc.add(weightVar);

   RooArgSet argsetPhys(argset);
   argsetPhys.add(weightVar);

   string TrigCutString = "("+opts.trigStr+"  > 0.9)";
   string MassCutString = "(B_plus_M > "+d2s(opts.minBMass)+") && (B_plus_M < "+d2s(opts.maxBMass)+")";
   string BDTCutString = "("+opts.BDTVar+">"+d2s(opts.BDTcut)+")";
   string BDT_forPrcCutString = "("+BDTVar_forPrc+">"+d2s(opts.BDTcut)+")";
   
   string ZeroBremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5)";
   string OneBremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5)";
   string TwoBremCut = "((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5)";

   ZeroBremCut += " && " + TrigCutString + " && " + MassCutString + " && " + BDTCutString;
   OneBremCut += " && " + TrigCutString + " && " + MassCutString + " && " + BDTCutString;
   TwoBremCut += " && " + TrigCutString + " && " + MassCutString + " && " + BDTCutString;
   string BkgCut = TrigCutString + " && " + MassCutString + " && " + BDTCutString;
   string BkgCut_prc = TrigCutString + " && " + MassCutString + " && " + BDT_forPrcCutString;

   // ZeroBremCut += " && " + TrigCutString + " && " + MassCutString;
   // OneBremCut += " && " + TrigCutString + " && " + MassCutString;
   // TwoBremCut += " && " + TrigCutString + " && " + MassCutString;
   // string BkgCut = TrigCutString + " && " + MassCutString;
   

   //***********Create workspace
   // TFile* fw(NULL);
   TFile fw(opts.workspacefile.c_str(),"UPDATE");
   RooWorkspace *workspace = new RooWorkspace("workspace_ee", "workspace_ee");
   workspace->import(B_plus_DTFM_M_zero);
   workspace->import(B_plus_M);
   workspace->import(misPT);


   //***********Add pdfs to workspace
   // double errSignalZero = prepare_component_PDF_ee(opts.signalfile, opts.signaltree,
   //                                              "SignalZeroGamma", ZeroBremCut,  
   //                                              argsetPhys, workspace,
   //                                              defaultbinningname);
   
   // double errSignalOne = prepare_component_PDF_ee(opts.signalfile, opts.signaltree,
   //                                             "SignalOneGamma", OneBremCut,  
   //                                             argsetPhys, workspace,
   //                                             defaultbinningname);

   
   // double errSignalTwo = prepare_component_PDF_ee(opts.signalfile, opts.signaltree,
   //                                             "SignalTwoGamma", TwoBremCut,  
   //                                             argsetPhys, workspace,
   //                                             defaultbinningname);

   
   double errPartReco = prepare_component_PDF_ee(opts.partrecofile, opts.partrecotree,
                                              "PartReco", BkgCut_prc,  
                                              argset_prc, workspace,
                                              defaultbinningname);

   double errJpsiLeak = prepare_component_PDF_ee(opts.JpsiLeakfile, opts.JpsiLeaktree,
                                              "JpsiLeak", BkgCut_prc,  
                                              argset_prc, workspace,
                                              // broaderbinningname);
                                              defaultbinningname);



   RooRealVar fractionalErrorJpsiLeak("fractionalErrorJpsiLeak","fractionalErrorJpsiLeak", errJpsiLeak);
   workspace->import(fractionalErrorJpsiLeak);
   

   cout<<"Error in jpsileak: "<< errJpsiLeak*opts.nGenJpsiLeak<<endl;

   //***************Create combinatorial from fit to data
   combinatorial.build_model_from_data(opts.combfile, opts.combtree, 
                                       BkgCut, &argset_comb, workspace,m_wantplot);


   if (opts.fitSimultaneous)
   {
     RooCategory category("category", "category");
     category.defineType("Kee");
     category.defineType("Kemu");
     workspace->import(category);   
     Kemu.copy_parameters(&combinatorial, workspace, true);
   }
   
   workspace->defineSet("observables",observables);
   // workspace->writeToFile(opts.workspacefile.c_str(),false);
   fw.cd();
   workspace->Write("",TObject::kOverwrite);

   fw.Close();
   delete workspace;

   
   
}


double FitterUtilsAna::prepare_component_PDF_mumu( string filename, string treename, 
					     string compname, string cuts,  
					     RooArgSet varset, RooWorkspace *workspace,
					     string binningName)
{

  bool binnedFit(true); //HARDCODED

  TFile fTrash(filename.c_str());
  TTree* t = (TTree*)fTrash.Get(treename.c_str());

   if(!t)
   {
     cerr<<"ERROR: in function , prepare_component_PDF_mumu no tree found in "<<filename<<endl;
      return 0.;
   }


   t->SetBranchStatus("*",0);

   TIterator *variter = varset.createIterator();
   RooRealVar *var = (RooRealVar*) variter->Next();
   RooRealVar *B_plus_M = (RooRealVar*) varset.find("B_plus_M");
   
   

   while (var)
   {
     t->SetBranchStatus(var->GetName(), 1);
     var = (RooRealVar*) variter->Next();
   }

   if (!B_plus_M)  cerr<<"No B_plus_M!"<<endl;

   // RooDataSet data(("dataSet"+compname).c_str(), ("dataSet"+compname).c_str(), varset, Import(*t), 
   //                                          Cut(cuts.c_str()), 
   //                                          WeightVar(opts.weightStr.c_str()));
   RooDataSet data(("dataSet"+compname).c_str(), ("dataSet"+compname).c_str(), varset, Import(*t), 
                                            Cut(cuts.c_str()));


   workspace->import(data);


   //prepare the fit function

   RooRealVar lambda(("lambda_"+compname).c_str(),("lambda_"+compname).c_str(),-2.39906e-03,-0.1, 0.1);
   RooRealVar meanB(("meanB_"+compname).c_str(),"meanB", 5283., 5280., 5290.);


   RooRealVar sigma(("sigma_"+compname).c_str(), "sigma",16. ,13., 20.);
   RooRealVar sigma1(("sigma1_"+compname).c_str(), "sigma1",28. ,22., 50.);
   RooRealVar n(("n_"+compname).c_str(), "n",20. ,0.5, 50.);
   RooRealVar alpha(("alpha_"+compname).c_str(), "alpha",1.5 ,1., 2.8);

   RooCBShape cb1(("cb1_"+compname).c_str(), "cb1", *B_plus_M,  meanB,  sigma,  alpha, n);
   RooCBShape cb2(("cb2_"+compname).c_str(), "cb2", *B_plus_M,  meanB,  sigma1,  alpha, n);

   RooRealVar vc1(("vc1_"+compname).c_str(), ("vc1_"+compname).c_str(),0.75, 0.15,0.9);
   RooAddPdf model_total(("model_"+compname).c_str(),"model",RooArgList( cb1, cb2),RooArgList( vc1)) ;

   //bin the data to be fast (for testing)

   B_plus_M->setBins(500);  
   RooDataHist* bdata;


   if(binnedFit) bdata = new RooDataHist("data_tau", "data_tau", RooArgSet( *B_plus_M), data) ;

   //fit the stuff


   RooAbsReal* nll;
   
   if(binnedFit) nll = model_total.createNLL(*bdata, RooFit::NumCPU(8) ) ;
   if(!binnedFit) nll = model_total.createNLL(data, RooFit::NumCPU(8) ) ;

   RooMinuit m(*nll);
   m.migrad();
   m.hesse();
   m.minos();


   //plot 
   
   if (m_wantplot)
     {
       RooPlot* frame = B_plus_M->frame();
       data.plotOn(frame);
       model_total.plotOn(frame, RooFit::LineColor(kRed) );
       model_total.plotOn(frame, RooFit::LineColor(kGreen), RooFit::Components("expo"));
       model_total.plotOn(frame, RooFit::LineColor(kBlue), RooFit::Components("sigCB0"));

       TCanvas canv("canv", "canv", 600, 600);
       frame->Draw();
       canv.Print((opts.plotsdir+"/plotKmumucalibration_"+compname+".pdf").c_str());
       canv.Print((opts.plotsdir+"/plotKmumucalibration_"+compname+".root").c_str());
     }
   //save the fit function

//   workspace->import(lambda, true);
//   workspace->import(meanB, true);
//   workspace->import(sigma, true );
//   workspace->import(sigma1, true);
//   workspace->import(n, true);
//   workspace->import(alpha, true);
//   //workspace->import(cb1);
//   //workspace->import(cb2);
//   //workspace->import(expo);
//   workspace->import(vc1, true);
//   workspace->import(sig, true);
//   workspace->import(bkg, true);
   //workspace->import(sigCB_B0);

   RooArgSet *par_set = model_total.getParameters(data);
   workspace->import(lambda);
   workspace->import(model_total);
   workspace->saveSnapshot("pars_for_generation", *par_set);


   // cout<<"************ Fit parameters: "<<endl;
   // cout<<"meanB: "<<meanB.getVal()<<endl;
   // cout<<"sigma: "<<sigma.getVal()<<endl;
   // cout<<"sigma1: "<<sigma1.getVal()<<endl;
   // cout<<"n: "<<n.getVal()<<endl;
   // cout<<"alpha: "<<alpha.getVal()<<endl;
   // cout<<"vc1: "<<vc1.getVal()<<endl;
   // cout<<"lambda: "<<lambda.getVal()<<endl;


   delete par_set;
   delete nll;
   return  1.;
}





void FitterUtilsAna::prepare_PDFs_mumu()
{

   //**********Define variables
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", opts.minBMass_mumu, opts.maxBMass_mumu, "MeV");
   RooRealVar trigVar(opts.trigStr.c_str(), opts.trigStr.c_str(), -10, 100);
   RooRealVar weightVar(opts.weightStr.c_str(), opts.weightStr.c_str(), -10, 10);


   //***********Cuts of each component
   RooArgSet observables(B_plus_M);
   
   RooArgSet argset(observables); 
   argset.add(trigVar);
   
   RooArgSet argsetPhys(argset);
   // argsetPhys.add(weightVar);

   string TrigCutString = "("+opts.trigStr+"  > 0.9)";
   string MassCutString = "(B_plus_M > "+d2s(opts.minBMass_mumu)+") && (B_plus_M < "+d2s(opts.maxBMass_mumu)+")";
   string BkgCut = TrigCutString + " && " + MassCutString;
   

   //***********Create workspace
   // TFile* fw(NULL);
   TFile fw(opts.workspacefile.c_str(),"UPDATE");
   RooWorkspace *workspace = new RooWorkspace("workspace_mumu", "workspace_mumu");
   // workspace->import(B_plus_M);

   //***********Add pdfs to workspace
   double errSignal = prepare_component_PDF_mumu(opts.signalfile_mumu, opts.signaltree_mumu,
   						 "Signal_mumu", BkgCut,  
   						 argsetPhys, workspace,
   						 "default");   


   workspace->defineSet("observables",observables);
   fw.cd();
   workspace->Write("",TObject::kOverwrite);

   fw.Close();
   delete workspace;

   
}







void FitterUtilsAna::generate_component(string compname, RooArgSet varset, int nGenComp, 
				     RooWorkspace* workspace, RooWorkspace* workspaceGen)
{

  
  RooHistPdf *histPdfComp = (RooHistPdf *) workspace->pdf(("histPdf"+compname).c_str());
  RooAbsPdf::GenSpec* GenSpecComp = histPdfComp->prepareMultiGen(varset, RooFit::Extended(1), NumEvents(nGenComp), AutoBinned(false)); 
  
  cout<<"Generating "<<compname<<endl;
  RooDataSet* dataGenComp = histPdfComp->generate(*GenSpecComp);
  dataGenComp->SetName(("dataGen"+compname).c_str()); dataGenComp->SetTitle(("dataGen"+compname).c_str());

  workspaceGen->import(*dataGenComp);


  if (m_wantplot) {
    RooDataSet* dataSetComp = (RooDataSet*)workspace->data(("dataSet"+compname).c_str());
    PlotShape(*dataSetComp, *dataGenComp, *histPdfComp, opts.plotsfile, ("c"+compname).c_str(), varset);
  }
  
  delete dataGenComp;
  delete GenSpecComp;
  
}


void FitterUtilsAna::generate_component_ana(string compname, RooArgSet varset, int nGenComp, 
                                            RooWorkspace* workspaceMCFit, RooWorkspace* workspaceGen,
                                            int trigCat, int PhotCat)
{

  

  string cat("Trig"+i2s(trigCat)+"Phot"+i2s(PhotCat));
  string trigCatS("Trig"+i2s(trigCat));
  string trigCatSPrc("Trig"+i2s(trigCat)+"Phot-1");
  
  RooRealVar *B_plus_M = (RooRealVar*) varset.find("B_plus_M");

   RooRealVar* mean = (RooRealVar*)workspaceMCFit->var(("mean"+cat).c_str());
   cout<<"Hello3: "<<mean<<endl;
   mean->setConstant(true);
   RooRealVar* mean2 = (RooRealVar*)workspaceMCFit->var(("mean2"+cat).c_str());
   cout<<"Hello3: "<<mean2<<endl;
   mean2->setConstant(true);
   RooRealVar* mean3 = (RooRealVar*)workspaceMCFit->var(("mean3"+cat).c_str());
   cout<<"Hello3: "<<mean3<<endl;
   mean3->setConstant(true);
   RooRealVar* mean4 = (RooRealVar*)workspaceMCFit->var(("mean4"+cat).c_str());
   cout<<"Hello3: "<<mean4<<endl;
   mean4->setConstant(true);
   RooRealVar* sigma = (RooRealVar*)workspaceMCFit->var(("sigma"+cat).c_str());
   cout<<"Hello4: "<<sigma<<endl;
   sigma->setConstant(true);
   RooRealVar* sigma2 = (RooRealVar*)workspaceMCFit->var(("sigma2"+cat).c_str());
   cout<<"Hello4: "<<sigma2<<endl;
   sigma2->setConstant(true);
   RooRealVar* sigma3 = (RooRealVar*)workspaceMCFit->var(("sigma3"+cat).c_str());
   cout<<"Hello4: "<<sigma3<<endl;
   sigma3->setConstant(true);
   RooRealVar* sigma4 = (RooRealVar*)workspaceMCFit->var(("sigma4"+cat).c_str());
   cout<<"Hello4: "<<sigma4<<endl;
   sigma4->setConstant(true);
//   RooRealVar* fracSigma0 = (RooRealVar*)workspaceMCFit0->var(("fracSigma"+cat0).c_str());
//   cout<<"Hello5: "<<fracSigma0<<endl;
//   fracSigma0->setConstant(true);
   RooRealVar* al = (RooRealVar*)workspaceMCFit->var(("alpha_left"+cat).c_str());
   cout<<"Hello6: "<<al<<endl;
   al->setConstant(true);
   RooRealVar* ar = (RooRealVar*)workspaceMCFit->var(("alpha_right"+cat).c_str());
   cout<<"Hello7: "<<ar<<endl;
   ar->setConstant(true);
   RooRealVar* nl = (RooRealVar*)workspaceMCFit->var(("n_left"+cat).c_str());
   cout<<"Hello8: "<<nl<<endl;
   nl->setConstant(true);
   RooRealVar* nr = (RooRealVar*)workspaceMCFit->var(("n_right"+cat).c_str());
   cout<<"Hello9: "<<nr<<endl;
   nr->setConstant(true);
   RooRealVar* fracCB2 = (RooRealVar*)workspaceMCFit->var(("fracCB2"+cat).c_str());
   cout<<"Hello1: "<<fracCB2<<endl;
   fracCB2->setConstant(true);
   RooRealVar* fracGaus = (RooRealVar*)workspaceMCFit->var(("fracGaus"+cat).c_str());
   cout<<"Hello1: "<<fracGaus<<endl;
   fracGaus->setConstant(true);
   RooRealVar* fracGaus2 = (RooRealVar*)workspaceMCFit->var(("fracGaus2"+cat).c_str());
   cout<<"Hello1: "<<fracGaus2<<endl;
   fracGaus2->setConstant(true);
   RooRealVar* fracCB3 = (RooRealVar*)workspaceMCFit->var(("fracCB3"+cat).c_str());
   cout<<"Hello1: "<<fracCB3<<endl;
   fracCB3->setConstant(true);

  
   // RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+cat).c_str(), ("sigmaScaleFactor"+cat).c_str(), 1.);
   // RooRealVar meanShift(("meanShift"+cat).c_str(), ("meanShift"+cat).c_str(), 0.);
   RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+cat).c_str(), ("sigmaScaleFactor"+cat).c_str(), 1.1931);//1.12, 0.9, 1.5);
   RooRealVar meanShift(("meanShift"+cat).c_str(), ("meanShift"+cat).c_str(), 5.4707);//0, -100, 100);

   if (!opts.massshift)
   {
     sigmaScaleFactor.setVal(1.);
     meanShift.setVal(0.);
   }

   RooFormulaVar shiftedMean(("shiftedMean"+cat).c_str(), ("shiftedMean"+cat).c_str(), "@0+@1", RooArgList(*mean, meanShift));
   RooFormulaVar shiftedMean2(("shiftedMean2"+cat).c_str(), ("shiftedMean2"+cat).c_str(), "@0+@1", RooArgList(*mean2, meanShift));
   RooFormulaVar shiftedMean3(("shiftedMean3"+cat).c_str(), ("shiftedMean3"+cat).c_str(), "@0+@1", RooArgList(*mean3, meanShift));
   RooFormulaVar shiftedMean4(("shiftedMean4"+cat).c_str(), ("shiftedMean4"+cat).c_str(), "@0+@1", RooArgList(*mean4, meanShift));
   RooFormulaVar scaledSigma(("scaledSigma"+cat).c_str(), ("scaledSigma"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma));
   RooFormulaVar scaledSigma2(("scaledSigma2"+cat).c_str(), ("scaledSigma2"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma2));
   RooFormulaVar scaledSigma4(("scaledSigma4"+cat).c_str(), ("scaledSigma4"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma4));
   RooFormulaVar scaledSigma3(("scaledSigma3"+cat).c_str(), ("scaledSigma3"+cat).c_str(), "@0*@1", RooArgList(sigmaScaleFactor, *sigma3));
   
   RooRealVar arScale("arScale", "arScale", 1,0,3);
   arScale.setConstant(true);
   RooFormulaVar scaledAr(("scaledAr"+cat).c_str(), ("scaledAr"+cat).c_str(), "@0*@1", RooArgList(*ar, arScale));

   RooCBShape cb1(("cb1"+cat).c_str(), ("cb1"+cat).c_str(), *B_plus_M, shiftedMean, scaledSigma, *al, *nl);
   // RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, shiftedMean, scaledSigma2, scaledAr, *nr);
   RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, shiftedMean2, scaledSigma2, scaledAr, *nr);
   
   RooGaussian gaus(("gaus"+cat).c_str(), ("gaus"+cat).c_str(), *B_plus_M, shiftedMean3, scaledSigma3 );
   RooGaussian gaus2(("gaus2"+cat).c_str(), ("gaus2"+cat).c_str(), *B_plus_M, shiftedMean4, scaledSigma4 );
  
   RooAbsPdf *modelComp;

   if (PhotCat==0)
   {
     
     modelComp = new RooCBShape(("modelComp"+cat).c_str(), ("modelComp"+cat).c_str(), *B_plus_M, shiftedMean, scaledSigma,*al, *nl);
     }
     if (PhotCat==1)
       // if (PhotCat==1 || PhotCat==2)
     {
       modelComp = new RooAddPdf(("modelComp"+cat).c_str(), ("modelComp"+cat).c_str(), RooArgList(cb1, gaus) ,
                                   RooArgList(*fracCB2), true);
     }
     if (PhotCat==2)
     {
       modelComp = new RooAddPdf(("modelComp"+cat).c_str(), ("modelComp"+cat).c_str(), RooArgList(cb1, gaus),
                                   RooArgList(*fracCB2),true);
     }


  // RooHistPdf *modelComp = (RooHistPdf *) workspace->pdf(("histPdf"+compname).c_str());
  RooAbsPdf::GenSpec* GenSpecComp = modelComp->prepareMultiGen(varset, RooFit::Extended(1), NumEvents(nGenComp), AutoBinned(false)); 
  
  cout<<"Generating "<<compname<<endl;
  RooDataSet* dataGenComp = modelComp->generate(*GenSpecComp);
  dataGenComp->SetName(("dataGen"+compname).c_str()); dataGenComp->SetTitle(("dataGen"+compname).c_str());

  workspaceGen->import(*dataGenComp);


  if (m_wantplot) {
    // RooDataSet* dataSetComp = (RooDataSet*)workspace->data(("dataSet"+compname).c_str());
    // PlotShape(*dataSetComp, *dataGenComp, *modelComp, opts.plotsfile, ("c"+compname).c_str(), varset);
    RooDataSet* dataSetComp = (RooDataSet*)workspaceMCFit->data("data");
    PlotShape(*dataSetComp, *dataGenComp, *modelComp, opts.plotsfile, ("c"+compname).c_str(), varset);
  }
  
  delete modelComp;
  delete dataGenComp;
  delete GenSpecComp;
  
}



void FitterUtilsAna::generate(int fitmode)
{
  if (fitmode==0) generate_ee();
  if (fitmode==1) generate_mumu();
  if (fitmode==2) 
    {
      generate_ee();
      generate_mumu();
    }

}




void FitterUtilsAna::generate_ee()
{

  // HARDCODED!!
  string workspaceMCFitName = "/vols/lhcb/palvare1/RK_analysis/FitWorkspaces/WorkspaceKee_Run"+to_string(opts.Run)+".root";
  string cat0("Trig"+i2s(opts.trigCat)+"Phot0");
  string cat1("Trig"+i2s(opts.trigCat)+"Phot1");
  string cat2("Trig"+i2s(opts.trigCat)+"Phot2");
  TFile fwMCFit(workspaceMCFitName.c_str());
  RooWorkspace* workspaceMCFit0 = (RooWorkspace*)fwMCFit.Get(("workspaceMCFit"+cat0).c_str());
  RooWorkspace* workspaceMCFit1 = (RooWorkspace*)fwMCFit.Get(("workspaceMCFit"+cat1).c_str());
  RooWorkspace* workspaceMCFit2 = (RooWorkspace*)fwMCFit.Get(("workspaceMCFit"+cat2).c_str());

  

   //***************Get the PDFs from the workspace
   TFile fw(opts.workspacefile.c_str(), "UPDATE");   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace_ee");
   RooArgSet *argset = (RooArgSet*) workspace->set("observables");

   //***************Generate some datasets and save
   RooWorkspace *workspaceGen = new RooWorkspace("workspaceGen_ee", "workspaceGen_ee");
   generate_component_ana("SignalZeroGamma", *argset, opts.nGenSignalZeroGamma, workspaceMCFit0, workspaceGen, opts.trigCat, 0);
   generate_component_ana("SignalOneGamma", *argset, opts.nGenSignalOneGamma, workspaceMCFit1, workspaceGen, opts.trigCat, 1);
   generate_component_ana("SignalTwoGamma", *argset, opts.nGenSignalTwoGamma, workspaceMCFit2, workspaceGen, opts.trigCat, 2);
   generate_component("PartReco", *argset, opts.nGenPartReco, workspace, workspaceGen);
   if(opts.nGenJpsiLeak>0) generate_component("JpsiLeak", *argset, opts.nGenJpsiLeak, workspace, workspaceGen);
   
   fw.cd();
   workspaceGen->Write("", TObject::kOverwrite);
   
   combinatorial.generate("Comb", argset,  opts.nGenComb, workspace, workspaceGen, m_wantplot);

   
   if (opts.fitSimultaneous) Kemu.generate("Kemu", argset, opts.nGenKemu, workspace, workspaceGen, m_wantplot);   
   
   fw.cd();
   workspaceGen->Write("", TObject::kOverwrite);
   fw.Close();
   delete workspaceGen;

}


void FitterUtilsAna::generate_mumu()
{
   //***************Get the PDFs parameters from the workspace

   TFile fw(opts.workspacefile.c_str(), "UPDATE");   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace_mumu");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");

   string compname = "Signal_mumu";
   RooRealVar *meanB = workspace->var(("meanB_"+compname).c_str()); 
   RooRealVar *sigma = workspace->var(("sigma_"+compname).c_str()); 
   RooRealVar *sigma1 = workspace->var(("sigma1_"+compname).c_str()); 
   RooRealVar *n = workspace->var(("n_"+compname).c_str()); 
   RooRealVar *alpha = workspace->var(("alpha_"+compname).c_str()); 
   RooRealVar *vc1 = workspace->var(("vc1_"+compname).c_str()); 
   RooRealVar *lambda = workspace->var(("lambda_"+compname).c_str()); 

   workspace->loadSnapshot("pars_for_generation");
   cout<<"************ Generation parameters: "<<endl;
   cout<<"meanB: "<<meanB->getVal()<<endl;
   cout<<"sigma: "<<sigma->getVal()<<endl;
   cout<<"sigma1: "<<sigma1->getVal()<<endl;
   cout<<"n: "<<n->getVal()<<endl;
   cout<<"alpha: "<<alpha->getVal()<<endl;
   cout<<"vc1: "<<vc1->getVal()<<endl;
   cout<<"lambda: "<<lambda->getVal()<<endl;

   
   //*****************Build pdfs
   RooCBShape cb1(("cb1_"+compname).c_str(), "cb1", *B_plus_M,  *meanB,  *sigma,  *alpha, *n);
   RooCBShape cb2(("cb2_"+compname).c_str(), "cb2", *B_plus_M,  *meanB,  *sigma1,  *alpha, *n);

   RooAddPdf model_signal(("model_"+compname).c_str(),"model",RooArgList( cb1, cb2),RooArgList( *vc1)) ;

   RooExponential expo("expo", "exponential PDF", *B_plus_M,  *lambda);

   //***************Generate some datasets and save
   RooArgSet varset(*B_plus_M);

   cout<<"Requested to generate: "<<opts.nGenSignal_mumu<<" signal and "<<opts.nGenBkg_mumu<<" bkg."<<endl;

   RooAbsPdf::GenSpec* GenSpecSignal = model_signal.prepareMultiGen(varset, Extended(1), NumEvents(opts.nGenSignal_mumu)); 
   // RooDataSet* dataGenSignal = model_signal.generate(*GenSpecSignal);
   RooDataSet* dataGenSignal = model_signal.generate(varset, Extended(1), NumEvents(opts.nGenSignal_mumu));
   dataGenSignal->SetName("dataGenSignal_mumu"); dataGenSignal->SetTitle("dataGenSignal_mumu");


   RooAbsPdf::GenSpec* GenSpecBkg = expo.prepareMultiGen(varset, Extended(1), NumEvents(opts.nGenBkg_mumu)); 
   RooDataSet* dataGenBkg = expo.generate(*GenSpecBkg);
   dataGenBkg->SetName("dataGenBkg_mumu"); dataGenBkg->SetTitle("dataGenBkg_mumu");

   cout<<"Generated: "<<dataGenSignal->sumEntries()<<" signal and "<<dataGenBkg->sumEntries()<<" bkg."<<endl;

   RooWorkspace *workspaceGen = new RooWorkspace("workspaceGen_mumu", "workspaceGen_mumu");
   workspaceGen->import(*dataGenSignal);
   workspaceGen->import(*dataGenBkg);
   workspaceGen->Write("", TObject::kOverwrite);
   //delete workspace;
   fw.Close();

}



void FitterUtilsAna::run_toy(double fracPartReco_const,
			     ofstream& out, TTree* t, bool update,
			     int fitmode)
{

  //***************Generate
  generate(fitmode);

  if (fitmode==0)
    {
      TFile fw(opts.workspacefile.c_str());

      //***************Retrieve generated samples electrons
      RooWorkspace* workspaceGen_ee = (RooWorkspace*)fw.Get("workspaceGen_ee");
      RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen_ee->data("dataGenSignalZeroGamma");
      RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen_ee->data("dataGenSignalOneGamma");
      RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen_ee->data("dataGenSignalTwoGamma");
      RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen_ee->data("dataGenPartReco");
      RooDataSet* dataGenComb = (RooDataSet*)workspaceGen_ee->data("dataGenComb");
      RooDataSet* dataGenJpsiLeak(0);
      RooDataSet* dataGenKemu(0);
      if(opts.nGenJpsiLeak>0) dataGenJpsiLeak = (RooDataSet*)workspaceGen_ee->data("dataGenJpsiLeak");
      if(opts.fitSimultaneous) dataGenKemu = (RooDataSet*)workspaceGen_ee->data("dataGenKemu");
      
      //***************Merge datasets
      RooDataSet* dataGenTot_ee(dataGenPartReco);
      dataGenTot_ee->append(*dataGenSignalZeroGamma);
      dataGenTot_ee->append(*dataGenSignalOneGamma);
      dataGenTot_ee->append(*dataGenSignalTwoGamma);
      dataGenTot_ee->append(*dataGenComb);
      if(opts.nGenJpsiLeak>0) dataGenTot_ee->append(*dataGenJpsiLeak);
      
      RooWorkspace* workspace_ee = (RooWorkspace*)fw.Get("workspace_ee");      
      fit_ee(dataGenTot_ee, dataGenKemu, workspace_ee, fracPartReco_const, out, t, update);

      fw.Close();
    }

  if (fitmode==1)
    {
  
      TFile fw(opts.workspacefile.c_str());

      //***************Retrieve generated samples muons
      RooWorkspace* workspaceGen_mumu = (RooWorkspace*)fw.Get("workspaceGen_mumu");
      RooDataSet* dataGenSignal_mumu = (RooDataSet*) workspaceGen_mumu->data("dataGenSignal_mumu"); 
      RooDataSet* dataGenBkg_mumu = (RooDataSet*) workspaceGen_mumu->data("dataGenBkg_mumu"); 

      //***************Merge datasets
      RooDataSet* dataGenTot_mumu(dataGenSignal_mumu);
      dataGenTot_mumu->append(*dataGenBkg_mumu);

      RooWorkspace* workspace_mumu = (RooWorkspace*)fw.Get("workspace_mumu");      
      fit_mumu(dataGenTot_mumu, workspace_mumu, fracPartReco_const, out, t, update);

      fw.Close();
    }


  if (fitmode==2)
    {

      TFile fw(opts.workspacefile.c_str());

      //***************Retrieve generated samples electrons
      RooWorkspace* workspaceGen_ee = (RooWorkspace*)fw.Get("workspaceGen_ee");
      RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen_ee->data("dataGenSignalZeroGamma");
      RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen_ee->data("dataGenSignalOneGamma");
      RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen_ee->data("dataGenSignalTwoGamma");
      RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen_ee->data("dataGenPartReco");
      RooDataSet* dataGenComb = (RooDataSet*)workspaceGen_ee->data("dataGenComb");
      RooDataSet* dataGenJpsiLeak(0);
      RooDataSet* dataGenKemu(0);
      if(opts.nGenJpsiLeak>0) dataGenJpsiLeak = (RooDataSet*)workspaceGen_ee->data("dataGenJpsiLeak");
      if(opts.fitSimultaneous) dataGenKemu = (RooDataSet*)workspaceGen_ee->data("dataGenKemu");
      
      //***************Merge datasets electrons
      RooDataSet* dataGenTot_ee(dataGenPartReco);
      dataGenTot_ee->append(*dataGenSignalZeroGamma);
      dataGenTot_ee->append(*dataGenSignalOneGamma);
      dataGenTot_ee->append(*dataGenSignalTwoGamma);
      dataGenTot_ee->append(*dataGenComb);
      if(opts.nGenJpsiLeak>0) dataGenTot_ee->append(*dataGenJpsiLeak);
      
      //***************Retrieve generated samples muons
      RooWorkspace* workspaceGen_mumu = (RooWorkspace*)fw.Get("workspaceGen_mumu");
      RooDataSet* dataGenSignal_mumu = (RooDataSet*) workspaceGen_mumu->data("dataGenSignal_mumu"); 
      RooDataSet* dataGenBkg_mumu = (RooDataSet*) workspaceGen_mumu->data("dataGenBkg_mumu"); 

      //***************Merge datasets
      RooDataSet* dataGenTot_mumu(dataGenSignal_mumu);
      dataGenTot_mumu->append(*dataGenBkg_mumu);

      RooWorkspace* workspace_mumu = (RooWorkspace*)fw.Get("workspace_mumu");      
      RooWorkspace* workspace_ee = (RooWorkspace*)fw.Get("workspace_ee");      

      fit_RK(dataGenTot_ee, dataGenTot_mumu, dataGenKemu, workspace_ee, workspace_mumu, fracPartReco_const, out, t, update);

    }

}

RooAbsPdf* FitterUtilsAna::build_component_ana(string compname, RooArgSet varset, 
                                               RooWorkspace* workspaceMCFit, 
                                               int trigCat, int PhotCat,
                                               RooRealVar* meanShift, RooRealVar* sigmaScaleFactor)
{

  

  string cat("Trig"+i2s(trigCat)+"Phot"+i2s(PhotCat));
  string trigCatS("Trig"+i2s(trigCat));
  string trigCatSPrc("Trig"+i2s(trigCat)+"Phot-1");
  
  RooRealVar *B_plus_M = (RooRealVar*) varset.find("B_plus_M");

   RooRealVar* mean = (RooRealVar*)workspaceMCFit->var(("mean"+cat).c_str());
   cout<<"Hello3: "<<mean<<endl;
   mean->setConstant(true);
   RooRealVar* mean2 = (RooRealVar*)workspaceMCFit->var(("mean2"+cat).c_str());
   cout<<"Hello3: "<<mean2<<endl;
   mean2->setConstant(true);
   RooRealVar* mean3 = (RooRealVar*)workspaceMCFit->var(("mean3"+cat).c_str());
   cout<<"Hello3: "<<mean3<<endl;
   mean3->setConstant(true);
   RooRealVar* mean4 = (RooRealVar*)workspaceMCFit->var(("mean4"+cat).c_str());
   cout<<"Hello3: "<<mean4<<endl;
   mean4->setConstant(true);
   RooRealVar* sigma = (RooRealVar*)workspaceMCFit->var(("sigma"+cat).c_str());
   cout<<"Hello4: "<<sigma<<endl;
   sigma->setConstant(true);
   RooRealVar* sigma2 = (RooRealVar*)workspaceMCFit->var(("sigma2"+cat).c_str());
   cout<<"Hello4: "<<sigma2<<endl;
   sigma2->setConstant(true);
   RooRealVar* sigma3 = (RooRealVar*)workspaceMCFit->var(("sigma3"+cat).c_str());
   cout<<"Hello4: "<<sigma3<<endl;
   sigma3->setConstant(true);
   RooRealVar* sigma4 = (RooRealVar*)workspaceMCFit->var(("sigma4"+cat).c_str());
   cout<<"Hello4: "<<sigma4<<endl;
   sigma4->setConstant(true);
//   RooRealVar* fracSigma0 = (RooRealVar*)workspaceMCFit0->var(("fracSigma"+cat0).c_str());
//   cout<<"Hello5: "<<fracSigma0<<endl;
//   fracSigma0->setConstant(true);
   RooRealVar* al = (RooRealVar*)workspaceMCFit->var(("alpha_left"+cat).c_str());
   cout<<"Hello6: "<<al<<endl;
   al->setConstant(true);
   RooRealVar* ar = (RooRealVar*)workspaceMCFit->var(("alpha_right"+cat).c_str());
   cout<<"Hello7: "<<ar<<endl;
   ar->setConstant(true);
   RooRealVar* nl = (RooRealVar*)workspaceMCFit->var(("n_left"+cat).c_str());
   cout<<"Hello8: "<<nl<<endl;
   nl->setConstant(true);
   RooRealVar* nr = (RooRealVar*)workspaceMCFit->var(("n_right"+cat).c_str());
   cout<<"Hello9: "<<nr<<endl;
   nr->setConstant(true);
   RooRealVar* fracCB2 = (RooRealVar*)workspaceMCFit->var(("fracCB2"+cat).c_str());
   cout<<"Hello1: "<<fracCB2<<endl;
   fracCB2->setConstant(true);
   RooRealVar* fracGaus = (RooRealVar*)workspaceMCFit->var(("fracGaus"+cat).c_str());
   cout<<"Hello1: "<<fracGaus<<endl;
   fracGaus->setConstant(true);
   RooRealVar* fracGaus2 = (RooRealVar*)workspaceMCFit->var(("fracGaus2"+cat).c_str());
   cout<<"Hello1: "<<fracGaus2<<endl;
   fracGaus2->setConstant(true);
   RooRealVar* fracCB3 = (RooRealVar*)workspaceMCFit->var(("fracCB3"+cat).c_str());
   cout<<"Hello1: "<<fracCB3<<endl;
   fracCB3->setConstant(true);

  
   // RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+cat).c_str(), ("sigmaScaleFactor"+cat).c_str(), 1.);
   // RooRealVar meanShift(("meanShift"+cat).c_str(), ("meanShift"+cat).c_str(), 0.);
   // RooRealVar sigmaScaleFactor(("sigmaScaleFactor"+cat).c_str(), ("sigmaScaleFactor"+cat).c_str(), 1.1931);//1.12, 0.9, 1.5);
   // RooRealVar meanShift(("meanShift"+cat).c_str(), ("meanShift"+cat).c_str(), 5.4707);//0, -100, 100);

   RooFormulaVar* shiftedMean = new RooFormulaVar(("shiftedMean"+cat).c_str(), ("shiftedMean"+cat).c_str(), "@0+@1", RooArgList(*mean, *meanShift));
   RooFormulaVar* shiftedMean2 = new RooFormulaVar(("shiftedMean2"+cat).c_str(), ("shiftedMean2"+cat).c_str(), "@0+@1", RooArgList(*mean2, *meanShift));
   RooFormulaVar* shiftedMean3 = new RooFormulaVar(("shiftedMean3"+cat).c_str(), ("shiftedMean3"+cat).c_str(), "@0+@1", RooArgList(*mean3, *meanShift));
   RooFormulaVar* shiftedMean4 = new RooFormulaVar(("shiftedMean4"+cat).c_str(), ("shiftedMean4"+cat).c_str(), "@0+@1", RooArgList(*mean4, *meanShift));
   RooFormulaVar* scaledSigma = new RooFormulaVar(("scaledSigma"+cat).c_str(), ("scaledSigma"+cat).c_str(), "@0*@1", RooArgList(*sigmaScaleFactor, *sigma));
   RooFormulaVar* scaledSigma2 = new RooFormulaVar(("scaledSigma2"+cat).c_str(), ("scaledSigma2"+cat).c_str(), "@0*@1", RooArgList(*sigmaScaleFactor, *sigma2));
   RooFormulaVar* scaledSigma4 = new RooFormulaVar(("scaledSigma4"+cat).c_str(), ("scaledSigma4"+cat).c_str(), "@0*@1", RooArgList(*sigmaScaleFactor, *sigma4));
   RooFormulaVar* scaledSigma3 = new RooFormulaVar(("scaledSigma3"+cat).c_str(), ("scaledSigma3"+cat).c_str(), "@0*@1", RooArgList(*sigmaScaleFactor, *sigma3));
   
   // RooRealVar arScale("arScale", "arScale", 1,0,3);
   // arScale.setConstant(true);
   // RooFormulaVar scaledAr(("scaledAr"+cat).c_str(), ("scaledAr"+cat).c_str(), "@0*@1", RooArgList(*ar, arScale));

   RooCBShape *cb1 = new RooCBShape(("cb1"+cat).c_str(), ("cb1"+cat).c_str(), *B_plus_M, *shiftedMean, *scaledSigma, *al, *nl);
   // RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, *shiftedMean, *scaledSigma2, scaledAr, *nr);
   RooCBShape cb2(("cb2"+cat).c_str(), ("cb2"+cat).c_str(), *B_plus_M, *shiftedMean2, *scaledSigma2, *ar, *nr);
   
   RooGaussian *gaus = new RooGaussian(("gaus"+cat).c_str(), ("gaus"+cat).c_str(), *B_plus_M, *shiftedMean3, *scaledSigma3 );
   RooGaussian gaus2(("gaus2"+cat).c_str(), ("gaus2"+cat).c_str(), *B_plus_M, *shiftedMean4, *scaledSigma4 );
  
   RooAbsPdf *modelComp;

   if (PhotCat==0)
   {
     
     modelComp = new RooCBShape((compname).c_str(), (compname).c_str(), *B_plus_M, *shiftedMean, *scaledSigma,*al, *nl);
     }
     if (PhotCat==1)
       // if (PhotCat==1 || PhotCat==2)
     {
       modelComp = new RooAddPdf((compname).c_str(), (compname).c_str(), RooArgList(*cb1, *gaus) ,
                                   RooArgList(*fracCB2), true);
     }
     if (PhotCat==2)
     {
       modelComp = new RooAddPdf((compname).c_str(), (compname).c_str(), RooArgList(*cb1, *gaus),
                                   RooArgList(*fracCB2),true);
     }


     return modelComp;  
}


void FitterUtilsAna::fit_ee(RooDataSet* dataset_ee, RooDataSet* dataset_Kemu, 
			 RooWorkspace* workspace, 
			 double fracPartReco_const,
			 ofstream& out, TTree* t, bool update)
{

  // HARDCODED!!
  string workspaceMCFitName = "/vols/lhcb/palvare1/RK_analysis/FitWorkspaces/WorkspaceKee_Run"+to_string(opts.Run)+".root";
  string cat0("Trig"+i2s(opts.trigCat)+"Phot0");
  string cat1("Trig"+i2s(opts.trigCat)+"Phot1");
  string cat2("Trig"+i2s(opts.trigCat)+"Phot2");
  string cat("Trig"+i2s(opts.trigCat));
  TFile fwMCFit(workspaceMCFitName.c_str());
  RooWorkspace* workspaceMCFit0 = (RooWorkspace*)fwMCFit.Get(("workspaceMCFit"+cat0).c_str());
  RooWorkspace* workspaceMCFit1 = (RooWorkspace*)fwMCFit.Get(("workspaceMCFit"+cat1).c_str());
  RooWorkspace* workspaceMCFit2 = (RooWorkspace*)fwMCFit.Get(("workspaceMCFit"+cat2).c_str());


   bool fixparscomb = true;

   //***************Get variables
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooArgSet *observables = (RooArgSet*) workspace->set("observables");
   RooArgSet *combparset = (RooArgSet*) workspace->set((combinatorial.get_name()+"_parset").c_str());
   RooArgSet *Kemuparset;
   if (opts.fitSimultaneous) Kemuparset = (RooArgSet*) workspace->set((Kemu.get_name()+"_parset").c_str());

   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");
   
   //***************Get the PDFs from the workspace
   // RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
   // RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
   // RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");

   RooRealVar* sigmaScaleFactor = new RooRealVar(("sigmaScaleFactor"+cat).c_str(), ("sigmaScaleFactor"+cat).c_str(), 1.1931);//1.12, 0.9, 1.5);
   RooRealVar* meanShift = new RooRealVar(("meanShift"+cat).c_str(), ("meanShift"+cat).c_str(), 5.4707);//0, -100, 100);

   if (!opts.massshift)
   {
     sigmaScaleFactor->setVal(1.);
     meanShift->setVal(0.);
   }
   

   RooAbsPdf *histPdfSignalZeroGamma = build_component_ana("histPdfSignalZeroGamma", *observables,workspaceMCFit0,opts.trigCat,0,
                                                           meanShift, sigmaScaleFactor);
   RooAbsPdf *histPdfSignalOneGamma = build_component_ana("histPdfSignalOneGamma", *observables,workspaceMCFit1,opts.trigCat,1,
                                                          meanShift, sigmaScaleFactor);
   RooAbsPdf *histPdfSignalTwoGamma = build_component_ana("histPdfSignalTwoGamma", *observables,workspaceMCFit2,opts.trigCat,2,
                                                          meanShift, sigmaScaleFactor);
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace->pdf("histPdfPartReco");
   RooHistPdf *histPdfJpsiLeak(0);
   if(opts.nGenJpsiLeak>0) histPdfJpsiLeak = (RooHistPdf *) workspace->pdf("histPdfJpsiLeak");

   RooAbsPdf *combPDF = combinatorial.build_model(observables, combparset, fixparscomb);
   RooAbsPdf *KemuPDF;
   if (opts.fitSimultaneous) KemuPDF = Kemu.build_model(observables, Kemuparset, fixparscomb);
   // expoConst->setVal(trueExp->getVal());

   cout<<"Integral zero: "<<histPdfSignalZeroGamma->createIntegral(RooArgSet(*B_plus_M))->getVal()<<endl;
   cout<<"Integral one: "<<histPdfSignalOneGamma->createIntegral(RooArgSet(*B_plus_M))->getVal()<<endl;
   cout<<"Integral two: "<<histPdfSignalTwoGamma->createIntegral(RooArgSet(*B_plus_M))->getVal()<<endl;
   cout<<"Integral partreco: "<<histPdfPartReco->createIntegral(RooArgSet(*B_plus_M))->getVal()<<endl;
   cout<<"Integral combinatorial: "<<combPDF->createIntegral(RooArgSet(*B_plus_M))->getVal()<<endl;
   

   //**************Prepare fitting function
   double minpartreco = 0.;
   if (opts.fit2D) minpartreco = -50;
   RooRealVar nSignal("nSignal", "#signal events", 0., opts.nGenSignal+100*sqrt(opts.nGenSignal)); //1.*opts.nGenSignal, opts.nGenSignal-10*sqrt(opts.nGenSignal), opts.nGenSignal+10*sqrt(opts.nGenSignal));
   RooRealVar nPartReco("nPartReco", "#nPartReco", minpartreco, opts.nGenPartReco+100*sqrt(opts.nGenPartReco)); //1.*opts.nGenPartReco, opts.nGenPartReco-10*sqrt(opts.nGenPartReco),opts.nGenPartReco+10*sqrt(opts.nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", minpartreco, opts.nGenComb+100*sqrt(opts.nGenComb)); //1.*opts.nGenComb, opts.nGenComb-10*sqrt(opts.nGenComb), opts.nGenComb+10*sqrt(opts.nGenComb));
   RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*opts.nGenJpsiLeak, 0.,//opts.nGenJpsiLeak-10*sqrt(opts.nGenJpsiLeak), 
                        opts.nGenJpsiLeak+100*sqrt(opts.nGenJpsiLeak));
   RooRealVar nKemu("nKemu", "#nKemu", 0., opts.nGenKemu+10*sqrt(opts.nGenKemu)); //1.*opts.nGenKemu, opts.nGenKemu-10*sqrt(opts.nGenKemu),opts.nGenKemu+10*sqrt(opts.nGenKemu));
   RooRealVar fracZero("fracZero", "fracZero",0.5,0,1);
   RooRealVar fracOne("fracOne", "fracOne",0.5, 0,1);
   RooFormulaVar fracPartReco("fracPartReco", "nPartReco/nSignal", RooArgList(nPartReco,nSignal));
   RooFormulaVar fracOneRec("fracOneRec", "(1-fracZero)*fracOne", RooArgList(fracZero, fracOne));


   RooAddPdf histPdfSignal("histPdfSignal", "histPdfSignal", 
                           RooArgList(*histPdfSignalZeroGamma, *histPdfSignalOneGamma, *histPdfSignalTwoGamma), 
                           RooArgList(fracZero, fracOneRec));

   RooArgList pdfList(histPdfSignal, *histPdfPartReco, *combPDF);
   RooArgList yieldList(nSignal, nPartReco, nComb);

   // nJpsiLeak.setConstant(1);

   if(opts.nGenJpsiLeak>0)
   {
      pdfList.add(*histPdfJpsiLeak);
      yieldList.add(nJpsiLeak); 
   }

   RooCategory category("category", "category");
   category.defineType("Kee");
   if (opts.fitSimultaneous) category.defineType("Kemu");

   RooAddPdf *totPdfee = new RooAddPdf("totPdfee", "totPdfee", pdfList, yieldList);   

   RooSimultaneous *totPdf = new RooSimultaneous("totPdf", "totPdf", category);
   totPdf->addPdf(*totPdfee, "Kee");

   RooExtendPdf *totKemuPdf;   
   if  (opts.fitSimultaneous) 
   {
     totKemuPdf = new RooExtendPdf("totKemuPdf", "totKemuPdf", *KemuPDF, nKemu);
     totPdf->addPdf(*totKemuPdf, "Kemu");
   }


   //**************** Constrain the fraction of zero and one photon
   int nGenSignalZeroGamma(floor(opts.nGenFracZeroGamma*opts.nGenSignal));
   int nGenSignalOneGamma(floor(opts.nGenFracOneGamma*opts.nGenSignal));
   int nGenSignalTwoGamma(floor(opts.nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));

   RooRealVar fracZeroMean("fracZeroMean", "fracZeroMean", nGenSignalZeroGamma*1./opts.nGenSignal);
   RooRealVar fracZeroSigma("fracZeroSigma", "fracZeroSigma", sqrt(nGenSignalZeroGamma)/opts.nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroMean, fracZeroSigma); 

   RooRealVar fracOneMean("fracOneMean", "fracOneMean", 
                               nGenSignalOneGamma*1./opts.nGenSignal/(1-fracZeroMean.getVal()));
   RooRealVar fracOneSigma("fracOneSigma", "fracOneSigma", 
                                sqrt(nGenSignalOneGamma)/opts.nGenSignal/(1-fracZeroMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneMean, fracOneSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", opts.nGenPartReco/(1.*opts.nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());
   RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

   RooRealVar nJpsiLeakMean("nJpsiLeakMean", "nJpsiLeakMean", opts.nGenJpsiLeak);
   // RooRealVar nJpsiLeakSigma("nJpsiLeakSigma", "nJpsiLeakSigma", opts.nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
   RooRealVar nJpsiLeakSigma("nJpsiLeakSigma", "nJpsiLeakSigma", 0.1);
   RooGaussian nJpsiLeakConst("nJpsiLeakConst", "nJpsiLeakConst", nJpsiLeak, nJpsiLeakMean, nJpsiLeakSigma); 

   // nJpsiLeak.setConstant(1);

   //Extra TEST CONSTRAINT


   //RooRealVar combConstMean("combConstMean", "combConstMean", opts.nGenComb);
   //RooRealVar combConstSigma("combConstSigma", "combConstSigma", 7.7);
   //RooGaussian combConst("combConst", "combConst", nComb, combConstMean, combConstSigma);

   //**************** fit
   RooDataSet *datasetTot;
   if (opts.fitSimultaneous) datasetTot =  new RooDataSet("datasetTot", "datasetTot", *observables, Index(category),
								    Import("Kee", *dataset_ee), Import("Kemu", *dataset_Kemu));
   else                      datasetTot =  new RooDataSet("datasetTot", "datasetTot", *observables, Index(category),
								   Import("Kee", *dataset_ee));


   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;

   RooArgSet *par_set = totPdf->getParameters(datasetTot);
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


   RooAbsReal* nll = totPdf->createNLL(*datasetTot, Extended(), ExternalConstraints(constraints));
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
      // initiateParams(*trueExp, nSignal, nPartReco, nComb, fracZero, fracOne, *expoConst, nJpsiLeak, opts.constPartReco, fracPartRecoSigma);
     initiateParams(*par_set, constraints, par_set_const);
     cout<<"FITTING: starting with nsignal = "<<nSignal.getValV()<<" refit nbr. "<<i<<endl;
     //if(fitRes != NULL && fitRes != 0) delete fitRes;
     
     migradRes = minuit.migrad();
     hesseRes = minuit.hesse();
     
     fitRes = minuit.save();
     edm = fitRes->edm();
     
     fitResVec.push_back(fitRes); 
     migradResVec.push_back(migradRes);
     hesseResVec.push_back(hesseRes);
     
     if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;
     
     ++nrefit;
     

     cout<<"Fitting nbr "<<i<<" done. Hesse: "
         <<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
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
      cout<<"Fit not converged, choose fit "<<minIndex
          <<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(t, fitRes,  update, migradRes, hesseRes, hasConverged);

   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);
   //totPdf->fitTo(*datasetTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(m_wantplot) 
   {
     // plot_fit_result(*totPdfee, *((RooDataSet*) datasetTot->reduce(Cut("category==category::Kee"))));
     plot_fit_result(*totPdfee, *dataset_ee);
     if (opts.fitSimultaneous) plot_kemu_fit_result(*totKemuPdf, *((RooDataSet*) datasetTot->reduce(Cut("category==category::Kemu"))));
   }

   //delete and return
   delete nll;
   delete par_set;
   delete combPDF;
   delete totPdfee;
   if (opts.fitSimultaneous) {
     delete totKemuPdf;
   }
   delete totPdf;
   
}

void FitterUtilsAna::fit_mumu(RooDataSet* datasetTot, RooWorkspace* workspace, 
                      double fracPartReco_const,
                      ofstream& out, TTree* t, bool update)
{

   //***************Get variables
  RooRealVar *B_plus_M = workspace->var("B_plus_M");
  RooArgSet *observables = (RooArgSet*) workspace->set("observables");
  
   string compname = "Signal_mumu";
   RooRealVar* meanB = (RooRealVar*)workspace->var(("meanB_"+compname).c_str());
   meanB->setConstant(true);
   cout<<"Hola1"<<endl;
   RooRealVar* sigma = (RooRealVar*)workspace->var(("sigma_"+compname).c_str());
   sigma->setConstant(true);
   cout<<"Hola2"<<endl;
   RooRealVar* sigma1 = (RooRealVar*)workspace->var(("sigma1_"+compname).c_str());
   sigma1->setConstant(true);
   // RooRealVar* sigma2 = (RooRealVar*)workspace->var(("sigma2_"+compname).c_str());
   // sigma2->setConstant(true);
   cout<<"Hola3"<<endl;
   RooRealVar* n = (RooRealVar*)workspace->var(("n_"+compname).c_str());
   n->setConstant(true);
   cout<<"Hola4"<<endl;
   RooRealVar* alpha = (RooRealVar*)workspace->var(("alpha_"+compname).c_str());
   alpha->setConstant(true);
   // RooRealVar* n1 = (RooRealVar*)workspace->var(("n1_"+compname).c_str());
   // n1->setConstant(true);
   // RooRealVar* alpha1 = (RooRealVar*)workspace->var(("alpha1_"+compname).c_str());
   // alpha1->setConstant(true);
   // RooRealVar* n2 = (RooRealVar*)workspace->var(("n2_"+compname).c_str());
   // n2->setConstant(true);
   // RooRealVar* alpha2 = (RooRealVar*)workspace->var(("alpha2_"+compname).c_str());
   // alpha2->setConstant(true);
   // RooRealVar* n1b = (RooRealVar*)workspace->var(("n1b_"+compname).c_str());
   // n1b->setConstant(true);
   // RooRealVar* alpha1b = (RooRealVar*)workspace->var(("alpha1b_"+compname).c_str());
   // alpha1b->setConstant(true);
   // RooRealVar* n2b = (RooRealVar*)workspace->var(("n2b_"+compname).c_str());
   // n2b->setConstant(true);
   // RooRealVar* alpha2b = (RooRealVar*)workspace->var(("alpha2b_"+compname).c_str());
   // alpha2b->setConstant(true);
   cout<<"Hola5"<<endl;
   RooRealVar* vc1 = (RooRealVar*)workspace->var(("vc1_"+compname).c_str());
   vc1->setConstant(true);
   cout<<"Hola6"<<endl;
   RooRealVar *lambda = workspace->var(("lambda_"+compname).c_str()); 



   //*****************Build pdfs
   RooCBShape cb1(("cb1_"+compname).c_str(), "cb1", *B_plus_M,  *meanB,  *sigma,  *alpha, *n);
   RooCBShape cb2(("cb2_"+compname).c_str(), "cb2", *B_plus_M,  *meanB,  *sigma1,  *alpha, *n);

   RooAddPdf model_signal("Signal_mumu","Signal_mumu",RooArgList( cb1, cb2),RooArgList( *vc1)) ;
   RooExponential model_bkg("Bkg_mumu", "Bkg_mumu", *B_plus_M,  *lambda);

   //**************Prepare fitting function
   RooRealVar nSignal("nSignal_mumu", "#signal events mumu", 1.*opts.nGenSignal_mumu, opts.nGenSignal_mumu-7*sqrt(opts.nGenSignal_mumu), 
                      opts.nGenSignal_mumu+7*sqrt(opts.nGenSignal_mumu));
   RooRealVar nComb("nComb_mumu", "#nComb mumu", 1.*opts.nGenBkg_mumu, opts.nGenBkg_mumu-7*sqrt(opts.nGenBkg_mumu), 
                    opts.nGenBkg_mumu+7*sqrt(opts.nGenBkg_mumu));

   RooArgList pdfList(model_signal, model_bkg);
   RooArgList yieldList(nSignal, nComb);

   RooAddPdf *totPdf = new RooAddPdf("totPdfmumu", "totPdfmumu", pdfList, yieldList);   


   //**************** Constrain the fraction of zero and one photon

   //Extra TEST CONSTRAINT


   //RooRealVar combConstMean("combConstMean", "combConstMean", opts.nGenComb);
   //RooRealVar combConstSigma("combConstSigma", "combConstSigma", 7.7);
   //RooGaussian combConst("combConst", "combConst", nComb, combConstMean, combConstSigma);

   //**************** fit
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;

   RooArgSet *par_set = totPdf->getParameters(datasetTot);
   initiateParams(*par_set);
   // initiateParams(*trueExp, nSignal, nPartReco, nComb, fracZero, fracOne, *expoConst, nJpsiLeak,  opts.constPartReco, fracPartRecoSigma);



   RooAbsReal* nll = totPdf->createNLL(*datasetTot, Extended());//, ExternalConstraints(constraints));
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
      // initiateParams(*trueExp, nSignal, nPartReco, nComb, fracZero, fracOne, *expoConst, nJpsiLeak, opts.constPartReco, fracPartRecoSigma);
     initiateParams(*par_set);
     cout<<"FITTING: starting with nsignal = "<<nSignal.getValV()<<" refit nbr. "<<i<<endl;
     //if(fitRes != NULL && fitRes != 0) delete fitRes;
     
     migradRes = minuit.migrad();
     hesseRes = minuit.hesse();
     
     fitRes = minuit.save();
     edm = fitRes->edm();
     
     fitResVec.push_back(fitRes); 
     migradResVec.push_back(migradRes);
     hesseResVec.push_back(hesseRes);
     
     if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;
     
     ++nrefit;
     

     cout<<"Fitting nbr "<<i<<" done. Hesse: "
         <<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
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
      cout<<"Fit not converged, choose fit "<<minIndex
          <<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(t, fitRes,  update, migradRes, hesseRes, hasConverged);

   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);
   //totPdf->fitTo(*datasetTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(m_wantplot) 
   {
     plot_mumu_fit_result(*totPdf, *datasetTot);
   }

   //delete and return
   delete nll;
   delete par_set;
   delete totPdf;
   
}



void FitterUtilsAna::fit_RK(RooDataSet* dataset_ee, RooDataSet* dataset_mumu, RooDataSet* dataset_Kemu,  
			 RooWorkspace* workspace_ee, RooWorkspace* workspace_mumu, 
			 double fracPartReco_const,
			 ofstream& out, TTree* t, bool update)
{


   bool fixparscomb = true;

   //***************Get variables ee
   RooRealVar *B_plus_M = workspace_ee->var("B_plus_M");
   RooArgSet *observables = (RooArgSet*) workspace_ee->set("observables");
   RooArgSet *combparset = (RooArgSet*) workspace_ee->set((combinatorial.get_name()+"_parset").c_str());
   RooArgSet *Kemuparset;
   if (opts.fitSimultaneous) Kemuparset = (RooArgSet*) workspace_ee->set((Kemu.get_name()+"_parset").c_str());

   RooRealVar *fractionalErrorJpsiLeak = workspace_ee->var("fractionalErrorJpsiLeak");
   

   //***************Get the PDFs from the workspace_ee
   RooHistPdf *histPdfSignalZeroGamma = (RooHistPdf *) workspace_ee->pdf("histPdfSignalZeroGamma");
   RooHistPdf *histPdfSignalOneGamma = (RooHistPdf *) workspace_ee->pdf("histPdfSignalOneGamma");
   RooHistPdf *histPdfSignalTwoGamma = (RooHistPdf *) workspace_ee->pdf("histPdfSignalTwoGamma");
   RooHistPdf *histPdfPartReco = (RooHistPdf *) workspace_ee->pdf("histPdfPartReco");
   RooHistPdf *histPdfJpsiLeak(0);
   if(opts.nGenJpsiLeak>0) histPdfJpsiLeak = (RooHistPdf *) workspace_ee->pdf("histPdfJpsiLeak");

   RooAbsPdf *combPDF = combinatorial.build_model(observables, combparset, fixparscomb);
   RooAbsPdf *KemuPDF;
   if (opts.fitSimultaneous) KemuPDF = Kemu.build_model(observables, Kemuparset, fixparscomb);
   // expoConst->setVal(trueExp->getVal());



   //***************Get the PDFs from the workspace_mumu
   string compname = "Signal_mumu";
   RooRealVar* meanB = (RooRealVar*)workspace_mumu->var(("meanB_"+compname).c_str());
   meanB->setConstant(true);
   RooRealVar* sigma = (RooRealVar*)workspace_mumu->var(("sigma_"+compname).c_str());
   sigma->setConstant(true);
   RooRealVar* sigma1 = (RooRealVar*)workspace_mumu->var(("sigma1_"+compname).c_str());
   sigma1->setConstant(true);
   // RooRealVar* sigma2 = (RooRealVar*)workspace_mumu->var(("sigma2_"+compname).c_str());
   // sigma2->setConstant(true);
   RooRealVar* n = (RooRealVar*)workspace_mumu->var(("n_"+compname).c_str());
   n->setConstant(true);
   RooRealVar* alpha = (RooRealVar*)workspace_mumu->var(("alpha_"+compname).c_str());
   alpha->setConstant(true);
   // RooRealVar* n1 = (RooRealVar*)workspace_mumu->var(("n1_"+compname).c_str());
   // n1->setConstant(true);
   // RooRealVar* alpha1 = (RooRealVar*)workspace_mumu->var(("alpha1_"+compname).c_str());
   // alpha1->setConstant(true);
   // RooRealVar* n2 = (RooRealVar*)workspace_mumu->var(("n2_"+compname).c_str());
   // n2->setConstant(true);
   // RooRealVar* alpha2 = (RooRealVar*)workspace_mumu->var(("alpha2_"+compname).c_str());
   // alpha2->setConstant(true);
   // RooRealVar* n1b = (RooRealVar*)workspace_mumu->var(("n1b_"+compname).c_str());
   // n1b->setConstant(true);
   // RooRealVar* alpha1b = (RooRealVar*)workspace_mumu->var(("alpha1b_"+compname).c_str());
   // alpha1b->setConstant(true);
   // RooRealVar* n2b = (RooRealVar*)workspace_mumu->var(("n2b_"+compname).c_str());
   // n2b->setConstant(true);
   // RooRealVar* alpha2b = (RooRealVar*)workspace_mumu->var(("alpha2b_"+compname).c_str());
   // alpha2b->setConstant(true);
   RooRealVar* vc1 = (RooRealVar*)workspace_mumu->var(("vc1_"+compname).c_str());
   vc1->setConstant(true);
   RooRealVar *lambda = workspace_mumu->var(("lambda_"+compname).c_str()); 

   //*****************Build pdfs
   RooCBShape cb1(("cb1_"+compname).c_str(), "cb1", *B_plus_M,  *meanB,  *sigma,  *alpha, *n);
   RooCBShape cb2(("cb2_"+compname).c_str(), "cb2", *B_plus_M,  *meanB,  *sigma1,  *alpha, *n);

   RooAddPdf model_mumu_signal("model_mumu_signal","model_mumu_signal",RooArgList( cb1, cb2),RooArgList( *vc1)) ;
   RooExponential model_mumu_bkg("model_mumu_bkg", "model_mumu_bkg", *B_plus_M,  *lambda);



   //**************Prepare fitting function
   RooRealVar nSignal_mumu("nSignal_mumu", "#signal events", 1.*opts.nGenSignal_mumu, opts.nGenSignal_mumu-7*sqrt(opts.nGenSignal_mumu), 
                      opts.nGenSignal_mumu+7*sqrt(opts.nGenSignal_mumu));

   RooRealVar nBkg_mumu("nBkg_mumu", "#signal events", 1.*opts.nGenBkg_mumu, opts.nGenBkg_mumu-7*sqrt(opts.nGenBkg_mumu), 
                      opts.nGenBkg_mumu+7*sqrt(opts.nGenBkg_mumu));

   RooRealVar RK("RK", "#signal events", 1. ,0, 2);

   RooFormulaVar nSignal("nSignal","nSignal_mumu/RK",RooArgList(nSignal_mumu,RK));

   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*opts.nGenPartReco, opts.nGenPartReco-7*sqrt(opts.nGenPartReco), 
                        opts.nGenPartReco+7*sqrt(opts.nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*opts.nGenComb, opts.nGenComb-7*sqrt(opts.nGenComb), 
                    opts.nGenComb+7*sqrt(opts.nGenComb));
   RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*opts.nGenJpsiLeak, opts.nGenJpsiLeak-7*sqrt(opts.nGenJpsiLeak), 
                        opts.nGenJpsiLeak+7*sqrt(opts.nGenJpsiLeak));
   RooRealVar nKemu("nKemu", "#nKemu", 1.*opts.nGenKemu, opts.nGenKemu-7*sqrt(opts.nGenKemu),opts.nGenKemu+7*sqrt(opts.nGenKemu));


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

   // Model mumu
   RooAddPdf *totPdfmumu = new RooAddPdf("totPdfmumu", "totPdfmumu", RooArgList(model_mumu_signal, model_mumu_bkg), RooArgList(nSignal_mumu, nBkg_mumu));   

   // Model ee
   RooAddPdf *totPdfee = new RooAddPdf("totPdfee", "totPdfee", pdfList, yieldList);   

   
   RooCategory category("category", "category");
   category.defineType("Kmumu");
   category.defineType("Kemu");
   if (opts.fitSimultaneous) category.defineType("Kee");

   RooSimultaneous *totPdf = new RooSimultaneous("totPdf", "totPdf", category);     
   totPdf->addPdf(*totPdfee, "Kee");
   totPdf->addPdf(*totPdfmumu, "Kmumu");

   RooExtendPdf *totKemuPdf;
   if  (opts.fitSimultaneous) 
   {
     totKemuPdf = new RooExtendPdf("totKemuPdf", "totKemuPdf", *KemuPDF, nKemu);
     totPdf->addPdf(*totKemuPdf, "Kemu");
   }



   //**************** Constrain the fraction of zero and one photon
   int nGenSignalZeroGamma(floor(opts.nGenFracZeroGamma*opts.nGenSignal));
   int nGenSignalOneGamma(floor(opts.nGenFracOneGamma*opts.nGenSignal));
   int nGenSignalTwoGamma(floor(opts.nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));

   RooRealVar fracZeroMean("fracZeroMean", "fracZeroMean", nGenSignalZeroGamma*1./opts.nGenSignal);
   RooRealVar fracZeroSigma("fracZeroSigma", "fracZeroSigma", sqrt(nGenSignalZeroGamma)/opts.nGenSignal);
   RooGaussian fracZeroConst("fracZeroConst", "fracZeroConst", fracZero, fracZeroMean, fracZeroSigma); 

   RooRealVar fracOneMean("fracOneMean", "fracOneMean", 
                               nGenSignalOneGamma*1./opts.nGenSignal/(1-fracZeroMean.getVal()));
   RooRealVar fracOneSigma("fracOneSigma", "fracOneSigma", 
                                sqrt(nGenSignalOneGamma)/opts.nGenSignal/(1-fracZeroMean.getVal()));
   RooGaussian fracOneConst("fracOneConst", "fracOneConst", fracOne, fracOneMean, fracOneSigma); 

   RooRealVar fracPartRecoMean("fracPartRecoMean", "fracPartRecoMean", opts.nGenPartReco/(1.*opts.nGenSignal));
   RooRealVar fracPartRecoSigma("fracPartRecoSigma", "fracPartRecoSigma", fracPartReco_const*fracPartRecoMean.getVal());
   RooGaussian fracPartRecoConst("fracPartRecoConst", "fracPartRecoConst", fracPartReco, fracPartRecoMean, fracPartRecoSigma);

   RooRealVar nJpsiLeakMean("nJpsiLeakMean", "nJpsiLeakMean", opts.nGenJpsiLeak);
   RooRealVar nJpsiLeakSigma("nJpsiLeakSigma", "nJpsiLeakSigma", opts.nGenJpsiLeak*fractionalErrorJpsiLeak->getVal());
   RooGaussian nJpsiLeakConst("nJpsiLeakConst", "nJpsiLeakConst", nJpsiLeak, nJpsiLeakMean, nJpsiLeakSigma); 

   //Extra TEST CONSTRAINT


   //RooRealVar combConstMean("combConstMean", "combConstMean", opts.nGenComb);
   //RooRealVar combConstSigma("combConstSigma", "combConstSigma", 7.7);
   //RooGaussian combConst("combConst", "combConst", nComb, combConstMean, combConstSigma);

   //**************** fit
   RooDataSet *datasetTot;
   if (!opts.fitSimultaneous) datasetTot = new RooDataSet("datasetTot", "datasetTot", *observables, Index(category),
								   Import("Kee", *dataset_ee), Import("Kmumu", *dataset_mumu));
   else                       datasetTot = new RooDataSet("datasetTot", "datasetTot", *observables, Index(category),
								   Import("Kee", *dataset_ee), Import("Kmumu", *dataset_mumu), Import("Kemu", *dataset_Kemu));
    
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;


   RooArgSet *par_set = totPdf->getParameters(datasetTot);
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


   B_plus_M->setRange("fit_Kee",opts.minBMass,opts.maxBMass);
   B_plus_M->setRange("fit_Kmumu",opts.minBMass_mumu,opts.maxBMass_mumu);
   if (opts.fitSimultaneous) B_plus_M->setRange("fit_Kemu",opts.minBMass,opts.maxBMass);
   
   totPdfmumu->fixAddCoefRange("fit_Kmumu") ;

   RooAbsReal* nll = totPdf->createNLL(*datasetTot, Extended(1), Range("fit"), SplitRange(), NormRange("fit"), ExternalConstraints(constraints));
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
      // initiateParams(*trueExp, nSignal, nPartReco, nComb, fracZero, fracOne, *expoConst, nJpsiLeak, opts.constPartReco, fracPartRecoSigma);
     initiateParams(*par_set,constraints, par_set_const);
     cout<<"FITTING: starting with nsignal = "<<nSignal.getValV()<<" refit nbr. "<<i<<endl;
     //if(fitRes != NULL && fitRes != 0) delete fitRes;
     
     migradRes = minuit.migrad();
     hesseRes = minuit.hesse();
     
     fitRes = minuit.save();
     edm = fitRes->edm();
     
     fitResVec.push_back(fitRes); 
     migradResVec.push_back(migradRes);
     hesseResVec.push_back(hesseRes);
     
     if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;
     
     ++nrefit;
     

     cout<<"Fitting nbr "<<i<<" done. Hesse: "
         <<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
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
      cout<<"Fit not converged, choose fit "<<minIndex
          <<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(t, fitRes,  update, migradRes, hesseRes, hasConverged);

   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);
   //totPdf->fitTo(*datasetTot, Extended(), Save(), Warnings(false));

   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(m_wantplot) 
   {

     plot_fit_result(*totPdfee, *((RooDataSet*) datasetTot->reduce(Cut("category==category::Kee"))));
     plot_mumu_fit_result(*totPdfmumu, *((RooDataSet*) datasetTot->reduce(Cut("category==category::Kmumu"))));

     if (opts.fitSimultaneous)   plot_kemu_fit_result(*totKemuPdf, *((RooDataSet*) datasetTot->reduce(Cut("category==category::Kemu"))));
   }

   //delete and return
   delete nll;
   delete par_set;
   delete combPDF;
   delete totPdfee;
   if (opts.fitSimultaneous) {
     delete totKemuPdf;
   }
   delete totPdf;
   
}





void FitterUtilsAna::plot_fit_result(RooAbsPdf &totPdf, RooDataSet dataGenTot)
{

   //**************Prepare TFile to save the plots

   TFile f2(opts.plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   while((var = (RooRealVar*) iter->Next()))
   {

      frame = var->frame();
      dataGenTot.plotOn(frame, Binning(defaultbinningname.c_str()));
      // dataGenTot.plotOn(frame);
      totPdf.plotOn(frame, Components("histPdfPartReco"), LineColor(kBlue));
      totPdf.plotOn(frame, Components("histPdfSignalZeroGamma"), LineColor(kGreen));
      totPdf.plotOn(frame, Components("histPdfSignalOneGamma"), LineColor(kMagenta));
      totPdf.plotOn(frame, Components("histPdfSignalTwoGamma"), LineColor(kOrange));
      totPdf.plotOn(frame, Components("histPdfJpsiLeak"), LineColor(14));
      totPdf.plotOn(frame, Components(("combPDF_"+combinatorial.get_name()).c_str()), LineColor(kBlack));
      totPdf.plotOn(frame, LineColor(kRed));

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

void FitterUtilsAna::plot_mumu_fit_result(RooAbsPdf &totPdf, RooDataSet dataGenTot)
{

   //**************Prepare TFile to save the plots

   TFile f2(opts.plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   while((var = (RooRealVar*) iter->Next()))
   {

     frame = var->frame(Range(opts.minBMass_mumu,opts.maxBMass_mumu));
     dataGenTot.plotOn(frame);
     
     totPdf.fixAddCoefRange("fit_Kmumu") ;

     totPdf.plotOn(frame, Range("fit_Kmumu"), NormRange("fit_Kmumu"), Components("model_mumu_signal"), LineColor(kBlue));
     totPdf.plotOn(frame, Range("fit_Kmumu"), NormRange("fit_Kmumu"), Components("model_mumu_bkg"), LineColor(kBlack));
     totPdf.plotOn(frame, Range("fit_Kmumu"), NormRange("fit_Kmumu"), LineColor(kRed));

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



void FitterUtilsAna::plot_kemu_fit_result(RooAbsPdf &totKemuPdf, RooDataSet datasetKemu)
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
    datasetKemu.plotOn(frame, Binning(defaultbinningname.c_str()));
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

void FitterUtilsAna::PlotShape(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, 
                            string canvName, RooArgSet varset)
{


   RooRealVar *B_plus_M = (RooRealVar*) varset.find("B_plus_M");
   RooRealVar *misPT = (RooRealVar*) varset.find("misPT");

  if(opts.fit2D) PlotShape2D(originDataSet, genDataSet, shape, plotsfile, canvName,*B_plus_M,*misPT);
  if(!opts.fit2D) PlotShape1D(originDataSet, genDataSet, shape, plotsfile, canvName,*B_plus_M);
}


void FitterUtilsAna::PlotShape2D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, 
                              string canvName, RooRealVar& B_plus_M, RooRealVar& misPT)
{
   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");

   //**************Plot Signal Zero Gamma

   TH2F* th2fKey = (TH2F*)shape.createHistogram("th2Shape", B_plus_M, Binning(20), YVar(misPT, Binning(20)));
   cout<<genDataSet.sumEntries()<<endl;
   TH2F* th2fGen = (TH2F*)genDataSet.createHistogram("th2fGen", B_plus_M, Binning(20), YVar(misPT, Binning(20)));

   RooPlot* plotM = B_plus_M.frame();
   originDataSet.plotOn(plotM,  Binning(defaultbinningname.c_str()));
   shape.plotOn(plotM);

   RooPlot* plotMisPT = misPT.frame();
   originDataSet.plotOn(plotMisPT,  Binning(defaultbinningname.c_str()));
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

void FitterUtilsAna::PlotShape1D(RooDataSet& originDataSet, RooDataSet& genDataSet, RooAbsPdf& shape, string plotsfile, 
                              string canvName, RooRealVar& B_plus_M)
{
   TFile f2(plotsfile.c_str(), "UPDATE");

   RooPlot* plotGen = B_plus_M.frame();
   genDataSet.plotOn(plotGen,  Binning(defaultbinningname.c_str()));

   RooPlot* plotM = B_plus_M.frame();
   originDataSet.plotOn(plotM, Binning(defaultbinningname.c_str()));
   shape.plotOn(plotM);

   TCanvas canv(canvName.c_str(), canvName.c_str(), 800, 800);
   canv.Divide(1,2);
   canv.cd(1); plotGen->Draw();
   canv.cd(2); plotM->Draw();

   canv.Write();

   f2.Close();
}




void FitterUtilsAna::initiateParams(RooRealVar const& expoConstGen, RooRealVar& nSignal, RooRealVar& nPartReco, 
                                 RooRealVar& nComb, RooRealVar& fracZero, RooRealVar& fracOne, 
                                 RooRealVar& expoConst, RooRealVar&  nJpsiLeak, 
                                 bool constPartReco, RooRealVar const& fracPartRecoSigma)
{
   TRandom rand;
   rand.SetSeed();

   int nGenSignalZeroGamma(floor(opts.nGenFracZeroGamma*opts.nGenSignal));
   int nGenSignalOneGamma(floor(opts.nGenFracOneGamma*opts.nGenSignal));
   int nGenSignalTwoGamma(floor(opts.nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma));

   int nGenSignal = nGenSignalZeroGamma + nGenSignalOneGamma + nGenSignalTwoGamma;

   double nGenSignal2;
   double nGenPartReco2;
   if(!constPartReco)
   {
      nGenSignal2 = rand.Uniform(opts.nGenSignal-5*sqrt(opts.nGenSignal), opts.nGenSignal+5*sqrt(opts.nGenSignal));
      nGenPartReco2 = rand.Uniform(opts.nGenPartReco-5*sqrt(opts.nGenPartReco), opts.nGenPartReco+5*sqrt(opts.nGenPartReco));
   }
   if(constPartReco)
   { 
      double nGenSigPartReco( opts.nGenSignal+opts.nGenPartReco );
      double nGenSigPartReco2( rand.Uniform( nGenSigPartReco-5*sqrt(nGenSigPartReco), nGenSigPartReco+5*sqrt(nGenSigPartReco) ) );
      double fracPartReco1( opts.nGenPartReco/(1.*opts.nGenSignal));
      double fracPartReco2( rand.Uniform(fracPartReco1-5*fracPartRecoSigma.getVal(), fracPartReco1+5*fracPartRecoSigma.getVal())); 

      nGenPartReco2 = fracPartReco2*nGenSigPartReco2 / (1+fracPartReco2); 
      nGenSignal2 = nGenSigPartReco2 / (1+fracPartReco2); 
   }
   double nGenComb2 = rand.Uniform(opts.nGenComb-5*sqrt(opts.nGenComb), opts.nGenComb+5*sqrt(opts.nGenComb));
   double nGenJpsiLeak2 = rand.Uniform(opts.nGenJpsiLeak-5*sqrt(opts.nGenJpsiLeak), opts.nGenJpsiLeak+5*sqrt(opts.nGenJpsiLeak));


   nSignal.setVal(nGenSignal2);
   nSignal.setRange(TMath::Max(0.,nGenSignal2-10.*sqrt(opts.nGenSignal)) , nGenSignal2+10*sqrt(opts.nGenSignal));

   nPartReco.setVal(nGenPartReco2);
   nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(opts.nGenPartReco)), nGenPartReco2+10*sqrt(opts.nGenPartReco));


   nComb.setVal(nGenComb2);
   nComb.setRange(TMath::Max(0.,nGenComb2-10.*sqrt(opts.nGenComb)), nGenComb2+10*sqrt(opts.nGenComb));

   nJpsiLeak.setVal(nGenJpsiLeak2);
   nJpsiLeak.setRange(TMath::Max(0., nGenJpsiLeak2-10*sqrt(opts.nGenJpsiLeak)), nGenJpsiLeak2+10*sqrt(opts.nGenJpsiLeak));

   double fracGenZero(nGenSignalZeroGamma/(1.*opts.nGenSignal));
   double fracGenOne(nGenSignalOneGamma/(1.*opts.nGenSignal));

   fracZero.setVal(rand.Gaus(fracGenZero, sqrt(nGenSignalZeroGamma)/(1.*opts.nGenSignal))) ;
   fracZero.setRange(0., 1.);
   fracOne.setVal(rand.Gaus(fracGenOne, sqrt(nGenSignalOneGamma)/(1.*opts.nGenSignal))) ;
   fracOne.setRange(0., 1.);

   expoConst.setVal(rand.Uniform( expoConstGen.getVal() - 5*expoConstGen.getError(), 
                                  expoConstGen.getVal() + 5*expoConstGen.getError() ) );
   expoConst.setRange( expoConstGen.getVal() - 10*expoConstGen.getError(), expoConstGen.getVal() + 10*expoConstGen.getError() );


   // if (opts.fitSimultaneous)
   // {
   //   nKemu.setVal(rand.Uniform(opts.nGenKemu-5*sqrt(opts.nGenKemu), opts.nGenKemu+5*sqrt(opts.nGenKemu)));
   //   nKemu.setRange(opts.nGenKemu-10*sqrt(opts.nGenKemu), opts.nGenKemu+10*sqrt(opts.nGenKemu));

   //   T.setVal(rand.Uniform( TMath::Max(0.,TGen.getVal() - 5*TGen.getError()), TGen.getVal() + 5*TGen.getError()));
   //   T.setRange(TMath::Max(TGen.getVal() - 10*TGen.getError(),0.), TGen.getVal() + 10*TGen.getError());

   //   n.setVal(rand.Uniform(TMath::Max( nGen.getVal() - 5*nGen.getError(),2.), nGen.getVal() + 5*nGen.getError()));
   //   n.setRange(TMath::Max(nGen.getVal() - 10*nGen.getError(),2.), nGen.getVal() + 10*nGen.getError());

   //   expoConstKemu.setVal(rand.Uniform( expoConstGen.getVal() - 5*expoConstGen.getError(),
   //                                      expoConstGen.getVal() + 5*expoConstGen.getError() ) );
   //   expoConstKemu.setRange( expoConstGen.getVal() - 10*expoConstGen.getError(),
   //                           expoConstGen.getVal() + 10*expoConstGen.getError() );
   // }

}
