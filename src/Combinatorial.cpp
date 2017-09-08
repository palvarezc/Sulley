// Include files 
#include "Combinatorial.h"
#include "RooChebychev.h"

Combinatorial::Combinatorial():
  name(""),
  fittype(FitType::RooFit1D),
  ownexp(0)
{
}


Combinatorial::Combinatorial( string name_, FitType fittype_, string plotsfile_ ):
  name(name_),
  fittype(fittype_),
  plotsfile(plotsfile_),
  ownexp(0)
{
}


void Combinatorial::create_parameters(RooWorkspace* workspace)
{


  RooArgSet parameters((name+"_parset").c_str()); 

  if (fittype==RooFit1D || fittype==HistFact1D)
  {
    RooRealVar expoConst("expoConst", "expoConst", -1e-3, -0.006, 0.006);
    parameters.add(expoConst);
    workspace->defineSet((name+"_parset").c_str(), parameters, true);

    // RooRealVar expoConst("expoConst", "expoConst", -1.,1.);
    // parameters.add(expoConst);
    // workspace->defineSet((name+"_parset").c_str(), parameters, true);
  }
  
  if (fittype==RooFit2D_RooPTMVis)
  {
    RooRealVar expoConst("expoConst", "expoConst", -1e-3, -0.006, 0.006);
    RooRealVar T("T", "T", 97, 0, 200);
    RooRealVar n("n", "n", 3.5, 1., 5.5);
    parameters.add(expoConst);
    parameters.add(T);
    parameters.add(n);
    workspace->defineSet((name+"_parset").c_str(), parameters, true);
}

  if (fittype==RooFit2D || fittype==HistFact2D)
  { 
    RooRealVar l1Kee("l1Kee", "l1Kee", +2.3e-3, -1e-1, 1e-1);    
    RooRealVar l2Kee("l2Kee", "l2Kee", +4.9e-3, -1e-1, 1e-1);
    RooRealVar l3Kee("l3Kee", "l3Kee", -4.1e-7, -1e-5, 1e-5);
    RooRealVar l4Kee("l4Kee", "l4Kee", -4e-7, -1e-5, 1e-5);
    RooRealVar l5Kee("l5Kee", "l5Kee", 1e-10, -1e-9, 1e-9);
    parameters.add(l1Kee);
    parameters.add(l2Kee);
    parameters.add(l3Kee);
    parameters.add(l4Kee);
    parameters.add(l5Kee);
    workspace->defineSet((name+"_parset").c_str(), parameters, true);
    cout<<"parameters created..."<<endl;
  }

  // workspace->Print();

}


void Combinatorial::copy_parameters(Combinatorial *other, RooWorkspace* workspace, bool _ownexp)
{

  use_own_exp(_ownexp);
  
  RooArgSet* parameters  = (RooArgSet*) workspace->set((other->get_name()+"_parset").c_str());
  RooRealVar expoConst(("expoConst_"+name).c_str(), ("expoConst_"+name).c_str(), -1e-3, -0.006, 0.006);
  // RooRealVar expoConst(("expoConst_"+name).c_str(), ("expoConst_"+name).c_str(), -1.,1.);

  if (ownexp) parameters->add(expoConst);

  workspace->defineSet((name+"_parset").c_str(), *parameters, true);

  // workspace->import(expoConst);
  // workspace->extendSet((name+"_parset").c_str(),("expoConst_"+name).c_str());

}



// void Combinatorial::delete_set(RooArgSet *parameters)
// {

//   RooRealVar *var;  
//   TIterator *iter = parameters->createIterator();

//   while((var = (RooRealVar*) iter->Next()))
//   {
//     delete var;
//   }
  
// }


void Combinatorial::build_model_from_data(string combfile, string combtree, string combcuts, 
                                          RooArgSet *varset,RooWorkspace *workspace,
                                          bool wantplot, string binningname)
{

  cout<<"Hello!"<<endl;

  create_parameters(workspace);
  RooArgSet* parameters = (RooArgSet*) workspace->set((name+"_parset").c_str());
  RooAbsPdf* model = build_model(varset, parameters, false, binningname); 
  
  TFile* fComb = new TFile(combfile.c_str());  
  TTree* tComb = (TTree*)fComb->Get(combtree.c_str());

  tComb->SetBranchStatus("*",0);
  

  TIterator *variter = varset->createIterator();
  RooRealVar *var = (RooRealVar*) variter->Next();
  while (var)
  {
    tComb->SetBranchStatus(var->GetName(), 1);
    var = (RooRealVar*) variter->Next();
  }
  

  RooDataSet dataSetComb("dataSetComb", "dataSetComb", tComb, *varset, combcuts.c_str());
  model->fitTo(dataSetComb);

  if(wantplot)
    {
      plot(dataSetComb, *model, "calibration");
    }

  workspace->saveSnapshot("pars_for_generation",*parameters,kTRUE) ;

  fComb->Close();
  delete model;
}



void Combinatorial::generate(string compname, RooArgSet *observables, 
                             int nGenComp, RooWorkspace* workspace, 
			     RooWorkspace* workspaceGen, bool wantplot, 
			     string binningname)
{


  RooArgSet *parset = (RooArgSet*) workspace->set((name+"_parset").c_str());
  workspace->loadSnapshot("pars_for_generationname");
  RooAbsPdf *model = build_model(observables, parset, false, binningname);

  RooAbsPdf::GenSpec* GenSpecComp = model->prepareMultiGen(*observables, RooFit::Extended(1), NumEvents(nGenComp));

  cout<<"Generating "<<compname<<endl;
  RooDataSet* dataGenComp = model->generate(*GenSpecComp);
  dataGenComp->SetName(("dataGen"+compname).c_str()); dataGenComp->SetTitle(("dataGen"+compname).c_str());
  
  // cout<<"Number of events requested: "<<nGenComp<<endl;
  // cout<<"Number of events generated: "<<dataGenComp->sumEntries()<<endl;

  if (wantplot)
  {
    plot(*dataGenComp, *model, "generated");
  }   

  workspaceGen->import(*dataGenComp);
  
  delete dataGenComp;
  delete GenSpecComp;
  delete model;
}


RooAbsPdf* Combinatorial::build_model(RooArgSet *varset, RooArgSet *parset, bool fixpars, string binningname)
{

  RooAbsPdf *model = NULL;
  string expoConstname = "expoConst";
  if (ownexp) expoConstname += "_"+name;

  // Single exponential
  if (fittype==RooFit1D || fittype==HistFact1D)
  {
    //Get observables
    RooRealVar * B_plus_M = (RooRealVar*) varset->find("B_plus_M");

    //Get parameters 
    RooRealVar * expoConst = (RooRealVar*) parset->find(expoConstname.c_str());

    if (!B_plus_M || !expoConst)
    {
      cout<<"Varset does not contain variables needed to build the model."<<endl;
      return model;
    }

    if (fittype==RooFit1D)
    {
      cout<<"Going for fittype: "<< "RooFit1D" <<endl;
      model =  new RooExponential(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *expoConst);
      // model =  new RooChebychev(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *expoConst);
    }
  
    if (fittype==HistFact1D)
    {
      cout<<"Going for fittype: "<< "HistFact1D" <<endl;
      model =  new RooExpBinned(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *expoConst, 
                                B_plus_M->getBinning(binningname.c_str()));
    }

    
  }
  
  
  // 2D - RooPTMVis
  if (fittype==RooFit2D_RooPTMVis)
  {

    //Get observables
    RooRealVar * B_plus_M = (RooRealVar*) varset->find("B_plus_M");
    RooRealVar * misPT = (RooRealVar*) varset->find("misPT");

    //Get parameters 
    RooRealVar * expoConst = (RooRealVar*) parset->find(expoConstname.c_str());
    RooRealVar * n = (RooRealVar*) parset->find("n");
    RooRealVar * T = (RooRealVar*) parset->find("T");

    if (!B_plus_M || !misPT || !expoConst || !n || !T)
    {
      cout<<"Varset does not contain variables needed to build the model."<<endl;
      return model;
    }

    if(fixpars)
      {
	n->setConstant(1);
	T->setConstant(1);
      }

    cout<<"Going for fittype: "<< "RooFit2D_RooPTMVis" <<endl;
    model =  new RooPTMVis(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *misPT, *B_plus_M, *T, *n, *expoConst);


  }



  // 2D - RooPolyTimesX
  if (fittype==RooFit2D || fittype==HistFact2D)
  {


    //Get observables
    RooRealVar * B_plus_M = (RooRealVar*) varset->find("B_plus_M");
    RooRealVar * misPT = (RooRealVar*) varset->find("misPT");

    //Get parameters 
    string l1Keename = "l1Kee";
    if (ownexp) l1Keename = "expoConst_"+name;
    RooRealVar * l1Kee = (RooRealVar*) parset->find(l1Keename.c_str());
    RooRealVar * l2Kee = (RooRealVar*) parset->find("l2Kee");
    RooRealVar * l3Kee = (RooRealVar*) parset->find("l3Kee");
    RooRealVar * l4Kee = (RooRealVar*) parset->find("l4Kee");
    RooRealVar * l5Kee = (RooRealVar*) parset->find("l5Kee");

    if (!B_plus_M || !misPT || !l1Kee || !l2Kee || !l3Kee || !l4Kee || !l5Kee)
    {
      cout<<"Varset does not contain variables needed to build the model."<<endl;
      return model;
    }    
    
    if(fixpars)
      {
	l2Kee->setConstant(1);
	l3Kee->setConstant(1);
	l4Kee->setConstant(1);
	l5Kee->setConstant(1);
      }

    if (fittype==RooFit2D)
    {  
      cout<<"Going for fittype: "<< "RooFit2D" <<endl;
      model =  new RooExpOfPolyTimesX(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *misPT, *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee);
    }
    
    if (fittype==HistFact2D)
    {
      // RooBinning m_binning = (RooBinning) B_plus_M->getBinning(m_binning_name.c_str()); 
      // RooBinning pT_binning = (RooBinning) misPT->getBinning(pT_binning_name.c_str()); 

      cout<<"Going for fittype: "<< "HistFact2D" <<endl;
      model =  new RooExpOfPolyTimesXBinned(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *misPT, *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee,
                                            B_plus_M->getBinning(binningname.c_str()), 
                                            misPT->getBinning(binningname.c_str()));
    }
  }

  
  return model;
  
}


void Combinatorial::plot(RooDataSet& dataset, RooAbsPdf& model, string plotname)
{

   TFile f2(plotsfile.c_str(), "UPDATE");

   RooArgSet *var_set = model.getObservables(dataset);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   while((var = (RooRealVar*) iter->Next()))
   {

      frame = var->frame();
      dataset.plotOn(frame);
      model.plotOn(frame);

      plots.push_back(frame);

   }  

   if (!(plots.size())) return;

   TCanvas cComb(("cComb_"+name+"_"+plotname).c_str(), ("cComb_"+name+"_"+plotname).c_str(), 600, 800);
   cComb.Divide(1,2);
   cComb.cd(1); plots[0]->Draw();
   if (plots.size()>1){ 
      cComb.cd(2); plots[1]->Draw();
   }

   cComb.Write();
   f2.Close();
}


