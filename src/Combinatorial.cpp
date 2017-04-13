// Include files 
#include "Combinatorial.h"


Combinatorial::Combinatorial( string name_, FitType fittype_ ):
  name(name_),
  fittype(fittype),
  parcreated(0)
{


}


RooArgSet* Combinatorial::create_parameters()
{


  // RooArgSet *parameters = new RooArgSet(parsetname.c_str());
  // RooArgSet *parameters = new RooArgSet();

  parameters = new RooArgSet(); 

  if (fittype==FitType::RooFit1D || fittype==FitType::HistFact1D)
  {
    expoConst = new RooRealVar("expoConst", "expoConst", -1e-3, -1, 1);
    parameters->add(*expoConst);
  }
  
  if (fittype==FitType::RooFit2D_RooPTMVis)
  {
    expoConst = new RooRealVar("expoConst", "expoConst", -1e-3, -1, 1);
    T = new RooRealVar("T", "T", 97, 0, 200);
    n = new RooRealVar("n", "n", 3.5, 1., 5.5);
    parameters->add(*expoConst);
    parameters->add(*T);
    parameters->add(*n);
  }

  if (fittype==FitType::RooFit2D || fittype==FitType::HistFact2D)
  {
    l1Kee = new RooRealVar("l1Kee", "l1Kee", +2.3e-3, -1e-1, 1e-1);    
    l2Kee = new RooRealVar("l2Kee", "l2Kee", +4.9e-3, -1e-1, 1e-1);
    l3Kee = new RooRealVar("l3Kee", "l3Kee", -4.1e-7, -1e-5, 1e-5);
    l4Kee = new RooRealVar("l4Kee", "l4Kee", -4e-7, -1e-5, 1e-5);
    l5Kee = new RooRealVar("l5Kee", "l5Kee", 1e-10, -1e-9, 1e-9);
    parameters->add(*l1Kee);
    parameters->add(*l2Kee);
    parameters->add(*l3Kee);
    parameters->add(*l4Kee);
    parameters->add(*l5Kee);
  }
  
  parcreated = true;
  return parameters;

}


void Combinatorial::delete_parameters()
{

  RooRealVar *var;  
  TIterator *iter = parameters->createIterator();

  while((var = (RooRealVar*) iter->Next()))
  {
    delete var;
  }

  
  delete parameters;
  parcreated = false;

}



void Combinatorial::set_observables(RooArgSet *varset)
{

  RooRealVar *var;  
  TIterator *iter = varset->createIterator();

  while((var = (RooRealVar*) iter->Next()))
  {
    if ( var->GetName()=="B_plus_M") B_plus_M = var;
    if ( var->GetName()=="misPT") misPT = var;
  }

}

void Combinatorial::set_parameters(RooArgSet *varset)
{

  RooRealVar *var;  
  TIterator *iter = varset->createIterator();

  while((var = (RooRealVar*) iter->Next()))
  {
    if ( var->GetName()=="expoConst") expoConst = var;
    if ( var->GetName()=="n") n = var;
    if ( var->GetName()=="T") T = var;
    if ( var->GetName()=="l1Kee") l1Kee = var;
    if ( var->GetName()=="l2Kee") l2Kee = var;
    if ( var->GetName()=="l3Kee") l3Kee = var;
    if ( var->GetName()=="l4Kee") l4Kee = var;
    if ( var->GetName()=="l5Kee") l5Kee = var;
  }

}


void Combinatorial::set_parameters_truevalue(RooArgSet *varset)
{

  RooRealVar *var;
  TIterator *iter = get_parameters()->createIterator();

  RooRealVar *vartrue;
  TIterator *itertrue = varset->createIterator();

  string varname, vartruename;

  while((var = (RooRealVar*) iter->Next()))
  {
    varname = var->GetName();

    while((vartrue = (RooRealVar*) itertrue->Next()))
    {

      vartruename = vartrue->GetName();
      if ( vartruename=="true"+varname) 
      {
        
        var->setVal(vartrue->getValV());
      }
    }
  }

}

// void Combinatorial::fit_to_data(string combfile, string combtree,
//                  string combcuts, RooArgSet *varset, RooWorkspace *workspace)
// {

//   RooBinning dumbbinning;
//   fit_to_data(combfile, combtree, combcuts,varset,workspace,dumbbinning,dumbbinning);

// }



void Combinatorial::build_model_from_data(string combfile, string combtree, string combcuts, 
                                          RooArgSet *varset,RooWorkspace *workspace,
                                          string binningname)
{

  RooAbsPdf* model = build_model(varset, binningname);
  
  
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
  

  RooDataSet *dataSetComb = new RooDataSet("dataSetComb", "dataSetComb", tComb, *varset, combcuts.c_str());
  model->fitTo(*dataSetComb);
  

  RooArgSet *parset = model->getParameters(dataSetComb);
  TIterator *pariter = parset->createIterator();
  RooRealVar *par = (RooRealVar*) pariter->Next();
  string parname;
  
  workspace->defineSet((name+"_parsettrue").c_str(),"");

  //Save values to workspace
  while (par)
  {
    parname = par->GetName();
    RooRealVar truepar(("true"+parname).c_str(), ("true"+parname).c_str(), par->getVal());
    workspace->import(truepar);
    workspace->extendSet((name+"_parsettrue").c_str(), ("true"+parname).c_str());
    workspace->import(*par);
    par = (RooRealVar*) pariter->Next();
  }

  workspace->defineSet((name+"_parset").c_str(),*parset);
  workspace->saveSnapshot((name+"_calibration").c_str(),*parset,kTRUE) ;

  delete parset;
  delete pariter;
  delete par;
  delete variter;
  delete var;
  delete fComb;  
  delete dataSetComb;
  delete model;
}



void Combinatorial::generate(string compname, RooArgSet *observables, RooArgSet *parameterstrue,
                             int nGenComp, RooWorkspace* workspaceGen, bool wantplot)
{

  RooAbsPdf *model = build_model(observables);

  set_parameters_truevalue(parameterstrue);

  RooAbsPdf::GenSpec* GenSpecComp = model->prepareMultiGen(*observables, RooFit::Extended(1), NumEvents(nGenComp));

  cout<<"Generating "<<compname<<endl;
  RooDataSet* dataGenComp = model->generate(*GenSpecComp);
  dataGenComp->SetName(("dataGen"+compname).c_str()); dataGenComp->SetTitle(("dataGen"+compname).c_str());

  // if (m_wantplot)
  // {
  //   RooDataSet* dataSetComp = (RooDataSet*)workspace->data("dataSetComp");    
  //   PlotShape(*dataSetComp, *dataGenComp, *histPdfComp, opts.plotsfile, ("c"+compname).c_str(), varset);
  // }
  

   

  workspaceGen->import(*dataGenComp);

  delete dataGenComp;
  delete GenSpecComp;
  delete model;
}





RooAbsPdf* Combinatorial::build_model(RooArgSet *varset, string binningname, RooArgSet *parset)
{

  RooAbsPdf *model = NULL;

  set_observables(varset);

  if (!parset && !parcreated) parset = create_parameters();
  if (parset) set_parameters(parset);

  // Single exponential
  if (fittype==FitType::RooFit1D || fittype==FitType::HistFact1D)
  {
    if (!B_plus_M || !expoConst)
    {
      cout<<"Varset does not contain variables needed to build the model."<<endl;
      return model;
    }

    if (fittype==FitType::RooFit1D)
    {
      model =  new RooExponential(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *expoConst);
    }
  
    if (fittype==FitType::RooFit1D)
    {
      model =  new RooExpBinned(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *B_plus_M, *expoConst, 
                                B_plus_M->getBinning(binningname.c_str()));
    }

    
  }
  
  
  // 2D - RooPTMVis
  if (fittype==FitType::RooFit2D_RooPTMVis)
  {

    if (!B_plus_M || !misPT || !expoConst || !n || !T)
    {
      cout<<"Varset does not contain variables needed to build the model."<<endl;
      return model;
    }

    model =  new RooPTMVis(("combPDF_"+name).c_str(), ("combPDF_"+name).c_str(), *misPT, *B_plus_M, *T, *n, *expoConst);


  }



  // 2D - RooPolyTimesX
  if (fittype==FitType::RooFit2D || fittype==FitType::HistFact2D)
  {

    if (!B_plus_M || !misPT || !l1Kee || !l2Kee || !l3Kee || !l4Kee || !l5Kee)
    {
      cout<<"Varset does not contain variables needed to build the model."<<endl;
      return model;
    }    
    
    parameters = new RooArgSet("parameters");
    parameters->add(*l1Kee);
    parameters->add(*l2Kee);
    parameters->add(*l3Kee);
    parameters->add(*l4Kee);
    parameters->add(*l5Kee);

    if (fittype==FitType::RooFit2D)
    {  
      model =  new RooExpOfPolyTimesX("combPDF", "combPDF", *B_plus_M, *misPT, *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee);
    }
    
    if (fittype==FitType::HistFact2D)
    {
      // RooBinning m_binning = (RooBinning) B_plus_M->getBinning(m_binning_name.c_str()); 
      // RooBinning pT_binning = (RooBinning) misPT->getBinning(pT_binning_name.c_str()); 

      model =  new RooExpOfPolyTimesXBinned("combPDF", "combPDF", *B_plus_M, *misPT, *l1Kee, *l2Kee, *l3Kee, *l4Kee, *l5Kee,
                                            B_plus_M->getBinning(binningname.c_str()), 
                                            misPT->getBinning(binningname.c_str()));
    }
  }

  
  return model;
  
}


