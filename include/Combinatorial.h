#ifndef INCLUDE_COMBINATORIAL_H 
#define INCLUDE_COMBINATORIAL_H 1

// Include files
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooPTMVis.h"
#include "RooExponential.h"
#include "RooExpBinned.h"
#include "RooAbsPdf.h"
#include "RooExpOfPolyTimesXBinned.h"
#include "RooExpOfPolyTimesX.h"
#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "Utils.h"

/** @class Combinatorial Combinatorial.h include/Combinatorial.h
 *  
 *
 *  @author Paula Alvarez Cartelle
 *  @date   2017-03-10
 */


using namespace std;
using namespace Sulley;
using namespace RooFit;


class Combinatorial {
public: 

  Combinatorial();
  Combinatorial(string name, FitType fittype );
  virtual ~Combinatorial( )
  {
    // if(model!=0 && model!=NULL) delete model;
    if(parcreated) delete_parameters();
    
  }

  void build_model_from_data(string combfile, string combtree,
                             string combcuts, RooArgSet *varset, RooWorkspace *workspace,
                             string binningname=0);
  
  // void fit_to_data(string combfile, string combtree,
  //                  string combcuts, RooArgSet *varset, RooWorkspace *workspace);

  RooArgSet* get_parameters()
  {
    return parameters;    
  }


  string get_name()
  {
    return name;    
  }

  // RooAbsPdf* get_model()
  // {
  //   return model;
  // }


  void delete_parameters();

  RooAbsPdf* build_model(RooArgSet *varset,  string binningname=0, RooArgSet *parset = NULL);
  void generate(string compname, RooArgSet *observables, RooArgSet *parameterstrue, 
                int nGenComp, RooWorkspace* workspaceGen, bool wantplot=0);
  

private:

  string name;
  FitType fittype;
  // RooAbsPdf *model;
  RooArgSet *parameters;
  bool parcreated;

  RooRealVar *B_plus_M;
  RooRealVar *misPT;
  RooRealVar *expoConst;
  RooRealVar *n;
  RooRealVar *T;
  RooRealVar *l1Kee;
  RooRealVar *l2Kee;
  RooRealVar *l3Kee;
  RooRealVar *l4Kee;
  RooRealVar *l5Kee;

  void set_observables(RooArgSet *observables);
  void set_parameters(RooArgSet *varset);
  void set_parameters_truevalue(RooArgSet *varset);
  RooArgSet* create_parameters();
  
  
};






#endif // INCLUDE_COMBINATORIAL_H  
