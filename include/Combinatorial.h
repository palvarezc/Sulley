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
#include "TCanvas.h"
#include "RooPlot.h"

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
  Combinatorial(string name_, FitType fittype_, string plotsfile_ );
  virtual ~Combinatorial(){};

  void build_model_from_data(string combfile, string combtree,
                             string combcuts, RooArgSet *varset, RooWorkspace *workspace,
                             bool wantplot=false, string binningname="default");
  
  string get_name()
  {
    return name;    
  }

  void use_own_exp(bool _ownexp)
  {
    ownexp = _ownexp;
  }

  /* void delete_set(RooArgSet *parameters); */

  RooAbsPdf* build_model(RooArgSet *varset, RooArgSet *parset, bool fixpars=false, string binningname="default");

  void generate(string compname, RooArgSet *observables, 
                int nGenComp, RooWorkspace* workspace, 
		RooWorkspace* workspaceGen, bool wantplot=0,
		string binningname="default");

  void create_parameters(RooWorkspace* workspace);  
  void copy_parameters(Combinatorial *other, RooWorkspace* workspace, bool _ownexp);
  void plot(RooDataSet& dataset, RooAbsPdf& model, string plotname);
  
  string name;
  string plotsfile;
  FitType fittype;
  bool ownexp;
  
};






#endif // INCLUDE_COMBINATORIAL_H  
