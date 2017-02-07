/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOEXPOFPOLYTIMESXBINNED
#define ROOEXPOFPOLYTIMESXBINNED

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TF1.h"
#include "TH2D.h"
#include "RooBinning.h"

class RooExpOfPolyTimesXBinned : public RooAbsPdf {
public:
  RooExpOfPolyTimesXBinned() {} ; 
  RooExpOfPolyTimesXBinned(const char *name, const char *title,
                           RooAbsReal& _m,
                           RooAbsReal& _pT,
                           RooAbsReal& _l1,
                           RooAbsReal& _l2,
                           RooAbsReal& _l3,
                           RooAbsReal& _l4,
                           RooAbsReal& _l5,
                           RooBinning& _m_binning,   
                           RooBinning& _pT_binning);
  RooExpOfPolyTimesXBinned(const RooExpOfPolyTimesXBinned& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpOfPolyTimesXBinned(*this,newname); }
  inline virtual ~RooExpOfPolyTimesXBinned() { 
    // delete grid;
  }

protected:

  RooRealProxy m ;
  RooRealProxy pT ;
  RooRealProxy l1 ;
  RooRealProxy l2 ;
  RooRealProxy l3 ;
  RooRealProxy l4 ;
  RooRealProxy l5 ;
  Int_t m_nbins;
  Int_t pT_nbins;
  Double_t m_lowBound;
  Double_t m_highBound;
  Double_t pT_lowBound;
  Double_t pT_highBound;
 

  TH2D *grid;
  TF1 *binnedFunction;
  
  Double_t evaluate() const ;
  // Double_t evaluateUnbinned() const ;
  Double_t evaluateBin(double *x, double *p) const ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName) const;


private:

  //ClassDef(RooExpOfPolyTimesXBinned,1) // Your description goes here...
};
 
#endif