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
#include "TH2D.h"

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
                           Int_t _m_nbins,   
                           Int_t _pT_nbins,
                           // Double_t _m_min,   
                           // Double_t _m_max,   
                           // Double_t _pT_min,
                           // Double_t _pT_max,
                           const char* rangeName="");
  RooExpOfPolyTimesXBinned(const RooExpOfPolyTimesXBinned& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpOfPolyTimesXBinned(*this,newname); }
  inline virtual ~RooExpOfPolyTimesXBinned() { }

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

  TH2D *grid;
  
  Double_t evaluate() const ;
  Double_t evaluateUnbinned() const ;
  Double_t evaluateBin(double *x, double *p) const ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName) const;


private:

  //ClassDef(RooExpOfPolyTimesXBinned,1) // Your description goes here...
};
 
#endif