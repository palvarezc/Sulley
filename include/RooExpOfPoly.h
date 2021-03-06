/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
  * 
  * Compute hte following expression:
  * Exp( l1*m + l2*pT + l3*m*pT + l4*m*pT*pT + l5*pT*pT*pT + l6*m*pT*pT*pT )
 *****************************************************************************/

#ifndef ROOEXPOFPOLY
#define ROOEXPOFPOLY

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooExpOfPoly : public RooAbsPdf {
public:
  RooExpOfPoly() {} ; 
  RooExpOfPoly(const char *name, const char *title,
	      RooAbsReal& _m,
	      RooAbsReal& _pT,
	      RooAbsReal& _l1,
	      RooAbsReal& _l2,
	      RooAbsReal& _l3,
	      RooAbsReal& _l4,
	      RooAbsReal& _l5);
  RooExpOfPoly(const RooExpOfPoly& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpOfPoly(*this,newname); }
  inline virtual ~RooExpOfPoly() { }

protected:

  RooRealProxy m ;
  RooRealProxy pT ;
  RooRealProxy l1 ;
  RooRealProxy l2 ;
  RooRealProxy l3 ;
  RooRealProxy l4 ;
  RooRealProxy l5 ;
  
  Double_t evaluate() const ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

private:

//  ClassDef(RooExpOfPoly,1) // Your description goes here...
};
 
#endif
