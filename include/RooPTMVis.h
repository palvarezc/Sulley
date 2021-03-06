/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOPTMVIS
#define ROOPTMVIS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooPTMVis : public RooAbsPdf {
public:
  RooPTMVis() {} ; 
  RooPTMVis(const char *name, const char *title,
	      RooAbsReal& _pT,
	      RooAbsReal& _mVis,
	      RooAbsReal& _T,
	      RooAbsReal& _n,
	      RooAbsReal& _Lambda);
  RooPTMVis(const RooPTMVis& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPTMVis(*this,newname); }
  inline virtual ~RooPTMVis() { }

protected:

  RooRealProxy pT ;
  RooRealProxy mVis ;
  RooRealProxy T ;
  RooRealProxy n ;
  RooRealProxy Lambda ;
  
  Double_t evaluate() const ;

  Double_t visMass1D() const;
  Double_t tsallis() const;


private:

};
 
#endif
