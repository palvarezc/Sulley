/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooMcorMvisTsallis.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 


RooMcorMvisTsallis::RooMcorMvisTsallis(const char *name, const char *title, 
      RooAbsReal& _Mcor,
      RooAbsReal& _Mvis,
      RooAbsReal& _T,
      RooAbsReal& _n,
      RooAbsReal& _Lambda) :
   RooAbsPdf(name,title), 
   Mcor("Mcor","Mcor",this,_Mcor),
   Mvis("Mvis","Mvis",this,_Mvis),
   T("T","T",this,_T),
   n("n","n",this,_n),
   Lambda("Lambda","Lambda",this,_Lambda)
{ 
} 


RooMcorMvisTsallis::RooMcorMvisTsallis(const RooMcorMvisTsallis& other, const char* name) :  
   RooAbsPdf(other,name), 
   Mcor("Mcor",this,other.Mcor),
   Mvis("Mvis",this,other.Mvis),
   T("T",this,other.T),
   n("n",this,other.n),
   Lambda("Lambda",this,other.Lambda)
{ 
} 



Double_t RooMcorMvisTsallis::evaluate() const 
{ 
   if(Mcor < Mvis) return 1e-13;
   return corMassAtFixedMTsallis()*visMass1D();
} 


Double_t RooMcorMvisTsallis::visMass1D() const
{
   return TMath::Exp(Lambda*Mvis);
}

Double_t RooMcorMvisTsallis::tsallis(Double_t& x) const
{
   double mt(sqrt(Mvis*Mvis+x*x));
   double C( 1. / (n*T+(n-2)*Mvis)  ) ;
   return x*C*TMath::Power( 1+((mt-Mvis)/(n*T)), -n);
}

Double_t RooMcorMvisTsallis::corMassAtFixedMTsallis() const
{
   if(Mcor < Mvis) return 1e-13;

   double arg( (Mcor*Mcor -  Mvis*Mvis) / (2.*Mcor ));
   double factor( (Mcor*Mcor +  Mvis*Mvis)  / ( 2.*Mcor*Mcor ));

   return factor*tsallis(arg);
}
