/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooExpOfPolyTimesXBinned.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
#include "TF1.h"

//ClassImp(RooExpOfPolyTimesXBinned) 

 RooExpOfPolyTimesXBinned::RooExpOfPolyTimesXBinned(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _pT,
                        RooAbsReal& _l1,
                        RooAbsReal& _l2,
                        RooAbsReal& _l3,
                        RooAbsReal& _l4,
                        RooAbsReal& _l5,
                        RooBinning& _m_binning,
                        RooBinning& _pT_binning                                                    
                        // Int_t _m_nbins,
                        // Int_t _pT_nbins,
                        // // Double_t _m_min,
                        // // Double_t _m_max,
                        // // Double_t _pT_min,
                        // // Double_t _pT_max
                        // const char* _rangeName                        
) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   pT("pT","pT",this,_pT),
   l1("l1","l1",this,_l1),
   l2("l2","l2",this,_l2),
   l3("l3","l3",this,_l3),
   l4("l4","l4",this,_l4),
   l5("l5","l5",this,_l5),
   m_nbins(_m_binning.numBins()),
   pT_nbins(_pT_binning.numBins()),
   m_lowBound(_m_binning.lowBound()),
   m_highBound(_m_binning.highBound()),
   pT_lowBound(_pT_binning.lowBound()),
   pT_highBound(_pT_binning.highBound())
 { 
   

   std::string gridname = "grid_";
   gridname+= name;

   grid = new TH2D(gridname.c_str(), gridname.c_str(),m_nbins, m_lowBound, m_highBound,
                   pT_nbins, pT_lowBound, pT_highBound);
   

 } 


 RooExpOfPolyTimesXBinned::RooExpOfPolyTimesXBinned(const RooExpOfPolyTimesXBinned& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   pT("pT",this,other.pT),
   l1("l1",this,other.l1),
   l2("l2",this,other.l2),
   l3("l3",this,other.l3),
   l4("l4",this,other.l4),
   l5("l5",this,other.l5),
   m_nbins(other.m_nbins),
   pT_nbins(other.pT_nbins),
   m_lowBound(other.m_lowBound),
   m_highBound(other.m_highBound),
   pT_lowBound(other.pT_lowBound),
   pT_highBound(other.pT_highBound)
 {

   std::string gridname = "grid_";
   if (name) gridname+= name;

   grid = new TH2D(gridname.c_str(), gridname.c_str(),m_nbins, m_lowBound, m_highBound,
                   pT_nbins, pT_lowBound, pT_highBound);


   // grid = (TH2D*) other.grid->Clone(gridname.c_str());
   

 }


// Double_t RooExpOfPolyTimesXBinned::evaluateUnbinned() const
// {
//    if(m<0 || pT<0) return 1e-13;
//    double ret(0);
//    ret = pT*TMath::Exp( - l1*m - l2*pT - l3*pT*pT - l4*m*pT - l5*m*pT*pT);   // + l6*m*pT*pT*pT );
//    return ret;
// }


Double_t RooExpOfPolyTimesXBinned::evaluateBin(double *x, double *p) const
{

  double pTvar = x[0];
  double mMin = p[0];
  double mMax = p[1];
  // double _l1 = p[2];
  // double _l2 = p[3];
  // double _l3 = p[4];
  // double _l4 = p[5];
  // double _l5 = p[6];
  

  double ret(0);
  if(mMin<=0 && mMax<=0) return 0;
  // if(mMin< 0) mMin = 0;
  
  // ret += -((TMath::Exp(-_l1*mMax - pTvar*(_l2 + _l4*mMax + _l3*pTvar + _l5*mMax*pTvar))*pTvar)/( _l1 + pTvar*(_l4 + _l5*pTvar)));
  // ret -= -((TMath::Exp(-_l1*mMin - pTvar*(_l2 + _l4*mMin + _l3*pTvar + _l5*mMin*pTvar))*pTvar)/( _l1 + pTvar*(_l4 + _l5*pTvar)));

  ret += -((TMath::Exp(-l1*mMax - pTvar*(l2 + l4*mMax + l3*pTvar + l5*mMax*pTvar))*pTvar)/( l1 + pTvar*(l4 + l5*pTvar)));      
  ret -= -((TMath::Exp(-l1*mMin - pTvar*(l2 + l4*mMin + l3*pTvar + l5*mMin*pTvar))*pTvar)/( l1 + pTvar*(l4 + l5*pTvar)));     
  
  return ret;
  
}




Double_t RooExpOfPolyTimesXBinned::evaluate() const
{
   if(m<m_lowBound || pT<pT_lowBound || m>m_highBound || pT>pT_highBound) return 1e-13;
   double ret(0);

   int bin_number_m = grid->GetXaxis()->FindBin(m);
   int bin_number_pT = grid->GetYaxis()->FindBin(pT);
   double m_min = grid->GetXaxis()->GetBinLowEdge(bin_number_m);
   double m_max = grid->GetXaxis()->GetBinLowEdge(bin_number_m+1);
   double pT_min = grid->GetYaxis()->GetBinLowEdge(bin_number_pT);
   double pT_max = grid->GetYaxis()->GetBinLowEdge(bin_number_pT+1);


   // cout<<"Nbins in mass: "<<grid->GetNbinsX()<<endl;
   // cout<<"Nbins in pT: "<<grid->GetNbinsY()<<endl;
   // // cout<<"Looking at bin..."<<bin_number<<endl;
   // cout<<"Looking at bin in mass: "<<grid->GetXaxis()->FindBin(m)<<endl;
   // cout<<"Looking at bin in pT: "<<grid->GetYaxis()->FindBin(pT)<<endl;
   // cout<<"Mass min: "<<m_min<<endl;
   // cout<<"Mass max: "<<m_max<<endl;
   // cout<<"misPT min: "<<pT_min<<endl;
   // cout<<"misPT max: "<<pT_max<<endl;
   

   TF1 *f = new TF1("f",this,&RooExpOfPolyTimesXBinned::evaluateBin,pT_min,pT_max,2);
   f->SetParameter(0,m_min);
   f->SetParameter(1,m_max);
   // f->SetParameter(2,l1);
   // f->SetParameter(3,l2);
   // f->SetParameter(4,l3);
   // f->SetParameter(5,l4);
   // f->SetParameter(6,l5);


   // ret = pT*TMath::Exp( - l1*m - l2*pT - l3*pT*pT - l4*m*pT - l5*m*pT*pT);   // + l6*m*pT*pT*pT );
   
   double binWidth = (m_max-m_min)*(pT_max-pT_min);
   ret = f->Integral(pT_min,pT_max)/binWidth;
   

   return ret;
}

Int_t RooExpOfPolyTimesXBinned::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
   // if(matchArgs(allVars, analVars, m)) return 1;
   // if(matchArgs(allVars, analVars, pT)) return 2;
  if(matchArgs(allVars, analVars, m, pT)) return 3;
   return 0;
}

Double_t RooExpOfPolyTimesXBinned::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
   {
      case 1:
         {
            // double mMin(m.min(rangeName));
            // double mMax(m.max(rangeName));
            // double ret(0);
            // if(mMin<=0 && mMax<=0) return 0;
            // if(mMin< 0) mMin = 0;

            // ret += -((TMath::Exp(-l1*mMax - pT*(l2 + l4*mMax + l3*pT + l5*mMax*pT))*pT)/( l1 + pT*(l4 + l5*pT)));      
            // ret -= -((TMath::Exp(-l1*mMin - pT*(l2 + l4*mMin + l3*pT + l5*mMin*pT))*pT)/( l1 + pT*(l4 + l5*pT)));     

            // return ret;

           

         }

      case 2:
         {
            double pTMin(pT.min(rangeName));
            double pTMax(pT.max(rangeName));
            double ret(0);

            if(pTMin<=0 && pTMax<=0) return 0;
            if(pTMin< 0) pTMin = 0;

            ret += (1./(4*TMath::Power((l3 + l5*m), 1.5))) * TMath::Exp(-l1*m - pTMax * (l2 + l4*m + l3*pTMax + l5*m*pTMax))*(-2*sqrt(l3 + l5*m) -  TMath::Exp(TMath::Power(l2 + l4*m + 2*(l3 + l5*m)*pTMax,2)/(4*(l3 + l5*m)))*(l2 + l4*m)*sqrt(M_PI)*TMath::Erf((l2 + l4*m + 2*(l3 + l5*m)*pTMax)/(2*sqrt(l3 + l5*m))));

            ret -= (1./(4*TMath::Power((l3 + l5*m), 1.5))) * TMath::Exp(-l1*m - pTMin * (l2 + l4*m + l3*pTMin + l5*m*pTMin))*(-2*sqrt(l3 + l5*m) -  TMath::Exp(TMath::Power(l2 + l4*m + 2*(l3 + l5*m)*pTMin,2)/(4*(l3 + l5*m)))*(l2 + l4*m)*sqrt(M_PI)*TMath::Erf((l2 + l4*m + 2*(l3 + l5*m)*pTMin)/(2*sqrt(l3 + l5*m))));

            return ret;
         }



   case 3:
     {

       double pTMin(pT.min(rangeName));
       double pTMax(pT.max(rangeName));
       
       double m_min, m_max;
       double pT_min, pT_max;
       
       double ret(0);
       
       TF1 *f = new TF1("f",this,&RooExpOfPolyTimesXBinned::evaluateBin,pTMin,pTMax,2);

       for (int i=1; i<m_nbins+1; i++)
       {
         for (int j=1; j<pT_nbins+1; j++)
         {
           
           m_min = grid->GetXaxis()->GetBinLowEdge(i);
           m_max = grid->GetXaxis()->GetBinLowEdge(i+1);
           pT_min = grid->GetYaxis()->GetBinLowEdge(j);
           pT_max = grid->GetYaxis()->GetBinLowEdge(j+1);
           
           f->SetParameter(0,m_min);
           f->SetParameter(1,m_max);

           ret += f->Integral(pT_min,pT_max);

         }
       }   

       return ret;

     }
     

   }

   assert(0);
   return 0;
}

