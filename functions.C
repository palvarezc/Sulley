#include <math.h>
#include "TMath.h"


Double_t RooExpOfPolyTimesXBinnedFunction(double *x, double *p)  
{

double pTvar = x[0];
double mMin = p[0];
double mMax = p[1];
double l1   = p[2];
double l2   = p[3];
double l3   = p[4];
double l4   = p[5];
double l5   = p[6];
 
// double _l1 = p[2];                                                                                                                  // double _l2 = p[3];                                                                                                                
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
