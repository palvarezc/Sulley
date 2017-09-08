#ifndef VALERROR_H
#define VALERROR_H 

#include<fstream>
#include<string>
#include<vector>
#include "usefulFunctions.h"

using namespace std;

class ValError
{
   public:

   double val;
   double err;

   ValError();
   ValError(double _val, double _err);
   ValError(ValError const& ve);


   //ValError& operator=(ValError const& ve);
   ValError& operator*=(ValError const& ve);
   ValError& operator*=(double alpha);
   ValError& operator/=(ValError const& ve);
   ValError operator*(ValError const& ve) const;
   ValError operator/(ValError const& ve) const;
   ValError& operator+=(ValError const& ve);
   ValError operator+(ValError const& ve);
   ValError& operator-=(ValError const& ve);
   ValError operator-(ValError const& ve);
   ValError& operator/=(double alpha);
   ValError operator/(double alpha);
   bool operator==(ValError const& ve) const;


   string roundToError(bool wantLatex = false);

};

ValError operator*(double scalar, ValError const& ve);

ostream& operator<<(ostream&, ValError const&);
void outLatex(ostream& out, ValError const&);


ValError sqrt(ValError ve);
string roundToError(ValError ve);
ValError averageValError(ValError va1, ValError va2 = ValError(), ValError va3 = ValError(), ValError va4 = ValError(),
      ValError va5 = ValError(), ValError va6 = ValError(), ValError va7 = ValError(), ValError va8 = ValError(), ValError va9 = ValError());
ValError averageValError(vector<ValError> vecVa);
ValError getRatio(ValError a, ValError b);
ValError getRatioBinomial(ValError a, ValError b);
ValError getRatioWeightedBinomial(ValError pass, ValError tot);
ValError getProduct(ValError a, ValError b);
ValError getMeanWithStdDev(vector<double> vals);

#endif
