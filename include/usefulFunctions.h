#ifndef USEFUL_FUNCTIONS_H
#define USEFUL_FUNCTIONS_H


#include<iostream>
#include<iomanip>
#include<sstream>
#include<string>
#include "TTree.h"
#include "RooFitResult.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "TRandom.h"
#include "TMath.h"
#include<fstream>

using namespace std;


struct valError{
   double err;
   double val;
};

string makestring(int sigma);
string dbl2str(double nbr, int nfixed );
string roundToError(valError& ve);
void makeTableResults(TTree* t, int nGenSignal, int nGenPartReco, int nGenComb, ostream& out = cout );
void fillTreeResult(TTree* t, RooFitResult* rfr, bool update);

#endif
