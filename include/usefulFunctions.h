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
string d2s(double nbr, int nfixed = 0 );
string roundToError(valError& ve);
void makeTableResults(TTree* t, int nGenSignal, int nGenPartReco, int nGenComb, ostream& out = cout );
void fillTreeResult(TTree* t, RooFitResult* rfr, bool update);
void makeTableResults(TTree* t, int nGenSignal, int nGenPartReco, int nGenComb, int nGenJpsiLeak, ostream& out );

#endif
