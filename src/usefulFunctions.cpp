#include "usefulFunctions.h"


string makestring(int sigma){
std::ostringstream strs;
std::setprecision(0); //std::fixed;
strs << setfill('0') << setw(3) << sigma;
std::string legs = strs.str();
return legs;
}

string d2s(double nbr, int nfixed )
{
   stringstream ss;
   if(nfixed>=1) ss<<fixed<<setprecision(nfixed)<<nbr;
   else ss<<nbr;
   return ss.str();
}

string roundToError(valError& ve)
{
   int power(floor(TMath::Log10(ve.err)));
   double todivide(TMath::Power(10, power-2));
   int err3dig(TMath::Nint(ve.err/todivide));
   int nfixed;
   if(err3dig<=354 || err3dig>=950)
   {
      todivide = TMath::Power(10, power-1);
      nfixed = 1-power;
   }
   if(err3dig>=355 && err3dig<=949)
   {
      todivide = TMath::Power(10, power);
      nfixed = 0-power;
   }
   ve.err = todivide*TMath::Nint(ve.err/todivide);
   ve.val = todivide*TMath::Nint(ve.val/todivide);
   string ret(d2s(ve.val, nfixed)+"+-"+d2s(ve.err, nfixed));
   return ret;
}


void makeTableResults(TTree* t, int nGenSignal, int nGenPartReco, int nGenComb, int nGenJpsiLeak, ostream& out )
{
   valError ve;
   TH1F* hStock(0);
   TCanvas canv("canv", "canv", 600, 600);

   t->Draw("nSignal>>hStock");
   hStock = (TH1F*)canv.GetPrimitive("hStock");
   ve.val = hStock->GetMean();
   ve.err = hStock->GetStdDev();
   double fracSig(100*ve.err/ve.val);
   string strSig(roundToError(ve));


   t->Draw("nComb>>h2");
   hStock = (TH1F*)canv.GetPrimitive("h2");
   ve.val = hStock->GetMean();
   ve.err = hStock->GetStdDev();
   double fracBkg(100*ve.err/ve.val);
   string strBkg = roundToError(ve);


   t->Draw("nPartReco>>h3");
   hStock = (TH1F*)canv.GetPrimitive("h3");
   ve.val = hStock->GetMean();
   ve.err = hStock->GetStdDev();
   double fracPartReco(100*ve.err/ve.val);
   string strPartReco = roundToError(ve);

   double fracJpsiLeak(0);
   string strJpsiLeak("0");
   if(nGenJpsiLeak > 0)
   {
      t->Draw("nJpsiLeak>>h4");
      hStock = (TH1F*)canv.GetPrimitive("h4");
      ve.val = hStock->GetMean();
      ve.err = hStock->GetStdDev();
      fracJpsiLeak = 100*ve.err/ve.val;
      strJpsiLeak = roundToError(ve);
   }

   int w(16);
   out<<setw(w)<<" "<<setw(w)<<"nSignal"<<setw(w)<<"nPartReco"<<setw(w)<<"nCombinatorial"<<setw(w)<<"nJpsiLeak"<<endl;
   out<<setw(w)<<"Generated"<<setw(w)<<nGenSignal<<setw(w)<<nGenPartReco<<setw(w)<<nGenComb<<setw(w)<<nGenJpsiLeak<<endl;
   out<<setw(w)<<"Fit"<<setw(w)<<strSig<<setw(w)<<strPartReco<<setw(w)<<strBkg<<setw(w)<<strJpsiLeak<<endl;
   out<<setw(w)<<"FracErr"<<setw(w)<<fracSig<<setw(w)<<fracPartReco<<setw(w)<<fracBkg<<setw(w)<<fracJpsiLeak<<endl;
}


void fillTreeResult(TTree* t, RooFitResult* rfr, bool update)
{

   RooArgList constPars = rfr->constPars();
   RooArgList floatPars = rfr->floatParsFinal();
   int nConst(constPars.getSize());
   int nFloat(floatPars.getSize());

   vector<double> vec(2*nFloat + nConst);

   string name;
   double value;

   t->ResetBranchAddresses();

   for(int i(0); i < nConst; ++i)
   {   
      name = ((RooRealVar*)constPars.at(i))->GetName();
      value = ((RooRealVar*)constPars.at(i))->getValV();
      vec[i] = value;
      if(!update) t->Branch( name.c_str(), &vec[i], (name+"/D").c_str() );
      if(update) t->SetBranchAddress(name.c_str(), &vec[i]);
   }

   for(int i(0); i < nFloat; ++i)
   {   
      name = ((RooRealVar*)floatPars.at(i))->GetName();
      value = ((RooRealVar*)floatPars.at(i))->getValV();
      vec[nConst+i] = value; 
      if(!update) t->Branch( name.c_str(), &vec[nConst+i], (name+"/D").c_str() );
      if(update) t->SetBranchAddress(name.c_str(), &vec[nConst+i]);
   }

   for(int i(0); i < nFloat; ++i)
   {   
      name = ((RooRealVar*)floatPars.at(i))->GetName();
      value = ((RooRealVar*)floatPars.at(i))->getError();
      vec[nConst+nFloat+i] = value;
      if(!update) t->Branch( (name+"_err").c_str(), &vec[nConst+nFloat+i], (name+"_err/D").c_str() );
      if(update) t->SetBranchAddress((name+"_err").c_str(), &vec[nConst+nFloat+i]);
   }

   double edm;
   edm = TMath::Log10(rfr->edm()); 
   if(!update) t->Branch("logedm", &edm, "logedm/D");
   if(update) t->SetBranchAddress("logedm", &edm);

   double status;
   status = rfr->status(); 
   if(!update) t->Branch("status", &status, "status/D");
   if(update) t->SetBranchAddress("status", &status);

   t->Fill();
}
