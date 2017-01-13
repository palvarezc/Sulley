#include "usefulFunctions.h"


string makestring(int sigma){
std::ostringstream strs;
std::setprecision(0); //std::fixed;
strs << setfill('0') << setw(3) << sigma;
std::string legs = strs.str();
return legs;
}


string i2s(int nbr)
{
   stringstream ss;
   ss<<nbr;
   return ss.str();
}


string d2s(double nbr, int nfixed )
{
   stringstream ss;
   if(nfixed>=1) ss<<fixed<<setprecision(nfixed)<<nbr;
   else ss<<nbr;
   return ss.str();
}

string roundToError(RooRealVar const& var, bool wantLatex)
{
  valError ve;
  ve.val = var.getVal();
  ve.err = var.getError();
  return roundToError(ve, wantLatex);
}

string roundToError(valError& ve, bool wantLatex)
{
   int power(floor(TMath::Log10(ve.err)));
   double todivide(TMath::Power(10, power-2));
   int err3dig(TMath::Nint(ve.err/todivide));
   int nfixed;
   // if(err3dig<=354 || err3dig>=950)
   // {
      todivide = TMath::Power(10, power-1);
      nfixed = 1-power;
  //   }
  // if(err3dig>=355 && err3dig<=949)
  //  {
  //     todivide = TMath::Power(10, power);
  //     nfixed = 0-power;
  //  }
   ve.err = todivide*TMath::Nint(ve.err/todivide);
   ve.val = todivide*TMath::Nint(ve.val/todivide);
   string ret(d2s(ve.val, nfixed)+"+-"+d2s(ve.err, nfixed));
   if(wantLatex) ret= "$"+d2s(ve.val, nfixed)+"\\pm"+d2s(ve.err, nfixed)+"$";
   return ret;
}

void makeTableResults(string filename, string treename, int nGenSignal, int nGenPartReco, int nGenComb, int nGenJpsiLeak, ostream& out, bool wantLatex )
{
   TFile f(filename.c_str());
   TTree* t = (TTree*)f.Get(treename.c_str());

   makeTableResults(t, nGenSignal, nGenPartReco,  nGenComb, nGenJpsiLeak, out, wantLatex);
   
   f.Close();
}

void makeTableResults(TTree* t, int nGenSignal, int nGenPartReco, int nGenComb, int nGenJpsiLeak, ostream& out, bool wantLatex)//, string titleString )
{
   valError ve;
   TH1F* hStock(0);
   TCanvas canv("canv", "canv", 600, 600);


   t->Draw("nSignal>>hStock");
   cout<<"nsignal: "<<t->GetEntries("nSignal>0")<<endl;
   hStock = (TH1F*)canv.GetPrimitive("hStock");
   cout<<"Mean is: "<<hStock->GetMean()<<endl;

   ve.val = hStock->GetMean();
   ve.err = hStock->GetStdDev();
   cout<<"I break here..."<<endl;
   double fracSig(100*ve.err/ve.val);
   string strSig(roundToError(ve, wantLatex));


   t->Draw("nComb>>h2");
   hStock = (TH1F*)canv.GetPrimitive("h2");
   ve.val = hStock->GetMean();
   ve.err = hStock->GetStdDev();
   double fracBkg(100*ve.err/ve.val);
   string strBkg = roundToError(ve, wantLatex);


   t->Draw("nPartReco>>h3");
   hStock = (TH1F*)canv.GetPrimitive("h3");
   ve.val = hStock->GetMean();
   ve.err = hStock->GetStdDev();
   double fracPartReco(100*ve.err/ve.val);
   string strPartReco = roundToError(ve, wantLatex);

   double fracJpsiLeak(0);
   string strJpsiLeak("0");
   if(nGenJpsiLeak > 0)
   {
      t->Draw("nJpsiLeak>>h4");
      hStock = (TH1F*)canv.GetPrimitive("h4");
      ve.val = hStock->GetMean();
      ve.err = hStock->GetStdDev();
      fracJpsiLeak = 100*ve.err/ve.val;
      strJpsiLeak = roundToError(ve, wantLatex);
   }

   int w(16);
   if(!wantLatex)
   {
      out<<setw(w)<<" "<<setw(w)<<"nSignal"<<setw(w)<<"nPartReco"<<setw(w)<<"nCombinatorial"<<setw(w)<<"nJpsiLeak"<<endl;
      out<<setw(w)<<"Generated"<<setw(w)<<nGenSignal<<setw(w)<<nGenPartReco<<setw(w)<<nGenComb<<setw(w)<<nGenJpsiLeak<<endl;
      out<<setw(w)<<"Fit"<<setw(w)<<strSig<<setw(w)<<strPartReco<<setw(w)<<strBkg<<setw(w)<<strJpsiLeak<<endl;
      out<<setw(w)<<"FracErr"<<setw(w)<<fracSig<<setw(w)<<fracPartReco<<setw(w)<<fracBkg<<setw(w)<<fracJpsiLeak<<endl;
   }

   if(wantLatex)
   {
      out<<"\\documentclass[11pt]{article}"<<endl;
      out<<"\\usepackage[left=0cm,right=0cm,top=0cm,bottom=0cm,paperwidth=12cm,paperheight=2.5cm]{geometry}"<<endl;
      out<<"\\usepackage{xcolor,colortbl}"<<endl;
      out<<"\\begin{document}"<<endl;
      out<<"\\begin{tabular}{r r r r r}"<<endl;
      out<<"  &    \\textbf{Signal}    &   \\textbf{PartReco} & \\textbf{Comb.}  & \\textbf{\\boldmath$J/\\psi$ leak}\\\\"<<endl;
      out<<"\\cellcolor{gray} &  \\multicolumn{4}{c}{\\cellcolor{gray} 1D $M$ fit:} \\\\"<<endl;
      out<<"Generated & $"<<nGenSignal<<"$ & $"<<nGenPartReco<<"$ & $"<<nGenComb<<"$ & $"<<nGenJpsiLeak<<"$ \\\\"<<endl;
      out<<"From fit & "<<strSig<<" & "<<strPartReco<<" & "<<strBkg<<" & "<<strJpsiLeak<<" \\\\"<<endl;
      out<<"Error (\\%) & \\textcolor{red}{$"<<setprecision(3)<<fracSig<<"$} & $"<<setprecision(3)<<fracPartReco<<"$ & $"<<setprecision(3)<<fracBkg<<"$ & $"<<setprecision(3)<<fracJpsiLeak<<"$ \\\\"<<endl;
      out<<"\\end{tabular}"<<endl;
      out<<"\\end{document}"<<endl;
   }
   cout<<"CACA "<<endl;
}


void fillTreeResult(TTree* t, RooFitResult* rfr, bool update, int migradRes, int hesseRes,bool  hasConverged)
{

   RooArgList constPars = rfr->constPars();
   RooArgList floatPars = rfr->floatParsFinal();
   int nConst(constPars.getSize());
   int nFloat(floatPars.getSize());

   vector<double> vec(2*nFloat + nConst);

   string name;
   double value;

   t->ResetBranchAddresses();

   double minNll(rfr->minNll());


   if(!update) t->Branch( "migradRes", &migradRes, "migradRes/I" );
   if(update) t->SetBranchAddress("migradRes", &migradRes);

   if(!update) t->Branch( "hesseRes", &hesseRes, "hesseRes/I" );
   if(update) t->SetBranchAddress("hesseRes", &migradRes);

   if(!update) t->Branch( "minNll", &minNll, "minNll/D" );
   if(update) t->SetBranchAddress("minNll", &minNll);

   if(!update) t->Branch( "hasConverged", &hasConverged, "hasConverged/O" );
   if(update) t->SetBranchAddress("hasConverged", &hasConverged);

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

//void makePull(string nameFile, string nameTree, string nameVar, double centre, string namePlot)
//{
//   TFile f(nameFile.c_str());
//   TTree* t = (TTree*)f.Get(nameTree.c_str());
//
//   TCanvas canv("canv", "canv", 600, 600);
//   TH1F* hpull(0);
//   t->Draw( ("("+nameVar+"-"+d2s(centre)+")/("+nameVar+"_err)>>hpull(100, -6, 6)").c_str());
//
//   hpull = (TH1F*)canv.GetPrimitive("hpull");
//
//   hpull->SetLineWidth(2);
//   
//   hpull->GetXaxis()->SetTitle( "(n-#lambda)/#sigma");
//   hpull->GetYaxis()->SetTitle( "#events");
//
//   hpull->Draw();
//   
//   canv.Print(namePlot.c_str());
//}
