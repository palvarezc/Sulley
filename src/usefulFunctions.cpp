#include "usefulFunctions.h"
#include "TH2I.h"
#include "THStack.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "RooHist.h"
#include "TF1.h"


double s2d(string str)
{
   double ret;
   istringstream is(str);
   is >> ret;
   return ret;
}

string d2s(double d)
{
   string ret;
   ostringstream os;
   os<<d;
   return os.str();
}

string i2s(int i)
{
   string ret;
   ostringstream os;
   os<<i;
   return os.str();
}

string dbl2str(double nbr, int nfixed )
{
   stringstream ss;
   if(nfixed>=1) ss<<fixed<<setprecision(nfixed)<<nbr;
   else ss<<nbr;
   return ss.str();
}


void addWeightBranch(string namefile, string nameweight, string formulaweight)
{
   TFile f(namefile.c_str());
   TTree* t = (TTree*)f.Get("DecayTree");

   t->SetBranchStatus(nameweight.c_str(), 0);

   TTree* t2 = t->CloneTree(0);

   TTreeFormula ttf("ttf", formulaweight.c_str(), t);
   double weight;
   t2->Branch(nameweight.c_str(), &weight, (nameweight+"/D").c_str());

   for(int i(0); i<t->GetEntries(); ++i)
   {
      t->GetEntry(i);
      weight = ttf.EvalInstance();
      t2->Fill();
   }

   t2->Write();
   f.Close();
}

vector<string> loadList(ifstream& inf)
{
   vector<string> ret;
   string stock;
   while(!inf.eof())
   {
      getline(inf, stock);
      if(stock != "" && stock != " " && stock.find("#") == string::npos) ret.push_back(stock);
   }
   return ret;
}

void removeBranchFromTree(string namefile, string namebranch)
{
   TFile f(namefile.c_str(), "update");
   TTree* t = (TTree*)f.Get("DecayTree");

   if(! t->GetBranch(namebranch.c_str()))
   {
      cout<<"WARNING: in function removeBranchFromTree, for file "<<namefile<<" , no branch called "<<namebranch<<endl;
      return;
   }

   t->SetBranchStatus(namebranch.c_str(), 0);
   TTree* t2 = t->CloneTree();
   t2->Write();
   f.Close();
}

double getEntries(string filename, string cut, string weightName, string nametree, bool resetBranches)
{
   TFile f(filename.c_str());
   TTree* t = (TTree*)f.Get(nametree.c_str());
   if(!t)
   {
      cout<<"WARNING: in getEntries, tree "<<nametree<<" in file "<<filename<<" not found, return 0"<<endl;
      f.Close();
      return 0;
   }
   double ret;
   ret = getEntries(t, cut, weightName, resetBranches);
   f.Close();
   return ret;
}

double getEntries(TTree* t, string cut, string weightName, bool resetBranches)
{
   if(cut == "") cut = "1";
   if(weightName == "") return t->GetEntries(cut.c_str());
   if(resetBranches)
   {
      t->SetBranchStatus("*", 0);
      setBranchStatusTTF(t, cut);
      setBranchStatusTTF(t, weightName);
   }

   double sumW(0);
   TTreeFormula ttf("ttf", cut.c_str(), t);
   TTreeFormula ttfW("ttfW", weightName.c_str(), t);

   for(int i(0); i<t->GetEntries(); ++i)
   {
      t->GetEntry(i);
      if(ttf.EvalInstance()) sumW+=ttfW.EvalInstance();
   }
   
   if(resetBranches) t->SetBranchStatus("*", 1);

   return sumW;
}



void getEntriesWithPoisson(vector<double>& entries, unsigned int nPoisson, string filename, string cut, string weightName, string nametree, bool resetBranches)
{
   TFile f(filename.c_str());
   TTree* t = (TTree*)f.Get(nametree.c_str());
   if(!t)
   {
      cout<<"WARNING: in getEntriesWithPoisson, tree "<<nametree<<" in file "<<filename<<" not found, return 0"<<endl;
      f.Close();
   }
   getEntriesWithPoisson(entries, nPoisson, t, cut, weightName, resetBranches);
   f.Close();
}

void getEntriesWithPoisson(vector<double>& entries, unsigned int nPoisson, TTree* t, string cut, string weightName, bool resetBranches)
{
   if(nPoisson > entries.size())
   {
      cerr<<"ERROR: in getEntriesWithPoisson, bad vector size"<<endl;
      return ;
   }

   for(unsigned int i(0); i<nPoisson; ++i) entries[i] = 0;

   if(cut == "") cut = "1";
   if(weightName == "") weightName = "(1)";
   if(resetBranches)
   {
      t->SetBranchStatus("*", 0);
      setBranchStatusTTF(t, cut);
      setBranchStatusTTF(t, weightName);
   }

   unsigned char tab[100];
   t->SetBranchStatus("poissonWeight", 1);
   t->SetBranchAddress("poissonWeight", &tab);

   double w;
   TTreeFormula ttf("ttf", cut.c_str(), t);
   TTreeFormula ttfW("ttfW", weightName.c_str(), t);

   for(int i(0); i<t->GetEntries(); ++i)
   {
      t->GetEntry(i);
      if(ttf.EvalInstance()) 
      {         
         w = ttfW.EvalInstance();
         for(unsigned int p(0); p<nPoisson; ++p)
         {
            entries[p] += w * tab[p];
         }
      }
   }
   
   if(resetBranches) t->SetBranchStatus("*", 1);
}


ValError getEntriesWithError(string filename, string cut, string weightName, string nametree)
{
   TFile f(filename.c_str());
   TTree* t = (TTree*)f.Get(nametree.c_str());
   if(!t)
   {
      cout<<"WARNING: in getEntriesWithError, tree "<<nametree<<" in file "<<filename<<" not found, return 0"<<endl;
      f.Close();
      return ValError(0,0);
   }
   ValError ret;
   ret = getEntriesWithError(t, cut, weightName);
   f.Close();
   return ret;
}


ValError getEntriesWithError(TTree* t, string cut, string weightName)
{
   if(cut == "") cut = "1";
   if(weightName == "") weightName = "1";
   t->SetBranchStatus("*", 0);
   setBranchStatusTTF(t, cut);
   setBranchStatusTTF(t, weightName);
   double w(0);
   double sumW(0);
   double sumW2(0);

   TTreeFormula ttf("ttf", cut.c_str(), t);
   TTreeFormula ttfWeight("ttfWeight", weightName.c_str(), t );

   for(int i(0); i<t->GetEntries(); ++i)
   {
      t->GetEntry(i);
      if(ttf.EvalInstance()) 
      {
         w = ttfWeight.EvalInstance();
         sumW+=w;
         sumW2+= w*w;
      }
   }

   t->SetBranchStatus("*", 1);
   return ValError(sumW, sqrt(sumW2));
}

void plotOneVariable(string filename, 
      string varname, 
      string outputname, 
      int nbin, 
      string cut, 
      string weightExpression, 
      string legEntry, 
      string outputlegname,
      string drawOptions, 
      bool wantLogPlot,
      string nameTree
      )
{
   vector<string> filenames;
   filenames.push_back(filename);
   vector<string> weightExpressions;
   weightExpressions.push_back(weightExpression);
   vector<string> cuts;
   cuts.push_back(cut);
   vector<string> legEntries;
   legEntries.push_back(legEntry);
   
   plotOneVariable(filenames, varname, outputname, nbin, cuts, weightExpressions, outputlegname, legEntries, drawOptions, wantLogPlot, nameTree);
}

void plotOneVariable(vector<string> const& filenames,
      string varname,
      string outputname,
      int nbin,
      string cut,
      vector<string> weightExpressions,
      string outputlegname,
      vector<string> legEntries,
      string drawOptions,
      bool wantLogPlot,
      string nameTree
      )
{
   vector<string> cuts;
   for(unsigned int i(0); i<filenames.size(); ++i) cuts.push_back(cut);
   plotOneVariable(filenames, varname, outputname, nbin, cuts, weightExpressions, outputlegname, legEntries, drawOptions, wantLogPlot, nameTree);
}


void plotOneVariable(vector<string> const& filenames,
      string varname,
      string outputname,
      int nbin,
      vector<string> cuts,
      vector<string> weightExpressions,
      string outputlegname,
      vector<string> legEntries,
      string drawOptions,
      bool wantLogPlot,
      string nameTree)
{
   vector<string> varnames(filenames.size(), varname);

   plotOneVariable(filenames,
      varnames,
      outputname,
      nbin,
      cuts,
      weightExpressions,
      outputlegname,
      legEntries,
      drawOptions,
      wantLogPlot,
      nameTree);


}

void plotOneVariable(vector<string> const& filenames,
      vector<string> varnames,
      string outputname,
      int nbin,
      vector<string> cuts,
      vector<string> weightExpressions,
      string outputlegname,
      vector<string> legEntries,
      string drawOptions,
      bool wantLogPlot,
      string nameTree)
{
   unsigned int n(filenames.size());
   if(n==0) return;

   if(varnames.size() != n)
   {
      cerr<<"ERROR in plotOneVariable: the number of variables does not match the number of files"<<endl;
      return;
   }

   vector<TFile*> vecf(n);
   vector<TTree*> vect(n);

   for(unsigned int i(0); i<n; ++i)
   {
      vecf[i] = new TFile(filenames[i].c_str());
      vect[i] = (TTree*)vecf[i]->Get(nameTree.c_str());
   }

   TH1F* htemp;

   TCanvas tc("tc", "tc", 600, 600);

   double max(0);
   double min(0);

   for(unsigned int i(1); i<n; ++i) separateString(varnames[i], min, max);

   if(!separateString(varnames[0], min, max))
   {
      double maxima[n];
      double minima[n];
      for(unsigned int i(0); i<n; ++i)
      {
         tc.cd();
         vect[i] -> Draw(varnames[i].c_str(), cuts[i].c_str());
         if(-1 != vect[i] -> Draw(varnames[i].c_str(), cuts[i].c_str()))
         {
            htemp = (TH1F*)gPad->GetPrimitive("htemp");
            if(htemp == 0)
            {
               cout<<"Warning: in plotonevariable: could not get a histogram (tree "<<i<<" )"<<endl;
               maxima[i] = -1e10;
               minima[i] = 1e10;
            }
            else
            {
               maxima[i] = htemp->GetXaxis()->GetXmax();
               minima[i] = htemp->GetXaxis()->GetXmin();
            }
         }
         else
         {
            cout<<"Warning: in plotonevariable: "<<varnames[i]<<" not found, plot aborted"<<endl;
            for(unsigned int i(0); i<n; ++i)
            {
               vecf[i]->Close();
               delete vecf[i];
            }
            return;
         }
      }

      max = *max_element(maxima, maxima+n);
      min = *min_element(minima, minima+n); 
   }

   string canvname(cleanNameString(varnames[0]));
   if(canvname.size() > 50) canvname = canvname.substr(0,50);

   for(unsigned int i(0); i<cuts.size(); ++i) if(cuts[i] == "") cuts[i] = "1";
   THStack ths("thstack", varnames[0].c_str());
   bool thsIsEmpty(true);

   vector<TH1F*> htemps(n);


   for(unsigned int i(0); i< weightExpressions.size(); ++i)
   {
      if(weightExpressions[i] == "") weightExpressions[i] = "1";
   }

   if(weightExpressions.size() != n)
   {
      cout<<"Warning: in PlotOnevariable: weightexpression does not have the same size as nbr of trees: all weights set to 1"<<endl;
      weightExpressions = vector<string>(n, "1");
   }

   for(unsigned int i(0); i<n; ++i)
   {
      tc.cd();
      if( -1 != vect[i]->Draw((varnames[i] + ">>htemp_" +i2s(i)+  "("+d2s(nbin)+","+d2s(min)+","+d2s(max)+")").c_str(),
               ( weightExpressions[i]+" * ((("+varnames[i]+") > "+d2s(min)+" && ("+varnames[i]+") < "+d2s(max)+") && ("+cuts[i]+"))" ).c_str()) )
      {
         htemps[i] = (TH1F*)tc.GetPrimitive(("htemp_"+i2s(i)).c_str());
         if(htemps[i] == 0)
         {
            cout<<"Warning: in plotonevariable: could not get a histogram"<<endl;
         }
         if(htemps[i] && htemps[i]->Integral() !=0)
         {
            htemps[i]->SetLineColor(i+1);
            htemps[i]->SetMarkerColor(i+1);
            if(i == 4) htemps[i]->SetLineColor(46);
            if(i == 4) htemps[i]->SetMarkerColor(46);
            htemps[i]->SetLineWidth(2);
            htemps[i]->Sumw2();
            htemps[i]->Scale(1./htemps[i]->Integral(), "width");
            ths.Add(htemps[i]);
            thsIsEmpty = false;
         }
      }
      else
      {
         cout<<"WARNING: in plotonevariable, "<<varnames[i]<<" could not be plotted"<<endl;
      }
   }

   TCanvas canv(canvname.c_str(), canvname.c_str(), 600, 600);
   if(wantLogPlot) canv.SetLogy(true);
   if(!thsIsEmpty)
   {
      ths.Draw(("nostack"+drawOptions).c_str());
      ths.GetXaxis()->SetTitle((varnames[0]).c_str());
      if(outputname.find(".pdf") != string::npos || outputname.find(".root") != string::npos) canv.Print(outputname.c_str());
      else
      {
         outputname = outputname + canvname + ".pdf";
         canv.Print(outputname.c_str());
      }
   }
   if(thsIsEmpty) cout<<"In plotOneVariable, nothing to print for "<<varnames[0]<<endl;

   if(legEntries.size() == n && outputlegname.size() > 0)
   {
      TLegend leg(0.1, 0.1, 0.9, 0.9);
      for(unsigned int i(0); i < n; ++i)
      {
         if(htemps[i] && htemps[i]->Integral() > 0) leg.AddEntry(htemps[i], legEntries[i].c_str(), "l");
      }
      leg.SetLineColor(0);
      leg.SetFillColor(0);
      TCanvas canvleg("legend", "legend", 600,600);
      leg.Draw();

      if(outputlegname.find(".pdf") != string::npos) canvleg.Print(outputlegname.c_str());
      else
      {
         outputlegname = outputlegname  + "A_legend.pdf";
         canvleg.Print(outputlegname.c_str());
      }
   }

   for(unsigned int i(0); i<n; ++i)
   {
      vecf[i]->Close();
      delete vecf[i];
   }

}


bool separateString(string& str1, double& min, double& max)
{
   size_t frst;
   frst = str1.find(";");
   if(frst == string::npos)
   {
       return false;
   }

   size_t scnd;
   scnd = str1.find(";", frst+1);
   if(scnd == string::npos) return false;

   string str2, str3;

   str2 = str1.substr(frst+1, scnd-frst-1);
   str3 = str1.substr(scnd+1);
   str1 = str1.substr(0, frst);

   min = s2d(str2);
   max = s2d(str3);
   return true;
}


bool separateString(string const& str1, string& part1, string& part2)
{
   size_t frst;
   frst = str1.find(";");
   if(frst == string::npos)
   {
       return false;
   }


   string str2, str3;

   part1 = str1.substr(0, frst);
   part2 = str1.substr(frst+1);

   return true;
}

double getCorrelationCoefficient(TTree*& t, string var1name, string var2name, string cuts, bool draw, string plotname)
{
   TH2F th2f("th2f", "th2f", 50,-1e5,1e5,50,-1e5,1e5);

   double var1min(-3);
   double var1max(-3); 
   double var2min(-3);
   double var2max(-3);

   bool isSeparable1;
   bool isSeparable2;

   t->SetBranchStatus("*",0);
   if(cuts != "") setBranchStatusTTF(t, cuts);
   if(cuts == "") cuts = "1"; 

   isSeparable1 = separateString(var1name, var1min, var1max);
   isSeparable2 = separateString(var2name, var2min, var2max);

   if(!isSeparable1 && !isSeparable2) th2f.SetCanExtend(TH1::kAllAxes);
   if(!isSeparable1 && isSeparable2)
   {
      if(draw) th2f.SetBins(50, -1e5, 1e5, 50, var2min, var2max);
      if(draw) th2f.SetCanExtend(TH1::kXaxis);
      cuts += (" && " + var2name + " > " + d2s(var2min) + " && " + var2name + " < " + d2s(var2max));
   }
   if(isSeparable1 && !isSeparable2)
   {
      if(draw) th2f.SetBins(50, var1min, var1max, 50, -1e5, 1e5);
      if(draw) th2f.SetCanExtend(TH1::kYaxis);
      cuts += (" && " + var1name + " > " + d2s(var1min) + " && " + var1name + " < " + d2s(var1max));
   }
   if(isSeparable1 && isSeparable2)
   {
      if(draw) th2f.SetBins(50, var1min, var1max, 50, var2min, var2max);
      cuts += (" && " + var2name + " > " + d2s(var2min) + " && " + var2name + " < " + d2s(var2max));
      cuts += (" && " + var1name + " > " + d2s(var1min) + " && " + var1name + " < " + d2s(var1max));
   }


   if(!setBranchStatusTTF(t, var1name)) 
   {
      cout<<"WARNING: IN getCorrelationCoefficient, variable "<<var1name<<" not found, return 0"<<endl;
      return 0;
   }
   if(!setBranchStatusTTF(t, var2name)) 
   {
      cout<<"WARNING: IN getCorrelationCoefficient, variable "<<var2name<<" not found, return 0"<<endl;
      return 0;
   }

   double sumx(0);
   double sumy(0);
   double sumxy(0);
   double sumxx(0);
   double sumyy(0);
   double var1(0);
   double var2(0);
   double tot(0);

   TTreeFormula* ttf = new TTreeFormula("ttf", cuts.c_str(), t);
   TTreeFormula* var1ttf = new TTreeFormula("var1ttf", var1name.c_str(), t);
   TTreeFormula* var2ttf = new TTreeFormula("var2ttf", var2name.c_str(), t);

   for(int i(0); i<t->GetEntries(); ++i)
   {
      t->GetEntry(i);
      if(ttf->EvalInstance())
      {
         var1 = var1ttf->EvalInstance();
         var2 = var2ttf->EvalInstance();
         sumx += var1;
         sumy += var2;
         sumxy += var1*var2;
         sumxx += var1*var1;
         sumyy += var2*var2;
         tot += 1;
         if(draw) th2f.Fill(var1, var2);
      }
   }

   if(draw) 
   {
      TCanvas canv("canv2", "canv2", 600, 600);
      th2f.GetXaxis()->SetTitle(var1name.c_str());
      th2f.GetYaxis()->SetTitle(var2name.c_str());
      th2f.SetTitle(( var1name+"_vs_"+var2name).c_str());
      th2f.Draw("colz");
      //canv.Print( (filename+cleanNameString(var1name+"_vs_"+var2name)+".pdf").c_str());
      canv.Print(plotname.c_str());
   }

   double correl = (sumxy - sumx*sumy/tot )/ sqrt( (sumxx -sumx*sumx/tot )*(sumyy - sumy*sumy/tot ) ); 

   t->SetBranchStatus("*",1);
   delete ttf;
   delete var1ttf;
   delete var2ttf;
   if(correl != correl)
   {
      cout<<"WARNING: for variables "<<var1name<<" and "<<var2name<<", correlation coef is nan. Return 1.02"<<endl;
      correl = 1.02;
   }
   return correl;
}

double getCorrelationCoefficient(string filename, string var1name, string var2name, string cuts, bool draw, string dirname)
{
   double ret;
   TFile f(filename.c_str());
   TTree* t = (TTree*)f.Get("DecayTree");
   ret = getCorrelationCoefficient(t, var1name, var2name, cuts, draw, dirname);
   f.Close();
   return ret;
}

void getCorrelationMatrix(vector<string> vars, string namefile, string nameoutput, string cuts, bool draw2d)
{
   TFile f(namefile.c_str());
   TTree* t = (TTree*)f.Get("DecayTree");
   unsigned int n(vars.size());
   TH2I megamatrix("megamatrix", "Mega Matrix", n,0,n,n,0,n);
   megamatrix.SetStats(0);

   double dummy1, dummy2;
   string stringdummy;

   for(unsigned int i(0); i<n; ++i)
   {
      stringdummy = vars[i];
      separateString(stringdummy, dummy1, dummy2);
      megamatrix.GetXaxis()->SetBinLabel(i+1, stringdummy.c_str()); 
      megamatrix.GetYaxis()->SetBinLabel(i+1, stringdummy.c_str()); 
   }

   megamatrix.GetXaxis()->SetLabelSize(0.015);
   megamatrix.GetYaxis()->SetLabelSize(0.015);
   megamatrix.GetXaxis()->LabelsOption("v");

   double coef(0);
   int count(0);
   int tot(0.5*n*(n+1));

   cout<<"computation of big matrix"<<endl;
   for(unsigned int i(0); i<n;++i)
   {
      for(unsigned int j(i); j<n; ++j)
      {
         if(count % (tot/10) == 0) cout<<100*count/tot<<" \% processed"<<endl;
         coef = getCorrelationCoefficient(t, vars[i], vars[j], cuts, draw2d, nameoutput);
         megamatrix.SetBinContent(i+1, j+1, 100*coef);
         megamatrix.SetBinContent(j+1, i+1, 100*coef);
         ++count;
      }
   }

   TCanvas canv("canv", "canv", 1200, 1200);
   canv.SetBottomMargin(0.2);
   canv.SetLeftMargin(0.2);
   megamatrix.SetMarkerSize(0.7);
   megamatrix.Draw("textcolz");

   canv.Print( nameoutput.c_str());
   f.Close();
}

bool setBranchStatusTTF(TTree* t, string cuts)
{
  TObjArray* array = t->GetListOfBranches();
  int n(array->GetEntries());
  string name;
  bool ret(false);
  for(int i(0); i<n; ++i)
  {
     name = ((*array)[i])->GetName();
     if(cuts.find(name) != string::npos)
     {
        t->SetBranchStatus(name.c_str(), 1);
        ret = true;
     }
  }
  return ret;
}

string cleanNameString(string str)
{
   size_t st(0);

   string toreplace(" */#.,[]{}()+-=:");

   for(unsigned int i(0); i<toreplace.size(); ++i)
   {
      while(st != string::npos)
      {
         st = str.find(toreplace[i]);
         if(st!=string::npos) str = str.replace(st, 1, "_");
      }
      st = 0;
   }

   return str;
}

double gamma(double* x, double* param)
{
   return param[0]*TMath::GammaDist(x[0], param[1], param[2], param[3]);
}

double gammaDataNorm(double SPDHits)
{
   double tab[4];
   tab[0] = 33.76956;
   tab[1] = 4.70346;
   tab[2] = -3.45522e-8;
   tab[3] = 123.049;
   return gamma(&SPDHits, tab);
}

double gammaMCNorm(double SPDHits)
{
   double tab[4];
   tab[0] = 20.15166;
   tab[1] = 4.15994;
   tab[2] = -6.73319;
   tab[3] = 62.8926;
   return gamma(&SPDHits, tab);
}


struct EntEvRun{
   unsigned int entry;
   unsigned long long int eventNumber;
   unsigned int runNumber;
};

bool compare(EntEvRun eer1, EntEvRun eer2) { return eer1.eventNumber>eer2.eventNumber; }

void getDoubleEvents(vector<pair<unsigned int, unsigned int> >& doubleEntries, vector<int>& notDoubleEntries, string nameTree )
{
   doubleEntries.clear();
   notDoubleEntries.clear();

   TFile f(nameTree.c_str());
   TTree* t = (TTree*)f.Get("DecayTree");   

   unsigned long long int currentEventNumber;
   unsigned int currentRunNumber;

   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("eventNumber", 1);
   t->SetBranchStatus("runNumber", 1);
   t->SetBranchAddress("eventNumber", &currentEventNumber);
   t->SetBranchAddress("runNumber", &currentRunNumber);

   vector<EntEvRun> toSort(t->GetEntries());
   EntEvRun currentEntEvRun;

   for(int currentEntry(0); currentEntry < t->GetEntries(); ++currentEntry)
   {
      t->GetEntry(currentEntry);
      currentEntEvRun.entry = currentEntry;
      currentEntEvRun.eventNumber = currentEventNumber;
      currentEntEvRun.runNumber = currentRunNumber;
      toSort[currentEntry] = currentEntEvRun;
   }

   sort(toSort.begin(), toSort.end(), compare);

   pair<unsigned int, unsigned int> currentPair;

   for(unsigned int i(0); i< toSort.size()-1; ++i)
   {
      if((toSort[i].eventNumber == toSort[i+1].eventNumber) && (toSort[i].runNumber == toSort[i+1].runNumber))
      {
         currentPair.first = toSort[i].entry;
         currentPair.second = toSort[i+1].entry;
         doubleEntries.push_back(currentPair);
      }
      else notDoubleEntries.push_back(toSort[i].entry);
   }
   sort(notDoubleEntries.begin(), notDoubleEntries.end());
}

void addBranchInTree(string namefile, vector<string> nameNewBranch, vector<string> formulaNewBranch, string nameTree)
{
   if(nameNewBranch.size() != formulaNewBranch.size() || nameNewBranch.size() == 0 || formulaNewBranch.size() == 0)
   {
      cout<<"ERROR: in addBranchInTree, invalid sizes"<<endl;
      return;
   }
   int n(formulaNewBranch.size());

   string namefileInter(namefile);
   namefileInter.insert(namefileInter.size()-5, "ToErase");

   TFile f(namefile.c_str());//,"update" );
   TTree* t = (TTree*)f.Get(nameTree.c_str());

   vector<TTreeFormula*> ttfs(n);
   vector<double> valuesNewBranch(n);

   cout<<endl<<"====================="<<endl<<"Errors from TTreeFormula don't matter from here"<<endl;
   for(int i(0); i<n; ++i)
   {
      t->SetBranchStatus((nameNewBranch[i]).c_str(), 0);
   }
   cout<<"====================="<<endl<<"From now on errors matter"<<endl<<endl;

   TFile f2(namefileInter.c_str(), "RECREATE");

   TTree* t2 = t->CloneTree(0);

   for(int i(0); i<n; ++i)
   {
      ttfs[i] = new TTreeFormula(("ttf"+i2s(i)).c_str(), (formulaNewBranch[i]).c_str(), t) ;
      t2->Branch( (nameNewBranch[i]).c_str(), &valuesNewBranch[i], (nameNewBranch[i]+"/D").c_str());
   }

   int nEntries(t->GetEntries());
   cout<<"Filling the tree with the new branch... "<<flush;
   for(int i(0); i<nEntries; ++i)
   {
      t->GetEntry(i);
      if(i % (nEntries/10) == 0) cout<<" "<<ceil(100*i/(1.*nEntries))<<"\% "<<flush;
      for(int j(0); j<n; ++j) valuesNewBranch[j] = ttfs[j]->EvalInstance();
      t2->Fill();
   }
   cout<<endl;

   for(int i(0); i<n; ++i) delete ttfs[i];

   t2->Write("", TObject::kOverwrite);
   f.Close();
   f2.Close();

   rename(namefileInter.c_str(), namefile.c_str() );
}


void addBranchInTree(string namefile, string nameNewBranch, string formulaNewBranch)
{
   TFile f(namefile.c_str(),"update" );
   TTree* t = (TTree*)f.Get("DecayTree");
   
   cout<<"====================="<<endl<<"Errors from TTreeFormula don't matter from here"<<endl;
   t->SetBranchStatus(nameNewBranch.c_str(), 0);
   cout<<"====================="<<endl<<"From now on errors matter"<<endl;

   TTreeFormula ttf("ttf", formulaNewBranch.c_str(), t);
   double valueNewBranch;
   TTree* t2 = t->CloneTree(0);

   t2->Branch(nameNewBranch.c_str(), &valueNewBranch, (nameNewBranch+"/D").c_str());


   int nEntries(t->GetEntries());
   cout<<"Filling the tree with the new branch... "<<flush;
   for(int i(0); i<nEntries; ++i)
   {
      t->GetEntry(i);
      if(i % (nEntries/10) == 0) cout<<" "<<ceil(100*i/(1.*nEntries))<<"\% "<<flush;
      valueNewBranch = ttf.EvalInstance();
      t2->Fill();
   }
   cout<<endl;
   t2->Write("", TObject::kOverwrite);
   f.Close();
}



void mergeTrees(string targetName, string treeFile1, string treeFile2, string treeFile3, string treeFile4, string treeFile5, string treeFile6)
{
   TChain chain("DecayTree", "DecayTree");
   if(treeFile1 != "") chain.Add( (treeFile1+"/DecayTree").c_str());
   if(treeFile2 != "") chain.Add( (treeFile2+"/DecayTree").c_str());
   if(treeFile3 != "") chain.Add( (treeFile3+"/DecayTree").c_str());
   if(treeFile4 != "") chain.Add( (treeFile4+"/DecayTree").c_str());
   if(treeFile5 != "") chain.Add( (treeFile5+"/DecayTree").c_str());
   if(treeFile6 != "") chain.Add( (treeFile6+"/DecayTree").c_str());

   cout<<"merging..."<<endl;
   chain.Merge(targetName.c_str());
}

void mergeTrees(string targetName, vector<string> fileNames, string tupleName)
{
   TChain chain("DecayTree", "DecayTree");
   for(unsigned int i(0); i<fileNames.size(); ++i) chain.Add( (fileNames[i]+"/"+tupleName).c_str() ); 

   cout<<"merging... "<<endl;
   chain.Merge(targetName.c_str());
}


double getCosThetaL(bool bzero,
      const TLorentzVector& muplus,
      const TLorentzVector& muminus,
      const TLorentzVector& kaon)
{
   //set up boost vectors                                                                                                                                                           
      TLorentzVector b = muplus + muminus + kaon;
   TLorentzVector mumu = muplus + muminus;
   TVector3 mumuboost(-mumu.BoostVector());
   TVector3 bboost(-b.BoostVector());

   //determine costhetal                                                                                                                                                            
   TLorentzVector muminusd(muminus);
   muminusd.Boost(mumuboost);
   TLorentzVector muplusd(muplus);
   muplusd.Boost(mumuboost);
   TLorentzVector bd( b );
   bd.Boost(mumuboost);
   TLorentzVector kd(kaon);
   kd.Boost(mumuboost);

   double costhetal(0.0);
   if (bzero)
      costhetal = cos(muminusd.Vect().Angle(kd.Vect()));
   else
      costhetal = cos(muplusd.Vect().Angle(kd.Vect()));
   return costhetal;
}

void divideByBinWidth(TH1* h)
{
   if(!h)
   {
      cerr<<"ERROR in function divideByBinWidth, histogram not valid"<<endl;
      return;
   }
   int nbins(h->GetNbinsX());
   double content;
   double error;

   for(int i(1); i<=nbins; ++i)
   {
      content = h->GetBinContent(i);
      error = h->GetBinError(i);
      h->SetBinContent( i, content / h->GetXaxis()->GetBinWidth(i) );          
      h->SetBinError( i, error / h->GetXaxis()->GetBinWidth(i) );          
   }
}



// all users - please change the name of this file to lhcbStyle.C
// // Commits to lhcbdocs svn of .C files are not allowed
void lhcbStyle(double fontsize)
{
   // define names for colours
//   Int_t black  = 1;
//   Int_t red    = 2;
//   Int_t green  = 3;
//   Int_t blue   = 4;
//   Int_t yellow = 5; 
//   Int_t magenta= 6;
//   Int_t cyan   = 7;
//   Int_t purple = 9;


   ////////////////////////////////////////////////////////////////////
   // PURPOSE:
   //
   // This macro defines a standard style for (black-and-white) 
   // "publication quality" LHCb ROOT plots. 
   //
   // USAGE:
   //
   // Include the lines
   //   gROOT->ProcessLine(".L lhcbstyle.C");
   //   lhcbStyle();
   // at the beginning of your root macro.
   //
   // Example usage is given in myPlot.C
   //
   // COMMENTS:
   //
   // Font:
   // 
   // The font is chosen to be 132, this is Times New Roman (like the text of
   //  your document) with precision 2.
   //
   // "Landscape histograms":
   //
   // The style here is designed for more or less square plots.
   // For longer histograms, or canvas with many pads, adjustements are needed. 
   // For instance, for a canvas with 1x5 histograms:
   //  TCanvas* c1 = new TCanvas("c1", "L0 muons", 600, 800);
   //  c1->Divide(1,5);
   //  Adaptions like the following will be needed:
   //  gStyle->SetTickLength(0.05,"x");
   //  gStyle->SetTickLength(0.01,"y");
   //  gStyle->SetLabelSize(0.15,"x");
   //  gStyle->SetLabelSize(0.1,"y");
   //  gStyle->SetStatW(0.15);
   //  gStyle->SetStatH(0.5);
   //
   // Authors: Thomas Schietinger, Andrew Powell, Chris Parkes, Niels Tuning
   // Maintained by Editorial board member (currently Niels)
   ///////////////////////////////////////////////////////////////////

   // Use times new roman, precision 2 
   Int_t lhcbFont        = 132;  // Old LHCb style: 62;
   // Line thickness
   Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
   // Text size
   Double_t lhcbTSize    = fontsize;//0.06; 

   // use plain black on white colors
   gROOT->SetStyle("Plain"); 
   TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");

   //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

   lhcbStyle->SetFillColor(1);
   lhcbStyle->SetFillStyle(1001);   // solid
   lhcbStyle->SetFrameFillColor(0);
   lhcbStyle->SetFrameBorderMode(0);
   lhcbStyle->SetPadBorderMode(0);
   lhcbStyle->SetPadColor(0);
   lhcbStyle->SetCanvasBorderMode(0);
   lhcbStyle->SetCanvasColor(0);
   lhcbStyle->SetStatColor(0);
   lhcbStyle->SetLegendBorderSize(0);
   lhcbStyle->SetLegendFont(132);

   // If you want the usual gradient palette (blue -> red)
   lhcbStyle->SetPalette(1);
   // If you want colors that correspond to gray scale in black and white:
//   int colors[8] = {0,5,7,3,6,2,4,1};
//   lhcbStyle->SetPalette(8,colors);

   // set the paper & margin sizes
   lhcbStyle->SetPaperSize(20,26);
   lhcbStyle->SetPadTopMargin(0.05);
   lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
   lhcbStyle->SetPadBottomMargin(0.16);
   lhcbStyle->SetPadLeftMargin(0.14);

   // use large fonts
   lhcbStyle->SetTextFont(lhcbFont);
   lhcbStyle->SetTextSize(lhcbTSize);
   lhcbStyle->SetLabelFont(lhcbFont,"x");
   lhcbStyle->SetLabelFont(lhcbFont,"y");
   lhcbStyle->SetLabelFont(lhcbFont,"z");
   lhcbStyle->SetLabelSize(lhcbTSize,"x");
   lhcbStyle->SetLabelSize(lhcbTSize,"y");
   lhcbStyle->SetLabelSize(lhcbTSize,"z");
   lhcbStyle->SetTitleFont(lhcbFont);
   lhcbStyle->SetTitleFont(lhcbFont,"x");
   lhcbStyle->SetTitleFont(lhcbFont,"y");
   lhcbStyle->SetTitleFont(lhcbFont,"z");
   lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
   lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
   lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

   // use medium bold lines and thick markers
   lhcbStyle->SetLineWidth(lhcbWidth);
   lhcbStyle->SetFrameLineWidth(lhcbWidth);
   lhcbStyle->SetHistLineWidth(lhcbWidth);
   lhcbStyle->SetFuncWidth(lhcbWidth);
   lhcbStyle->SetGridWidth(lhcbWidth);
   lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
   lhcbStyle->SetMarkerStyle(20);
   lhcbStyle->SetMarkerSize(1.0);

   // label offsets
   lhcbStyle->SetLabelOffset(0.010,"X");
   lhcbStyle->SetLabelOffset(0.010,"Y");

   // by default, do not display histogram decorations:
   lhcbStyle->SetOptStat(0);  
   //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
   // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
   lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
   lhcbStyle->SetOptTitle(0);
   lhcbStyle->SetOptFit(0);
   //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
   //titles
   lhcbStyle->SetTitleOffset(0.95,"X");
   lhcbStyle->SetTitleOffset(0.95,"Y");
   lhcbStyle->SetTitleOffset(1.2,"Z");
   lhcbStyle->SetTitleFillColor(0);
   lhcbStyle->SetTitleStyle(0);
   lhcbStyle->SetTitleBorderSize(0);
   lhcbStyle->SetTitleFont(lhcbFont,"title");
   lhcbStyle->SetTitleX(0.0);
   lhcbStyle->SetTitleY(1.0); 
   lhcbStyle->SetTitleW(1.0);
   lhcbStyle->SetTitleH(0.05);

   // look of the statistics box:
   lhcbStyle->SetStatBorderSize(0);
   lhcbStyle->SetStatFont(lhcbFont);
   lhcbStyle->SetStatFontSize(0.05);
   lhcbStyle->SetStatX(0.9);
   lhcbStyle->SetStatY(0.9);
   lhcbStyle->SetStatW(0.25);
   lhcbStyle->SetStatH(0.15);

   // put tick marks on top and RHS of plots
//   lhcbStyle->SetPadTickX(1);
//   lhcbStyle->SetPadTickY(1);

   // histogram divisions: only 5 in x to avoid label overlaps
   lhcbStyle->SetNdivisions(505,"x");
   lhcbStyle->SetNdivisions(510,"y");

   gROOT->SetStyle("lhcbStyle");
   gROOT->ForceStyle();

   // add LHCb label
   TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
         0.87 - gStyle->GetPadTopMargin(),
         gStyle->GetPadLeftMargin() + 0.20,
         0.95 - gStyle->GetPadTopMargin(),
         "BRNDC");
   lhcbName->AddText("LHCb");
   lhcbName->SetFillColor(0);
   lhcbName->SetTextAlign(12);
   lhcbName->SetBorderSize(0);

   TText *lhcbLabel = new TText();
   lhcbLabel->SetTextFont(lhcbFont);
   lhcbLabel->SetTextColor(1);
   lhcbLabel->SetTextSize(lhcbTSize);
   lhcbLabel->SetTextAlign(12);

   TLatex *lhcbLatex = new TLatex();
   lhcbLatex->SetTextFont(lhcbFont);
   lhcbLatex->SetTextColor(1);
   lhcbLatex->SetTextSize(lhcbTSize);
   lhcbLatex->SetTextAlign(12);

   cout << "-------------------------" << endl;  
   cout << "Set LHCb Style - Feb 2012" << endl;
   cout << "-------------------------" << endl;  

}

void makePull(TH1* pull, TH1 const& h1, TH1 const& h2, bool wantUniform)
{
   if(!pull)
   {
      cerr<<"ERROR: in makePull, pointer not valid"<<endl;
      return;
   }
   
   //set the bins correctly

   pull->Reset();
   
   int nBins(h1.GetNbinsX());
   if(wantUniform) pull->SetBins(nBins, 0, 1);

   if(!wantUniform)
   {
      double binTabs[nBins+1];

      for(int i(0); i<=nBins; ++i) binTabs[i] = h1.GetXaxis()->GetBinUpEdge(i);

      pull->SetBins(nBins, binTabs);
   }

   //fill the histo
   
   ValError ve1;
   ValError ve2;
   ValError delta;

   for(int i(1); i<=nBins; ++i)
   {
      ve1.val = h1.GetBinContent(i);
      ve2.val = h2.GetBinContent(i);
      ve1.err = h1.GetBinError(i);
      ve2.err = h2.GetBinError(i);
      delta = ve2-ve1;
      if(delta.err != 0)
      {
         pull->SetBinContent(i, delta.val/delta.err);
         pull->SetBinError(i, 1.);
      }
      else
      {
         pull->SetBinContent(i, 1);
         pull->SetBinError(i, 0);
      }
   }
    
}



void copyTreeWithNewVar(TTree*& tNew, TTree* tOld, string cut, string formula, string newVarName)
{
   vector<string> formulas;
   vector<string> newVarNames;

   if(formula != "" && newVarName != "")
   {
      formulas.push_back(formula);
      newVarNames.push_back(newVarName);
   }
   
   copyTreeWithNewVars(tNew, tOld, cut, formulas, newVarNames);

//   setBranchStatusTTF(tOld, cut);
//   setBranchStatusTTF(tOld, formula);
//
//   TTreeFormula newVarTTF("newVarTTF", formula.c_str(), tOld);
//   if(cut == "") cut = "1";
//   TTreeFormula cutTTF("cutTTF", cut.c_str(), tOld);
//   TTreeFormula varTTF("varTTF", formula.c_str(), tOld);
//
//
//   tNew = tOld->CloneTree(0);
//
//   double newVarVal(0);
//   if(formula != "") tNew->Branch(newVarName.c_str(), &newVarVal, (newVarName+"/D").c_str() );
//
//   cout<<"Filling tree with new var... ";
//
//   for(int i(0); i<tOld->GetEntries(); ++i)
//   {
//      tOld->GetEntry(i);
//
//      if(cutTTF.EvalInstance() )
//      {
//         if(formula != "") newVarVal = newVarTTF.EvalInstance();
//         tNew->Fill();
//      }
//   }
//
//   cout<<"done! New tree has: "<<tNew->GetEntries()<<" entries."<<endl;
}


void copyTreeWithNewVars(string fileNewName, string fileOldName, string cut, string formula, string newVarName, vector<string> varsToSwitchOn, string tOldName)
{
   vector<string> formulas(1, formula);
   vector<string> newVarNames(1, newVarName);

   copyTreeWithNewVars(fileNewName, fileOldName, cut, formulas, newVarNames, varsToSwitchOn, tOldName);
}

void copyTreeWithNewVars(string fileNewName, string fileOldName, string cut, vector<string> formulas, vector<string> newVarNames, vector<string> varsToSwitchOn, string tOldName)
{
   TFile fOld(fileOldName.c_str());
   TTree* tOld = (TTree*)fOld.Get(tOldName.c_str());

   if(!tOld)
   {
      cerr<<"ERROR: in copyTreeWithNewVars, tOld not found"<<endl;
      return;
   }

   if(varsToSwitchOn.size() > 0)
   {
      tOld->SetBranchStatus("*", 0);
      for(unsigned int i(0); i<varsToSwitchOn.size(); ++i) setBranchStatusTTF(tOld, varsToSwitchOn[i]);
   }

   TFile fNew(fileNewName.c_str(), "RECREATE");
   TTree* tNew(NULL);
   
   copyTreeWithNewVars(tNew, tOld, cut, formulas, newVarNames);

   fNew.cd();
   tNew->Write("", TObject::kOverwrite);

   fNew.Close();
   fOld.Close();

   cout<<"Tree saved in file "<<fileNewName<<endl;
}

void copyTreeWithNewVars(TTree*& tNew, TTree* tOld, string cut, vector<string> formulas, vector<string> newVarNames)
{
   if(formulas.size() != newVarNames.size())
   {
      cerr<<"ERROR: in copyTreeWithNewVar, mismatch. Abort."<<endl;
      return;
   }

   int nVars(formulas.size());

   setBranchStatusTTF(tOld, cut);
   for(int i(0); i<nVars; ++i) setBranchStatusTTF(tOld, formulas.at(i));

   if(cut == "") cut = "1";
   TTreeFormula cutTTF("cutTTF", cut.c_str(), tOld);

   vector<TTreeFormula*> varTTFs(nVars, NULL);
   for(int i(0); i<nVars; ++i) varTTFs[i] = new TTreeFormula(("varTTF"+i2s(i)).c_str(), formulas.at(i).c_str(), tOld);

   cout<<endl<<"======================================================================="<<endl;
   cout<<"Errors from TTreeFormula do not matter from here"<<endl;
   for(int i(0); i<nVars; ++i) tOld->SetBranchStatus(newVarNames.at(i).c_str(), 0);
   cout<<"From now on errors matter"<<endl;
   cout<<"======================================================================="<<endl<<endl;


   tNew = tOld->CloneTree(0);

   vector<double> newVarVals(nVars, 0.);
   for(int i(0); i<nVars; ++i) tNew->Branch(newVarNames.at(i).c_str(), &newVarVals.at(i), (newVarNames.at(i)+"/D").c_str() );

   cout<<"Filling tree with new var... ";

   for(int i(0); i<tOld->GetEntries(); ++i)
   {
      tOld->GetEntry(i);

      if(cutTTF.EvalInstance() )
      {
         for(unsigned int i(0); i<varTTFs.size(); ++i) newVarVals.at(i) = varTTFs.at(i)->EvalInstance();
         tNew->Fill();
      }
   }

   cout<<"done! New tree has: "<<tNew->GetEntries()<<" entries."<<endl;

   for(int i(0); i<nVars; ++i)
   {
      if(varTTFs.at(i) != NULL) delete varTTFs.at(i);
   }
}



void reduceTreeByFactor(unsigned int factor, string nameFile, string nameReducedFile, string nameTree)
{
   TFile f(nameFile.c_str());
   TTree* t = (TTree*)f.Get(nameTree.c_str());

   if(!t)
   {
      cerr<<"ERROR in reduceTreeByFactor: no tree found, do nothing"<<endl;
   }

   TRandom rand;

   TFile f2(nameReducedFile.c_str(), "RECREATE");
   TTree* t2 = t->CloneTree(0);
   for(int i(0); i<t->GetEntries(); ++i)
   {
      t->GetEntry(i);
      if(rand.Integer(factor) == 0) t2->Fill();
   }

   t2->Write("", TObject::kOverwrite);

   f2.Close();
   f.Close();
}

void roundErrorIn2DHistogram(TH2* h)
{
   int nBinsX(h->GetXaxis()->GetNbins());
   int nBinsY(h->GetYaxis()->GetNbins());

   ValError stock;

   for(int i(1); i<=nBinsX; ++i)
   {
      for(int j(1); j<=nBinsY; ++j)
      {
         stock.val = h->GetBinContent(i,j);
         stock.err = h->GetBinError(i,j);

         stock.roundToError();

         h->SetBinContent(i,j,stock.val);
         h->SetBinError(i,j,stock.err);
      }
   }
}

int getKFoldFromNumbers(unsigned long long int eventNumber, int nFolds)
{
   return eventNumber % nFolds;
}


void savePullPlot(RooPlot& graph, string fileName)
{
   RooHist* hist = graph.pullHist(0,0,true);
   TGraphAsymmErrors tgae(hist->GetN());
   tgae.SetTitle("");
   double x(0);
   double y(0);
   for(int i(0); i<hist->GetN(); ++i)
   {
      hist->GetPoint(i,x,y);
      tgae.SetPoint(i,x,y); 
      tgae.SetPointError(i, hist->GetErrorXlow(i), hist->GetErrorXhigh(i),
            hist->GetErrorYlow(i), hist->GetErrorYhigh(i));
   }

   TCanvas c_pullplot("pullplot", "pullplot", 600, 150);
   tgae.SetMarkerStyle(20);
   tgae.SetMarkerSize(1);
   tgae.SetMarkerColor(1);

   tgae.GetXaxis()->SetLimits(graph.GetXaxis()->GetXmin(), graph.GetXaxis()->GetXmax());
   tgae.GetYaxis()->SetRangeUser(-5,5);
   tgae.GetYaxis()->SetLabelSize(0.);
   tgae.GetXaxis()->SetNdivisions(0);
   tgae.GetYaxis()->SetNdivisions(0);//(503);
   tgae.GetYaxis()->SetLabelSize(0);//(0.133);
   tgae.GetXaxis()->SetLabelSize(0);//(0.133);

   TLine line1(graph.GetXaxis()->GetXmin(), -3, graph.GetXaxis()->GetXmax(), -3);
   line1.SetLineColor(2);
   TLine line2(graph.GetXaxis()->GetXmin(), 3, graph.GetXaxis()->GetXmax(), 3);
   line2.SetLineColor(2);
   TLine line3(graph.GetXaxis()->GetXmin(), 0, graph.GetXaxis()->GetXmax(), 0);
   line3.SetLineColor(1);

   tgae.Draw("AP");
   line1.Draw();
   line2.Draw();
   line3.Draw();

   TFile f(fileName.c_str(), "RECREATE");

   c_pullplot.Write();

   f.Close();
}


void saveFitInfo(ofstream& out, RooPlot* frame, int floatingParams, RooAbsPdf* pdf)
{
   string separation("====================================");
   if(frame)
   {
      double chi2;
      int ndof;
      out<<separation<<endl;
      ComputeChi2 computeChi2;
      computeChi2.getChi2(frame, floatingParams, chi2, ndof);
      out<<"Chi2 = "<<chi2<<" / "<<ndof-floatingParams<<endl;
      out<<"Prob = "<<TMath::Prob(chi2, ndof-floatingParams)<<endl;

   }

   if(pdf)
   {
      out<<separation<<endl;
      RooArgSet* vars = pdf->getVariables();
      vars->writeToStream(out, false);
   }

}

void superImposeEfficiencyPlots(vector<string> fileNamesPass, vector<string> histNamesPass, vector<string> fileNamesTot, vector<string> histNamesTot, 
   string plotName, string yLabel, string xLabel, double min, double max, string plotOption, vector<string> legEntries, vector<double> legPos)
{
   unsigned int nTot(fileNamesPass.size());

   if(nTot == 0)
   {
      cerr<<"WARNING: in superImposeEfficiencyPlots, nothing to plot!"<<endl;
      return;
   }

   if(histNamesPass.size() != nTot || fileNamesTot.size() != nTot || histNamesTot.size() != nTot)
   {
      cerr<<"ERROR: in superImposeEfficiencyPlots, mismatch in the number of histograms. "<<endl;
      return;
   }

   vector<TFile*> filesPass(nTot, NULL);
   vector<TFile*> filesTot(nTot, NULL);
   vector<TH1*> histPass(nTot, NULL);
   vector<TH1*> histTot(nTot, NULL);

   for(unsigned int i(0); i < nTot; ++i)
   {
      filesPass.at(i) = new TFile( (fileNamesPass.at(i)).c_str() );
      histPass.at(i) = (TH1*)filesPass.at(i)->Get( (histNamesPass.at(i)).c_str() );
      filesTot.at(i) = new TFile( (fileNamesTot.at(i)).c_str() );
      histTot.at(i) = (TH1*)filesTot.at(i)->Get( (histNamesTot.at(i)).c_str() );
   }

   superImposeEfficiencyPlots(histPass, histTot, plotName, yLabel, xLabel, min, max, plotOption, legEntries, legPos);

   for(unsigned int i(0); i < nTot; ++i)
   {
      if(filesPass.at(i) != NULL)
      {
         filesPass.at(i)->Close();
         delete filesPass.at(i);
      }

      if(filesTot.at(i) != NULL)
      {
         filesTot.at(i)->Close();
         delete filesTot.at(i);
      }
   }
}



void makeStdDevPlot(vector<string> fileNamesPass, vector<string> histNamesPass, vector<string> fileNamesTot, vector<string> histNamesTot, 
   string plotName, string yLabel, string xLabel, double min, double max, string plotOption)
{
   unsigned int nTot(fileNamesPass.size());

   if(nTot == 0)
   {
      cerr<<"WARNING: in superImposeEfficiencyPlots, nothing to plot!"<<endl;
      return;
   }

   if(histNamesPass.size() != nTot || fileNamesTot.size() != nTot || histNamesTot.size() != nTot)
   {
      cerr<<"ERROR: in superImposeEfficiencyPlots, mismatch in the number of histograms. "<<endl;
      return;
   }

   vector<TFile*> filesPass(nTot, NULL);
   vector<TFile*> filesTot(nTot, NULL);
   vector<TH1*> histPass(nTot, NULL);
   vector<TH1*> histTot(nTot, NULL);

   for(unsigned int i(0); i < nTot; ++i)
   {
      filesPass.at(i) = new TFile( (fileNamesPass.at(i)).c_str() );
      histPass.at(i) = (TH1*)filesPass.at(i)->Get( (histNamesPass.at(i)).c_str() );
      filesTot.at(i) = new TFile( (fileNamesTot.at(i)).c_str() );
      histTot.at(i) = (TH1*)filesTot.at(i)->Get( (histNamesTot.at(i)).c_str() );
   }

   makeStdDevPlot(histPass, histTot, plotName, yLabel, xLabel, min, max, plotOption);

   for(unsigned int i(0); i < nTot; ++i)
   {
      if(filesPass.at(i) != NULL)
      {
         filesPass.at(i)->Close();
         delete filesPass.at(i);
      }

      if(filesTot.at(i) != NULL)
      {
         filesTot.at(i)->Close();
         delete filesTot.at(i);
      }
   }
}



void makeStdDevPlot(vector<TH1*> passHists, vector<TH1*> totHists, string plotName, string yLabel, 
   string xLabel, double min, double max, string plotOption)
{
   unsigned int nPlots(passHists.size());

   if(nPlots == 0)
   {
      cout<<"WARNING: in makeStdDevPlot: nothing to print!"<<endl;
      return;
   }

   if(totHists.size() != nPlots)
   {
      cerr<<"ERROR: in makeStdDevPlot, mismatch in the number of plots"<<endl;
      return;
   }

   //check if the histograms are fine

   for(unsigned int i(0); i < nPlots; ++i)
   {
      if(passHists.at(i) == NULL || totHists.at(i) == NULL)
      {
         cerr<<"ERROR: in makeStdDevPlot, "<<i<<"th histogram is NULL "<<passHists.at(i)<<" "<<totHists.at(i)<<endl;
         return;
      }
   }

   //make fraction

   vector<TH1*> effHists(nPlots, NULL); 

   for(unsigned int i(0); i < nPlots; ++i)
   {
      effHists.at(i) = (TH1*)passHists.at(i)->Clone(("effHist"+i2s(i)).c_str()); 
      effHists.at(i)->Divide(passHists.at(i), totHists.at(i), 1., 1., "B");
   }

   int nBins( effHists.at(0)->GetNbinsX() );

   TH1F* h = (TH1F*)effHists.at(0)->Clone();

   vector<double> means(nPlots, 0);
   ValError binVal(0,0);

   for(int i(0); i<nBins; ++i)
   {
      for(unsigned int p(0); p<nPlots; ++p) means.at(p) = effHists.at(p)->GetBinContent(i+1);

      binVal = getMeanWithStdDev(means);
      cout<<binVal<<endl;
      h->SetBinContent(i+1, binVal.val);
      h->SetBinError(i+1, binVal.err);
   }

   //beautify and print 

   if(min > -10. || max > -10. ) h->GetYaxis()->SetRangeUser(min, max);
   if(xLabel != "") h->GetXaxis()->SetTitle(xLabel.c_str());
   if(yLabel != "") h->GetYaxis()->SetTitle(yLabel.c_str());

   TCanvas canv("canv", "canv", 600, 600);
   h->Draw(plotOption.c_str());
   canv.Print( plotName.c_str() );

}

void superImposeEfficiencyPlots(vector<TH1*> passHists, vector<TH1*> totHists, string plotName, string yLabel, 
      string xLabel, double min, double max, string plotOption, vector<string> legEntries, vector<double> legPos)
{
   unsigned int nPlots(passHists.size());

   if(totHists.size() != nPlots)
   {
      cerr<<"ERROR: in superImposeEfficiencyPlots, mismatch in the number of plots"<<endl;
      return;
   }

   //check if the histograms are fine

   for(unsigned int i(0); i < nPlots; ++i)
   {
      if(passHists.at(i) == NULL || totHists.at(i) == NULL)
      {
         cerr<<"ERROR: in superImposeEfficiencyPlots, "<<i<<"th histogram is NULL "<<passHists.at(i)<<" "<<totHists.at(i)<<endl;
         return;
      }
   }

   //make fraction

   vector<TH1*> effHists(nPlots, NULL); 

   for(unsigned int i(0); i < nPlots; ++i)
   {
      effHists.at(i) = (TH1*)passHists.at(i)->Clone(("effHist"+i2s(i)).c_str()); 
      effHists.at(i)->Divide(passHists.at(i), totHists.at(i), 1., 1., "B");
   }

   //see if need to set max and min

   if(min < -10. && max < -10. )
   {
      min = numeric_limits<double>::max();
      max = numeric_limits<double>::min();
      double dummy(0);

      //in that case, set max and min yourself
      for(unsigned int i(0); i < nPlots; ++i)
      {
         for(int j(1); j <= effHists[i]->GetNbinsX(); ++j)
         {
            dummy = effHists[i]->GetBinContent(j) + effHists[i]->GetBinError(j);
            if(max < dummy) max = dummy;

            dummy = effHists[i]->GetBinContent(j) - effHists[i]->GetBinError(j);
            if(min > dummy) min = dummy;
         }
      }

      double shift(0.05*(max-min)); 
      max += shift;
      min -= shift;
   }
   
   //beautify  

   double nColor(1);  

   for(unsigned int i(0); i<nPlots; ++i)
   {
      if(nColor == 5) nColor = 6;
      if(nColor == 10) nColor = 30;
      effHists.at(i)->SetLineColor(nColor); 
      effHists.at(i)->SetMarkerColor(nColor);
      effHists.at(i)->GetYaxis()->SetRangeUser(min, max);
      if(xLabel != "") effHists.at(i)->GetXaxis()->SetTitle(xLabel.c_str());
      if(yLabel != "") effHists.at(i)->GetYaxis()->SetTitle(yLabel.c_str());
      ++nColor;
   }



   //plot

   TCanvas canv("canv", "canv", 600, 600);

   TLegend leg(legPos.at(0), legPos.at(1), legPos.at(2), legPos.at(3)); 
   leg.SetLineColor(0);
   leg.SetFillColor(0);

   if(nPlots > 0) effHists.at(0)->Draw(plotOption.c_str() ); 
   for(unsigned int i(1); i<nPlots; ++i) effHists.at(i)->Draw( (plotOption+"same").c_str() );

   if(legEntries.size() == nPlots )
   {
      for(unsigned int i(0); i<nPlots; ++i) leg.AddEntry(effHists.at(i), (legEntries.at(i)).c_str() );
      leg.Draw();
   }


   canv.Print(plotName.c_str() );
}







void projectEfficiencyPlots(vector<TH3*> passHists, vector<TH3*> totHists, string plotNameStart, string yLabel,
   vector<string> xLabel, double min, double max, string plotOption, vector<string> legEntries, vector<double> legPos)
{
   unsigned int nPlots(passHists.size());       

   if(totHists.size() != nPlots) 
   {
      cerr<<"ERROR: in function projectEfficiencyPlots, mismatch in pass and tot hists"<<endl;
      return;
   }

   if( xLabel.size() != 3) xLabel = vector<string>(3, "");

   //check if the histograms are fine


   for(unsigned int i(0); i < nPlots; ++i)
   {
      if(passHists.at(i) == NULL || totHists.at(i) == NULL)
      {
         cerr<<"ERROR: in projectEfficiencyPlots, "<<i<<"th histogram is NULL "<<passHists.at(i)<<" "<<totHists.at(i)<<endl;
         return;
      }
   }

   //make projections

   vector<TH1*> passProjX(nPlots, NULL); 
   vector<TH1*> passProjY(nPlots, NULL); 
   vector<TH1*> passProjZ(nPlots, NULL); 

   vector<TH1*> totProjX(nPlots, NULL); 
   vector<TH1*> totProjY(nPlots, NULL); 
   vector<TH1*> totProjZ(nPlots, NULL); 

   int nBinsX;
   int nBinsY;
   int nBinsZ;

   for(unsigned int i(0); i < nPlots; ++i)
   {
      nBinsX = passHists.at(i) -> GetNbinsX(); 
      nBinsY = passHists.at(i) -> GetNbinsY(); 
      nBinsZ = passHists.at(i) -> GetNbinsZ(); 

      passProjX.at(i) = passHists.at(i)->ProjectionX( ("passProjX"+i2s(i)).c_str(), 1, nBinsY, 1, nBinsZ );
      passProjY.at(i) = passHists.at(i)->ProjectionY( ("passProjY"+i2s(i)).c_str(), 1, nBinsX, 1, nBinsZ );
      passProjZ.at(i) = passHists.at(i)->ProjectionZ( ("passProjZ"+i2s(i)).c_str(), 1, nBinsX, 1, nBinsY );

      totProjX.at(i) = totHists.at(i)->ProjectionX( ("totProjX"+i2s(i)).c_str(), 1, nBinsY, 1, nBinsZ );
      totProjY.at(i) = totHists.at(i)->ProjectionY( ("totProjY"+i2s(i)).c_str(), 1, nBinsX, 1, nBinsZ );
      totProjZ.at(i) = totHists.at(i)->ProjectionZ( ("totProjZ"+i2s(i)).c_str(), 1, nBinsX, 1, nBinsY );
   }

   //plot the projections

   superImposeEfficiencyPlots(passProjX, totProjX, plotNameStart+"_projX.pdf", yLabel, 
         xLabel.at(0), min, max, plotOption, legEntries, legPos);

   superImposeEfficiencyPlots(passProjY, totProjY, plotNameStart+"_projY.pdf", yLabel, 
         xLabel.at(1), min, max, plotOption, legEntries, legPos);

   superImposeEfficiencyPlots(passProjZ, totProjZ, plotNameStart+"_projZ.pdf", yLabel, 
         xLabel.at(2), min, max, plotOption, legEntries, legPos);
}



void projectEfficiencyPlots(vector<string> fileNames, vector<string> totHistNames, vector<string> passHistNames, string plotNameStart, string yLabel,
   vector<string> xLabel, double min, double max, string plotOption, vector<string> legEntries, vector<double> legPos)
{
   unsigned int nPlots(fileNames.size()); 

   if(totHistNames.size() != nPlots || passHistNames.size() != nPlots)
   {
      cerr<<"ERROR: in function projectEfficiencyPlots, mismatch in the number of histograms"<<endl;
      return; 
   }

   //recover histograms from files

   vector<TFile*> files(nPlots, NULL);
   vector<TH3*> totHists(nPlots, NULL);
   vector<TH3*> passHists(nPlots, NULL);

   for(unsigned int i(0); i < nPlots; ++i)
   {
      files.at(i) = new TFile((fileNames.at(i)).c_str());
      totHists.at(i) = (TH3*)files.at(i)->Get(totHistNames.at(i).c_str());
      passHists.at(i) = (TH3*)files.at(i)->Get(passHistNames.at(i).c_str());
   }

   //plot everything

   projectEfficiencyPlots(passHists, totHists, plotNameStart, yLabel,
         xLabel, min, max, plotOption, legEntries, legPos);

   //close all

   for(unsigned int i(0); i < nPlots; ++i)
   {
      files.at(i)->Close();
      if(files.at(i) != NULL) delete files.at(i);
   }
}




void projectEfficiencyPlots2D(vector<string> fileNames, vector<string> totHistNames, vector<string> passHistNames, string plotNameStart, string yLabel,
   vector<string> xLabel, double min, double max, string plotOption, vector<string> legEntries, vector<double> legPos)
{
   unsigned int nPlots(fileNames.size()); 

   if(totHistNames.size() != nPlots || passHistNames.size() != nPlots)
   {
      cerr<<"ERROR: in function projectEfficiencyPlots2D, mismatch in the number of histograms"<<endl;
      return; 
   }

   //recover histograms from files

   vector<TFile*> files(nPlots, NULL);
   vector<TH2*> totHists(nPlots, NULL);
   vector<TH2*> passHists(nPlots, NULL);

   for(unsigned int i(0); i < nPlots; ++i)
   {
      files.at(i) = new TFile((fileNames.at(i)).c_str());
      totHists.at(i) = (TH2*)files.at(i)->Get(totHistNames.at(i).c_str());
      passHists.at(i) = (TH2*)files.at(i)->Get(passHistNames.at(i).c_str());
   }

   //plot everything

   projectEfficiencyPlots(passHists, totHists, plotNameStart, yLabel,
         xLabel, min, max, plotOption, legEntries, legPos);

   //close all

   for(unsigned int i(0); i < nPlots; ++i)
   {
      files.at(i)->Close();
      if(files.at(i) != NULL) delete files.at(i);
   }
}






void projectEfficiencyPlots(vector<TH2*> passHists, vector<TH2*> totHists, string plotNameStart, string yLabel,
   vector<string> xLabel, double min, double max, string plotOption, vector<string> legEntries, vector<double> legPos)
{
   unsigned int nPlots(passHists.size());       

   if(totHists.size() != nPlots) 
   {
      cerr<<"ERROR: in function projectEfficiencyPlots, mismatch in pass and tot hists"<<endl;
      return;
   }

   if( xLabel.size() != 2) xLabel = vector<string>(2, "");

   //check if the histograms are fine


   for(unsigned int i(0); i < nPlots; ++i)
   {
      if(passHists.at(i) == NULL || totHists.at(i) == NULL)
      {
         cerr<<"ERROR: in projectEfficiencyPlots, "<<i<<"th histogram is NULL "<<passHists.at(i)<<" "<<totHists.at(i)<<endl;
         return;
      }
   }

   //make projections

   vector<TH1*> passProjX(nPlots, NULL); 
   vector<TH1*> passProjY(nPlots, NULL); 

   vector<TH1*> totProjX(nPlots, NULL); 
   vector<TH1*> totProjY(nPlots, NULL); 

   int nBinsX;
   int nBinsY;

   for(unsigned int i(0); i < nPlots; ++i)
   {
      nBinsX = passHists.at(i) -> GetNbinsX(); 
      nBinsY = passHists.at(i) -> GetNbinsY(); 

      passProjX.at(i) = passHists.at(i)->ProjectionX( ("passProjX"+i2s(i)).c_str(), 1, nBinsY);
      passProjY.at(i) = passHists.at(i)->ProjectionY( ("passProjY"+i2s(i)).c_str(), 1, nBinsX);

      totProjX.at(i) = totHists.at(i)->ProjectionX( ("totProjX"+i2s(i)).c_str(), 1, nBinsY );
      totProjY.at(i) = totHists.at(i)->ProjectionY( ("totProjY"+i2s(i)).c_str(), 1, nBinsX );
   }

   //plot the projections

   superImposeEfficiencyPlots(passProjX, totProjX, plotNameStart+"_projX.pdf", yLabel, 
         xLabel.at(0), min, max, plotOption, legEntries, legPos);

   superImposeEfficiencyPlots(passProjY, totProjY, plotNameStart+"_projY.pdf", yLabel, 
         xLabel.at(1), min, max, plotOption, legEntries, legPos);
}

void gaussOscillateHisto(TH3* hist)
{
   TRandom rand;

   double val(0); 
   double err(0);
   double newVal(0);

   for(int i(1); i<=hist->GetNbinsX(); ++i)
   {
      for(int j(1); j<=hist->GetNbinsX(); ++j)
      {
         for(int k(1); k<=hist->GetNbinsX(); ++k)
         {
            val = hist->GetBinContent(i,j,k);
            err = hist->GetBinError(i,j,k);

            newVal = rand.Gaus(val, err); 

            hist->SetBinContent(i, j, k, newVal);
         }
      }
   }
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

void fillTreeResultSimple(TTree* t, RooArgSet* par_set, RooFitResult* rfr, bool update, int migradRes, int hesseRes,bool  hasConverged)
{


  t->ResetBranchAddresses();

  int nVars(par_set->getSize());

  vector<double> vec(nVars);
  vector<double> errs(nVars);

  string name;
  double value, error;

  TIterator *pariter = par_set->createIterator();
  RooRealVar *pari = (RooRealVar*) pariter->Next();

  int i=0;

  while (pari)
  {
    name = pari->GetName();
    value = pari->getValV();
    error = 0;

    if(!(pari->isConstant())) error = pari->getError();

    vec[i] = value;
    errs[i] = error;

    if(!update) t->Branch( name.c_str(), &vec[i], (name+"/D").c_str() );
    if(update) t->SetBranchAddress(name.c_str(), &vec[i]);

    if(!update) t->Branch( (name+"_err").c_str(), &errs[i], (name+"_err/D").c_str() );
    if(update) t->SetBranchAddress((name+"_err").c_str(), &errs[i]);

    pari = (RooRealVar*) pariter->Next();
    i++;

  }



  double minNll(rfr->minNll());


  if(!update) t->Branch( "migradRes", &migradRes, "migradRes/I" );
  if(update) t->SetBranchAddress("migradRes", &migradRes);

  if(!update) t->Branch( "hesseRes", &hesseRes, "hesseRes/I" );
  if(update) t->SetBranchAddress("hesseRes", &migradRes);

  if(!update) t->Branch( "minNll", &minNll, "minNll/D" );
  if(update) t->SetBranchAddress("minNll", &minNll);

  if(!update) t->Branch( "hasConverged", &hasConverged, "hasConverged/O" );
  if(update) t->SetBranchAddress("hasConverged", &hasConverged);


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



void makeTableResults(string filename, string treename, int nGenSignal, int nGenPartReco, int nGenComb, int nGenJpsiLeak, ostream& out, bool wantLatex )
{
  TFile f(filename.c_str());
  TTree* t = (TTree*)f.Get(treename.c_str());

  makeTableResults(t, nGenSignal, nGenPartReco,  nGenComb, nGenJpsiLeak, out, wantLatex);

  f.Close();
}



void makeTableResults(TTree* t, int nGenSignal, int nGenPartReco, int nGenComb, int nGenJpsiLeak, ostream& out, bool wantLatex)//, string titleString )
{

  gStyle->SetOptFit(1);

  ValError ve;
  TH1F* hStock(0);
  TCanvas canv("canv", "canv", 1800, 600);
  canv.Divide(4,1);
  TCanvas canv_pull("canv_pull", "canv_pull", 1800, 600);
  canv_pull.Divide(4,1);
  TVirtualPad *pad;

  TCanvas canv_sig("canv_sig", "canv_sig",600,600);
  TCanvas canv_sig_pull("canv_sig_pull","canv_sig_pull",600,600);


  std::string cut = "hasConverged>0";

  pad = canv.cd(1);
  t->Draw("nSignal>>hStock(50)",cut.c_str(),"e1");
  cout<<"nsignal: "<<t->GetEntries("nSignal>0")<<endl;
  // hStock = (TH1F*)canv.GetPrimitive("hStock");
  hStock = (TH1F*)pad->GetPrimitive("hStock");
  hStock->Fit("gaus");
  cout<<"Mean is: "<<hStock->GetMean()<<endl;
  canv_sig.cd();
  hStock->GetXaxis()->SetTitle("nSignal");
  hStock->Draw();
  canv_sig.SaveAs("signal_distribution.pdf");

  // ve.val = hStock->GetMean();
  // ve.err = hStock->GetStdDev();
  ve.val = hStock->GetFunction("gaus")->GetParameter(1);
  ve.err = hStock->GetFunction("gaus")->GetParameter(2);
  double fracSig(100*ve.err/ve.val);
  string strSig(roundToError(ve, wantLatex));

  pad = canv_pull.cd(1);
  t->Draw(("((nSignal-"+i2s(nGenSignal)+")/nSignal_err)>>pullsignal").c_str(),cut.c_str(),"e1");
  hStock = (TH1F*)pad->GetPrimitive("pullsignal");
  hStock->Fit("gaus");
  canv_sig_pull.cd();
  hStock->GetXaxis()->SetTitle("nSignal pull");
  hStock->Draw();
  canv_sig_pull.SaveAs("signal_pull.pdf");
  ve.val = hStock->GetFunction("gaus")->GetParameter(1);
  ve.err = hStock->GetFunction("gaus")->GetParError(1);
  string strpullmeanSig(roundToError(ve, wantLatex));
  ve.val = hStock->GetFunction("gaus")->GetParameter(2);
  ve.err = hStock->GetFunction("gaus")->GetParError(2);
  string strpullsigmaSig(roundToError(ve, wantLatex));


  pad = canv.cd(2);
  t->Draw("nComb>>h2",cut.c_str(),"e1");
  // hStock = (TH1F*)canv.GetPrimitive("h2");
  hStock = (TH1F*)pad->GetPrimitive("h2");
  hStock->Fit("gaus");
  // ve.val = hStock->GetMean();
  // ve.err = hStock->GetStdDev();
  ve.val = hStock->GetFunction("gaus")->GetParameter(1);
  ve.err = hStock->GetFunction("gaus")->GetParameter(2);
  double fracBkg(100*ve.err/ve.val);
  string strBkg = roundToError(ve, wantLatex);
  pad = canv_pull.cd(2);
  t->Draw(("((nComb-"+i2s(nGenComb)+")/nComb_err)>>pullcomb").c_str(),cut.c_str(),"e1");
  hStock = (TH1F*)pad->GetPrimitive("pullcomb");
  hStock->Fit("gaus");
  ve.val = hStock->GetFunction("gaus")->GetParameter(1);
  ve.err = hStock->GetFunction("gaus")->GetParError(1);
  string strpullmeanBkg(roundToError(ve, wantLatex));
  ve.val = hStock->GetFunction("gaus")->GetParameter(2);
  ve.err = hStock->GetFunction("gaus")->GetParError(2);
  string strpullsigmaBkg(roundToError(ve, wantLatex));


  pad = canv.cd(3);
  t->Draw("nPartReco>>h3",cut.c_str(),"e1");
  // hStock = (TH1F*)canv.GetPrimitive("h3");
  hStock = (TH1F*)pad->GetPrimitive("h3");
  hStock->Fit("gaus");
  // ve.val = hStock->GetMean();
  // ve.err = hStock->GetStdDev();
  ve.val = hStock->GetFunction("gaus")->GetParameter(1);
  ve.err = hStock->GetFunction("gaus")->GetParameter(2);
  double fracPartReco(100*ve.err/ve.val);
  string strPartReco = roundToError(ve, wantLatex);
  pad = canv_pull.cd(3);
  t->Draw(("((nPartReco-"+i2s(nGenPartReco)+")/nPartReco_err)>>pullpartreco").c_str(),cut.c_str(),"e1");
  hStock = (TH1F*)pad->GetPrimitive("pullpartreco");
  hStock->Fit("gaus");
  ve.val = hStock->GetFunction("gaus")->GetParameter(1);
  ve.err = hStock->GetFunction("gaus")->GetParError(1);
  string strpullmeanPartReco(roundToError(ve, wantLatex));
  ve.val = hStock->GetFunction("gaus")->GetParameter(2);
  ve.err = hStock->GetFunction("gaus")->GetParError(2);
  string strpullsigmaPartReco(roundToError(ve, wantLatex));

  double fracJpsiLeak(0);
  string strJpsiLeak("0");
  string strpullmeanJpsiLeak("0");
  string strpullsigmaJpsiLeak("0");
  
  if(nGenJpsiLeak > 0)
  {
    pad = canv.cd(4);
    t->Draw("nJpsiLeak>>h4",cut.c_str(),"e1");
    // hStock = (TH1F*)canv.GetPrimitive("h4");
    hStock = (TH1F*)pad->GetPrimitive("h4");
    hStock->Fit("gaus");
    // ve.val = hStock->GetMean();
    // ve.err = hStock->GetStdDev();
    ve.val = hStock->GetFunction("gaus")->GetParameter(1);
    ve.err = hStock->GetFunction("gaus")->GetParameter(2);
    fracJpsiLeak = 100*ve.err/ve.val;
    strJpsiLeak = roundToError(ve, wantLatex);
    pad = canv_pull.cd(4);
    t->Draw(("((nJpsiLeak-"+i2s(nGenJpsiLeak)+")/nJpsiLeak_err)>>pulljpsi").c_str(),cut.c_str(),"e1");
    hStock = (TH1F*)pad->GetPrimitive("pulljpsi");
    hStock->Fit("gaus");
    ve.val = hStock->GetFunction("gaus")->GetParameter(1);
    ve.err = hStock->GetFunction("gaus")->GetParError(1);
    strpullmeanJpsiLeak = roundToError(ve, wantLatex);
    ve.val = hStock->GetFunction("gaus")->GetParameter(2);
    ve.err = hStock->GetFunction("gaus")->GetParError(2);
    strpullsigmaJpsiLeak = roundToError(ve, wantLatex);

  }

  canv.SaveAs("par_distribution.pdf");
  canv_pull.SaveAs("par_pull.pdf");


  int w(16);
  if(!wantLatex)
  {
    out<<setw(w)<<" "<<setw(w)<<"nSignal"<<","<<setw(w)<<"nPartReco"<<","<<setw(w)<<"nCombinatorial"<<","<<setw(w)<<"nJpsiLeak"<<endl;
    out<<setw(w)<<"Generated"<<","<<setw(w)<<nGenSignal<<","<<setw(w)<<nGenPartReco<<","<<setw(w)<<nGenComb<<","<<setw(w)<<nGenJpsiLeak<<endl;
    out<<setw(w)<<"Fit"<<","<<setw(w)<<strSig<<","<<setw(w)<<strPartReco<<","<<setw(w)<<strBkg<<","<<setw(w)<<strJpsiLeak<<endl;
    out<<setw(w)<<"FracErr"<<","<<setw(w)<<fracSig<<","<<setw(w)<<fracPartReco<<","<<setw(w)<<fracBkg<<","<<setw(w)<<fracJpsiLeak<<endl;
    out<<setw(w)<<"PullMean"<<","<<setw(w)<<strpullmeanSig<<","<<setw(w)<<strpullmeanPartReco<<","<<setw(w)<<strpullmeanBkg<<","<<setw(w)<<strpullmeanJpsiLeak<<endl;
    out<<setw(w)<<"PullSigma"<<","<<setw(w)<<strpullsigmaSig<<","<<setw(w)<<strpullsigmaPartReco<<","<<setw(w)<<strpullsigmaBkg<<","<<setw(w)<<strpullsigmaJpsiLeak<<endl;

  
}

  if(wantLatex)
  {
    out<<"\\documentclass[11pt]{article}"<<endl;
    out<<"\\usepackage[left=0cm,right=0cm,top=0cm,bottom=0cm,paperwidth=15cm,paperheight=5cm]{geometry}"<<endl;
    out<<"\\usepackage{xcolor,colortbl}"<<endl;
    out<<"\\begin{document}"<<endl;
    out<<"\\begin{tabular}{r r r r r}"<<endl;
    out<<"  &    \\textbf{Signal}    &   \\textbf{PartReco} & \\textbf{Comb.}  & \\textbf{\\boldmath$J/\\psi$ leak}\\\\"<<endl;
    out<<"\\cellcolor{gray} &  \\multicolumn{4}{c}{\\cellcolor{gray} 1D $M$ fit:} \\\\"<<endl;
    out<<"Generated & $"<<nGenSignal<<"$ & $"<<nGenPartReco<<"$ & $"<<nGenComb<<"$ & $"<<nGenJpsiLeak<<"$ \\\\"<<endl;
    out<<"From fit & "<<strSig<<" & "<<strPartReco<<" & "<<strBkg<<" & "<<strJpsiLeak<<" \\\\"<<endl;
    out<<"Error (\\%) & \\textcolor{red}{$"<<setprecision(3)<<fracSig<<"$} & $"<<setprecision(3)<<fracPartReco<<"$ & $"<<setprecision(3)<<fracBkg<<"$ & $"<<setprecision(3)<<fracJpsiLeak<<"$ \\\\"<<endl;
    out<<"Bias & "<<strpullmeanSig<<" & "<<strpullmeanPartReco<<" & "<<strpullmeanBkg<<" & "<<strpullmeanJpsiLeak<<" \\\\"<<endl;
    out<<"Cover & "<<strpullsigmaSig<<" & "<<strpullsigmaPartReco<<" & "<<strpullsigmaBkg<<" & "<<strpullsigmaJpsiLeak<<" \\\\"<<endl;
    out<<"\\end{tabular}"<<endl;
    out<<"\\end{document}"<<endl;
  }
  cout<<"CACA "<<endl;

}

string roundToError(RooRealVar const& var, bool wantLatex)
{
  ValError ve;
  ve.val = var.getVal();
  ve.err = var.getError();
  return roundToError(ve, wantLatex);
}

string roundToError(ValError& ve, bool wantLatex)
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
  string ret(dbl2str(ve.val, nfixed)+"+-"+dbl2str(ve.err, nfixed));
  if(wantLatex) ret= "$"+dbl2str(ve.val, nfixed)+"\\pm"+dbl2str(ve.err, nfixed)+"$";
  return ret;
}

