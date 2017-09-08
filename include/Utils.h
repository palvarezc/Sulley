#ifndef UTILS_H 
#define UTILS_H 

// Include files
#include <vector>
#include <string>
#include <map>
#include <ostream>
#include <istream>
#include <TMatrixD.h>
#include <TMath.h>
#include <TObject.h>

/* enum class FitType : int {RooFit1D, RooFit2D_RooPTMVis, RooFit2D, HistFact1D, HistFact2D}; */
enum FitType{RooFit1D, RooFit2D_RooPTMVis, RooFit2D, HistFact1D, HistFact2D};

/* 
 *  @author Paula Alvarez Cartelle
 *  @date   2017-03-09
 *
 *  Options of the fit/toyMC
 */

using namespace std;
namespace Sulley {

  class Options 
  {
    

  public:

    int nGenSignal;  
    int nGenPartReco;
    int nGenComb;
    int nGenJpsiLeak;
    int nGenKemu;
    
    int nGenSignalZeroGamma;
    int nGenSignalOneGamma;
    int nGenSignalTwoGamma;

    double nGenFracZeroGamma;
    double nGenFracOneGamma;

    int nGenSignal_mumu;
    int nGenBkg_mumu;
  
    FitType fittype;
    bool fitSimultaneous;
    bool fit2D;
    bool templatestat;
    int interpolate;
    bool constPartReco;
    bool massshift;

    string trigStr;
    string weightStr; 
    string BDTVar;
    double BDTcut;
    
    string signalfile_mumu, signaltree_mumu;
    string signalfile, signaltree;
    string partrecofile, partrecotree;
    string JpsiLeakfile, JpsiLeaktree;
    string combfile, combtree;

    string name;
    string workspacefile;
    string workspacename;

    string plotsfile;
    string plotsdir;
    
    double minBMass, maxBMass;
    double minBMass_mumu, maxBMass_mumu;
    double minmisPT, maxmisPT;

    int Run;
    int trigCat;
    
    
    Options(string arg=""):
      name(arg),
      nGenSignal(0), 
      nGenPartReco(0),
      nGenComb(0), 
      nGenJpsiLeak(0),
      nGenKemu(0),
      nGenFracZeroGamma(0.), 
      nGenFracOneGamma(0.),
      nGenSignal_mumu(0.),
      nGenBkg_mumu(0.),
      minBMass(4880.),
      maxBMass(5700.),
      minmisPT(0.),
      maxmisPT(5000.),
      minBMass_mumu(5180.),
      maxBMass_mumu(5700.),
      trigStr(""),
      weightStr(""),
      BDTVar(""),
      BDTcut(0.),
      signalfile_mumu(""),
      signalfile(""),
      partrecofile(""),
      JpsiLeakfile(""),
      combfile(""),
      signaltree_mumu("DecayTree"),
      signaltree("DecayTree"),
      partrecotree("DecayTree"),
      JpsiLeaktree("DecayTree"),
      combtree("DecayTree"),
      fittype(RooFit1D), 
      fit2D(0),
      fitSimultaneous(0),
      interpolate(0),
      templatestat(0),
      constPartReco(0),    
      massshift(0),    
      workspacefile("myworkspace.root"),
      workspacename("myworkspace"),
	plotsfile("plots.root"),
      plotsdir(""),
      Run(1),
      trigCat(-1)
  {

    nGenSignalZeroGamma = floor(nGenFracZeroGamma*nGenSignal);
    nGenSignalOneGamma = floor(nGenFracOneGamma*nGenSignal);
    nGenSignalTwoGamma = floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma);


  }


    Options(const Options &other):
      name(""),
      nGenSignal(other.nGenSignal), 
      nGenPartReco(other.nGenPartReco),
      nGenComb(other.nGenComb), 
      nGenJpsiLeak(other.nGenJpsiLeak),
      nGenKemu(other.nGenKemu),
      nGenFracZeroGamma(other.nGenFracZeroGamma), 
      nGenFracOneGamma(other.nGenFracOneGamma),
	nGenSignal_mumu(other.nGenSignal_mumu),
	nGenBkg_mumu(other.nGenBkg_mumu),
      minBMass(other.minBMass),
      maxBMass(other.maxBMass),
      minmisPT(other.minmisPT),
      maxmisPT(other.maxmisPT),
      minBMass_mumu(other.minBMass_mumu),
      maxBMass_mumu(other.maxBMass_mumu),
      trigStr(other.trigStr),
      weightStr(other.weightStr),
      BDTVar(other.BDTVar),
      BDTcut(other.BDTcut),
      signalfile_mumu(other.signalfile_mumu),
      signalfile(other.signalfile),
      partrecofile(other.partrecofile),
      JpsiLeakfile(other.JpsiLeakfile),
      combfile(other.combfile),
      signaltree(other.signaltree),
      signaltree_mumu(other.signaltree_mumu),
      partrecotree(other.partrecotree),
      JpsiLeaktree(other.JpsiLeaktree),
      combtree(other.combtree),
      fittype(other.fittype), 
      fit2D(other.fit2D),
      fitSimultaneous(other.fitSimultaneous),
      interpolate(other.interpolate),
      templatestat(other.templatestat),
      constPartReco(other.constPartReco),
      massshift(other.massshift),
	workspacefile(other.workspacefile),
	workspacename(other.workspacename),
	plotsfile(other.plotsfile),
      plotsdir(other.plotsdir),
      Run(other.Run),
      trigCat(other.trigCat)
  {
    nGenSignalZeroGamma = floor(nGenFracZeroGamma*nGenSignal);
    nGenSignalOneGamma = floor(nGenFracOneGamma*nGenSignal);
    nGenSignalTwoGamma = floor(nGenSignal-nGenSignalZeroGamma-nGenSignalOneGamma);
  }
  


  istream &operator>>(std::istream &in)
  {
    
    
    in >> this->nGenSignal
       >> this->nGenPartReco
       >> this->nGenComb
       >> this->nGenJpsiLeak
       >> this->nGenFracZeroGamma
       >> this->nGenFracOneGamma
       >> this->fit2D
       >> this->fitSimultaneous
       >> this->interpolate
       >> this->templatestat
       >> this->signalfile_mumu;
    

    return in;
    

  }
  

  ostream &operator<<(std::ostream &out)
  {
    

    out.precision(15);
    

    out << this->nGenSignal
        << "\t" << this->nGenPartReco
        << "\t" << this->nGenComb
        << "\t" << this->nGenJpsiLeak
        << "\t" << this->nGenFracZeroGamma
        << "\t" << this->nGenFracOneGamma
        << "\t" << this->fit2D
        << "\t" << this->fitSimultaneous
        << "\t" << this->interpolate
        << "\t" << this->templatestat
        << "\t" << this->signalfile_mumu;
    


    return out;
    

  }

  };
  

}


#endif // UTILS_H
