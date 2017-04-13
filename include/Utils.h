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

enum class FitType : int {RooFit1D, RooFit2D_RooPTMVis, RooFit2D, HistFact1D, HistFact2D};

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
  
    FitType fittype;
    bool fit2D;
    bool templatestat;
    int interpolate;
    
    string trigStr;
    string weightStr; 
    string BDTVar;
    double BDTcut;
    
    string signalfile, signaltree;
    string partrecofile, partrecotree;
    string JpsiLeakfile, JpsiLeaktree;
    string combfile, combtree;

    string name;

    string plotsfile;
    
    double minBMass, maxBMass;

  
    Options(string arg=""):
      name(arg),
      nGenSignal(0), 
      nGenPartReco(0),
      nGenComb(0), 
      nGenJpsiLeak(0),
      nGenKemu(0),
      nGenFracZeroGamma(0.), 
      nGenFracOneGamma(0.),
      minBMass(4880.),
      maxBMass(5700.),
      trigStr(""),
      weightStr(""),
      BDTVar(""),
      BDTcut(0.),
      signalfile(""),
      partrecofile(""),
      JpsiLeakfile(""),
      combfile(""),
      signaltree("DecayTree"),
      partrecotree("DecayTree"),
      JpsiLeaktree("DecayTree"),
      combtree("DecayTree"),
      fittype(FitType::RooFit1D), 
      fit2D(0),
      interpolate(2),
      templatestat(0),
      plotsfile("plots.root")
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
      minBMass(other.minBMass),
      maxBMass(other.maxBMass),
      trigStr(other.trigStr),
      weightStr(other.weightStr),
      BDTVar(other.BDTVar),
      BDTcut(other.BDTcut),
      signalfile(other.signalfile),
      partrecofile(other.partrecofile),
      JpsiLeakfile(other.JpsiLeakfile),
      combfile(other.combfile),
      signaltree(other.signaltree),
      partrecotree(other.partrecotree),
      JpsiLeaktree(other.JpsiLeaktree),
      combtree(other.combtree),
      fittype(other.fittype), 
      fit2D(other.fit2D),
      interpolate(other.interpolate),
      templatestat(other.templatestat),
      plotsfile(other.plotsfile)
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
       >> this->interpolate
       >> this->templatestat;
    

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
        << "\t" << this->interpolate
        << "\t" << this->templatestat;
    


    return out;
    

  }

  };
  

}


#endif // UTILS_H
