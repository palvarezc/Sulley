// $Id: $
#ifndef GENERIC_COMPUTECHI2_HPP 
#define GENERIC_COMPUTECHI2_HPP 
 
// Include files
 
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"


class RooPlot;
 
/** @class ComputeChi2 ComputeChi2.hpp Generic/ComputeChi2.hpp
 * Class that computes a \f$ \chi^2 \f$ from a RooPlot containing
 * one dataset and one fit.
 *
 *  @author Johan Luisier
 *  @date   2011-04-07
 */
class ComputeChi2
{
public: 
  /// Standard constructor
  ComputeChi2( ); 
 
  virtual ~ComputeChi2( ); ///< Destructor
 
  double operator()( RooPlot * dataFrame, const int& fixedParam,
                     double& chi2Value, int& nDof );
  double operator()( RooPlot * dataFrame, const int& fixedParam );
 
  double getChi2( RooPlot * dataFrame, const int& fixedParam,
                  double& chi2Value, int& nDof );
  double getChi2( RooPlot * dataFrame, const int& fixedParam );

  void getChi22D(RooDataSet const& data, RooAbsPdf const& pdf, RooRealVar const& x, RooRealVar const& y, double& chi2Value, int& nDof, int xBins = 20, int yBins = 20);


protected:
 
private:
 
};



#endif // GENERIC_COMPUTECHI2_HPP
 
