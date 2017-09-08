// $Id: $
// Include files 
 
#include <RooPlot.h>
#include <RooHist.h>
 
#include <TMath.h>

#include "TH2F.h"
 
// local
#include "ComputeChi2.h"
 
using namespace std;
using namespace RooFit;
 
//-----------------------------------------------------------------------------
// Implementation file for class : ComputeChi2
//
// 2011-04-07 : Johan Luisier
//-----------------------------------------------------------------------------
 
//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ComputeChi2::ComputeChi2(  )
{}
 
//=============================================================================
// Destructor
//=============================================================================
ComputeChi2::~ComputeChi2()
{}
 
/**
 * Do the computation of the \f$ \chi^2 \f$ as follows : the residual histogram
 * is extracted from the RooPlot, and then a loop on the histogram bins is made,      
 * where the residual value and its error are extracted. Checks are performed to
 * be sure that the received values are meaningful, and then the \f$ \chi^2 \f$
 * is computed.
 *
 * The \f$ \chi^2 \f$ is obtained as :
 * \f[ \chi^2 = \sum_{i = 0}^n \left( \frac{y_i}{\sigma_i] \right)^2 \f]
 *
 * @param dataFrame pointer on the RooPlot that contains the data and fit.
 * @param fixedParam number of fixed parameters (i.e. number of fit parameters).
 * @þaram chi2Value value of the \f$ \chi^2 \f$.
 * @param nDof number of degrees of freedom.
 *
 * @return reduced \f$ \chi^2 \f$ value, -1. if an error occured.
 */
double ComputeChi2::operator()( RooPlot * dataFrame, const int& fixedParam,
                                double& chi2Value, int& nDof )
{
  chi2Value = -1.;
  nDof      = -1;
 
  RooHist * hist = dataFrame -> residHist(); 
  if ( hist == 0 )
    return -1.;
 
  Int_t i, nbr( hist -> GetN() ); //return num of points
 
  double result( 0. ), temp, x, y, sigma;
 
  int dof( 0 );
 
  for ( i = 0; i < nbr; i++ )
  {
    hist -> GetPoint( i, x, y );

    sigma = hist -> GetErrorY( i );
 
    //  cout<<"i  "<<i<<"  x: "<<x<<"   y:  "<<y<<"  sigma: "<<sigma<<endl;

    if ( TMath::Abs( sigma ) < 1. )
      continue;
 
    if ( x != x || y != y ) // x or y is NaN
      continue;
 
    if ( sigma != sigma )
      sigma = hist -> GetErrorYlow( i );
 
    if ( sigma != sigma )
      sigma = hist -> GetErrorYhigh( i );
 
    if ( sigma != sigma )
      continue;
 
    temp = TMath::Power( y / sigma, 2 );
 
    if ( temp != temp )
      cerr << "Error : chi2 went to NaN. Point i = " << i << ", x = "
           << x << ", y = " << y << " and sigma = " << sigma << endl;
 
    result += temp;
    dof++;

   //cout<<"NSIGMAS at x ="<<x<<": "<<y/sigma<<endl;
    //  cout<<"Chi2 = "<<temp<<endl;
  }
 
   chi2Value = result;
  nDof = dof - fixedParam;
 
  //  return result / double( nDof );

  return result;
}
 
/**
 * Same processing as ComputeChi2::operator()( RooPlot * dataFrame, const int& fixedParam,
 * double& chi2Value, int& nDof ), but only the reduced \f$ \chi^2 \f$ is accessible as output.
 * 
 * @param dataFrame pointer on the RooPlot that contains the data and fit.
 * @param fixedParam number of fixed parameters (i.e. number of fit parameters).
 */
double ComputeChi2::operator()( RooPlot * dataFrame, const int& fixedParam )
{
  double dTemp;
  int iTemp;
 
  return (*this)( dataFrame, fixedParam, dTemp, iTemp );
}
 
/**
 * Alias for ComputeChi2::operator()( RooPlot * dataFrame, const int& fixedParam,
 * double& chi2Value, int& nDof ).
 *
 * @param dataFrame pointer on the RooPlot that contains the data and fit.
 * @param fixedParam number of fixed parameters (i.e. number of fit parameters).
 * @þaram chi2Value value of the \f$ \chi^2 \f$.
 * @param nDof number of degrees of freedom.
 *
 * @return reduced \f$ \chi^2 \f$ value, -1. if an error occured.
 */
double ComputeChi2::getChi2( RooPlot * dataFrame, const int& fixedParam,
                             double& chi2Value, int& nDof )
{
  return (*this)( dataFrame, fixedParam, chi2Value, nDof );
}
 
/**
 * Alias for ComputeChi2::operator()( RooPlot * dataFrame, const int& fixedParam ).
 * 
 * @param dataFrame pointer on the RooPlot that contains the data and fit.
 * @param fixedParam number of fixed parameters (i.e. number of fit parameters).
 */
double ComputeChi2::getChi2( RooPlot * dataFrame, const int& fixedParam )
{
  double dTemp;
  int iTemp;
 
  return (*this)( dataFrame, fixedParam, dTemp, iTemp );
}


void getChi22D(RooDataSet const& data, RooAbsPdf const& pdf, RooRealVar const& x, RooRealVar const& y, double& chi2Value, int& nDof, int xBins, int yBins)
{
   TH2F* dataHist = (TH2F*)data.createHistogram("dataHist", x, Binning(xBins), YVar(y, Binning(yBins)));
   TH2F* pdfHist = (TH2F*)pdf.createHistogram("dataHist", x, Binning(xBins), YVar(y, Binning(yBins)));

   double chi2(0);
   double currentPull;

   for(int xi(0); xi<xBins; ++xi)
   {
      for(int yi(0); yi<yBins; ++yi)
      {
         currentPull = (dataHist->GetBinContent(xi+1, yi+1) - pdfHist->GetBinContent(xi+1, yi+1)) / dataHist->GetBinError(xi+1, yi+1);
         if(currentPull == currentPull) chi2 += currentPull*currentPull;
      }
   }

   int nFloatPars;
   nFloatPars = ( pdf.getParameters(data)->selectByAttrib("Constant", false) ) -> getSize();

   chi2Value = chi2;
   nDof = xBins*yBins - nFloatPars; 

}
 
//=============================================================================
 
