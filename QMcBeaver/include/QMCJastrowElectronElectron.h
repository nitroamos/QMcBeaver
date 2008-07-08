
//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef QMCJastrowElectronElectron_H
#define QMCJastrowElectronElectron_H

#include "Array1D.h"
#include "Array2D.h"
#include "QMCInput.h"
#include "QMCJastrowParameters.h"
#include "QMCWalkerData.h"

/**
  This class calculates the value of the electron-electron part of the 
  Jastrow function and its first two derivatives.

  The wavefunction is assumed to be of the form
  \f[
  \Psi_{QMC} = \Psi_{Trial}J
  \f]
  where \f$\Psi_{Trial}\f$ is a wavefunction calculated using a standard QM
  method and
  \f[
  J=exp(\sum{u_{i,j}(r_{i,j})})
  \f]
  is a Jastrow type correlation function.  \f$u_{ij}(r_{ij})\f$ are 
  QMCCorrelationFunction describing the interactions of particles 
  \f$i\f$ and \f$j\f$.  The sum can be broken up into electron-electron and
  electron-nuclear components.
 */

class QMCJastrowElectronElectron
{
public:
  QMCJastrowElectronElectron();

  ~QMCJastrowElectronElectron();
  
  /**
    Initializes the class with the data controlling the calculation. 
    
    @param input input data for the calculation
    */
  void initialize(QMCInput * input);

  /**
    Evaluates the electron-electron Jastrow function and its derivatives at X 
    using a given set of QMCJastrowParameters.

    @param JP Jastrow parameters to use during the evaluation
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */
  void evaluate( QMCJastrowParameters & JP,
		 QMCWalkerData * wData,
		 Array2D<double> & X);

  /**
     Partial derivative of the natural log of the function with respect to
     parameter ai.
  */
  double get_p_a_ln(int ai);

  /**
     Second Partial derivative of the natural log of the function with respect to
     parameters x and ai.
  */
  Array2D<double> * get_p2_xa_ln(int ai);

  /**
     Third partial derivative of the natural log of the function with respect to
     parameters x, x, and ai.
  */
  double get_p3_xxa_ln(int ai);

  /**
     When calculating the overlap integral for pseudopotentials, we need to
     evaluate the ratio of the Jastrow at each grid point to the Jastrow at
     its original position.
     
     @param E is the index of the electron being moved to the different grid points
     @param R the array of all electron positions
     @param grid the grid points for electron E
     @param integrand the value of the Jastrow at each of the grid points
     @return the Jastrow at the original position (the denominator)
  */
  double jastrowOnGrid(QMCJastrowParameters & JP,
		       int E,
		       Array2D<double> & R,
		       Array2D<double> & grid,
		       Array1D<double> & integrand);

protected:
  Array1D<double>  p_a;
  Array1D< Array2D<double> > p2_xa;
  Array1D<double> p3_xxa;

private:  
  QMCWalkerData * wd;

  /**
     The QMCJastrowElectronElectron::evaluate redirects here if
     we are only evaluating 1 electron.
     
     @param JP Jastrow parameters to use during the evaluation
     @param X \f$3N\f$ dimensional configuration of electrons represented by 
     a \f$N \times 3\f$ matrix
    */
  void updateOne( QMCJastrowParameters & JP,
		  Array2D<double> & X);

  /**
     The QMCJastrowElectronElectron::evaluate redirects here if
     we are evaluating all electrons.
     
     @param JP Jastrow parameters to use during the evaluation
     @param X \f$3N\f$ dimensional configuration of electrons represented by 
     a \f$N \times 3\f$ matrix
    */
  void updateAll( QMCJastrowParameters & JP,
		  Array2D<double> & X);
  
  void collectForPair(int el1, int el2,
		      QMCCorrelationFunction *U_Function,
		      Array2D<double> & X,
		      int index, int numP);
};

#endif
