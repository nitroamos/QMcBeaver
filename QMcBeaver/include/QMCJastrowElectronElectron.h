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
  /**
    Initializes the class with the data controling the calculation. 
    
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
  void evaluate( QMCJastrowParameters & JP, Array2D<double> & X);

  /**
    Gets the value of the natural log of the electron-electron Jastrow 
    function for the last 
    evaluated electronic configuration and parameter set.  
    \f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$

    @return natural log of the electron-electron Jastrow function 
    (\f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$)
    */
  double getLnJastrow();
  
  /**
    Gets the gradient of the natural log of the electron-electron 
    Jastrow function with
    respect to the cartesian electronic coordinates for the last 
    evaluated electronic configuration and parameter set.  
    \f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$

    @return gradient natural log of the electron-electron Jastrow function 
    (\f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$)
    */
  Array2D<double> * getGradientLnJastrow();

  /**
    Gets the laplacian of the natural log of the electron-electron 
    Jastrow function with
    respect to the cartesian electronic coordinates for the last 
    evaluated electronic configuration and parameter set.  
    \f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$

    @return gradient natural log of the electron-electron Jastrow function 
    (\f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$)
    */
  double getLaplacianLnJastrow();
  
private:    
  void calculateDistanceAndUnitVector(Array2D<double> & X1, int x1particle, 
				      Array2D<double> &X2, int x2particle, 
				      double & r, 
				      Array1D<double> & UnitVector);

  double sum_U;
  Array2D<double> grad_sum_U;
  double laplacian_sum_U;
  QMCInput* Input;
};

#endif
