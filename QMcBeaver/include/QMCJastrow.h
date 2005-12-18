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

#ifndef QMCJastrow_H
#define QMCJastrow_H

#include <iostream>
#include <math.h>

#include "Array2D.h"
#include "QMCInput.h"
#include "QMCJastrowParameters.h"
#include "QMCJastrowElectronNuclear.h"
#include "QMCJastrowElectronElectron.h"

using namespace std;


/**
  This class calculates the value of the Jastrow function and its first two 
  derivatives.

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
  \f$i\f$ and \f$j\f$.
*/

class QMCJastrow
{
public:
  /**
    Initializes the class with the data controlling the calculation. 
    
    @param input input data for the calculation
  */
  void initialize(QMCInput * input);

  /**
    Evaluates the Jastrow function and its derivatives at X using the 
    QMCJastrowParameters stored in the QMCInput class.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
  */
  void evaluate(Array1D<Array2D<double>*> &X, int num);

  /**
    Evaluates the Jastrow function and its derivatives at X using a
    given set of QMCJastrowParameters.

    @param JP Jastrow parameters to use during the evaluation
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
  */
  void evaluate( QMCJastrowParameters & JP, Array1D<Array2D<double>*> &X, int num);

  /**
    Gets the value of the Jastrow function for the last evaluated
    electronic configuration and parameter set.  
    \f$J=exp(\sum{u_{i,j}(r_{i,j})})\f$

    @return Jastrow function value (\f$J=exp(\sum{u_{i,j}(r_{i,j})})\f$).
  */
  double getJastrow(int which);
  
  /**
    Gets the value of the natural log of the Jastrow function for the last 
    evaluated electronic configuration and parameter set.  
    \f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$

    @return natural log of the Jastrow function 
    (\f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$)
  */
  double getLnJastrow(int which);
  
  /**
    Gets the gradient of the natural log of the Jastrow function with
    respect to the cartesian electronic coordinates for the last 
    evaluated electronic configuration and parameter set.  
    \f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$

    @return gradient natural log of the Jastrow function 
    (\f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$)
  */
  Array2D<double> * getGradientLnJastrow(int which);

  /**
    Gets the laplacian of the natural log of the Jastrow function with
    respect to the cartesian electronic coordinates for the last 
    evaluated electronic configuration and parameter set.  
    \f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$

    @return gradient natural log of the Jastrow function 
    (\f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$)
  */
  double getLaplacianLnJastrow(int which);
  
private:
  Array1D<double> sum_U;
  Array1D< Array2D<double> > grad_sum_U;
  Array1D<double> laplacian_sum_U;
  QMCJastrowElectronNuclear JastrowElectronNuclear;
  QMCJastrowElectronElectron JastrowElectronElectron;
  QMCInput* Input;
};

#endif
