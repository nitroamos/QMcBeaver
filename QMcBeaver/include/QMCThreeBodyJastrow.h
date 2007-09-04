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

#ifndef QMCThreeBodyJastrow_H
#define QMCThreeBodyJastrow_H

#include <string>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCInput.h"
#include "QMCJastrowParameters.h"
#include "QMCWalkerData.h"

using namespace std;

/**
  This class calculates the value of the three body part of the Jastrow 
  function and its first two derivatives.

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
  This class contains functions describing correlations between a nucleus and
  two electrons.
 */

class QMCThreeBodyJastrow
{
public:
  QMCThreeBodyJastrow();

  ~QMCThreeBodyJastrow();

  /**
    Initializes the class with the data controlling the calculation. 
    
    @param input input data for the calculation
    */
  void initialize(QMCInput * input);

  /**
    Evaluates the three body Jastrow function and its derivatives at X using a 
    given set of QMCJastrowParameters.

    @param JP Jastrow parameters to use during the evaluation
    @param X \f$3N\f$ dimensional configuration of electrons represented by a 
    \f$N \times 3\f$ matrix
    */
  void evaluate( QMCJastrowParameters & JP, QMCWalkerData * wData,
		 Array2D<double> & X);

  /**
     Partial derivative of the natural log of this function with respect to 
     parameter ai.
  */
  double get_p_a_ln(int ai);

  /**
     Second partial derivative of the natural log of this function with respect
     to parameters x and ai.
  */
  Array2D<double> * get_p2_xa_ln(int ai);

  /**
     Third partial derivative of the natural log of this function with respect
     to parameters x, x, and ai.
  */
  double get_p3_xxa_ln(int ai);

  /**
    Gets the value of the natural log of the three body Jastrow function for 
    the last evaluated electronic configuration and parameter set.  
    \f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$

    @return natural log of the three body Jastrow function 
    (\f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$)
    */
  double getLnJastrow();
  
  /**
    Gets the gradient of the natural log of the three body Jastrow function 
    with respect to the cartesian electronic coordinates for the last evaluated
    electronic configuration and parameter set.  
    \f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$

    @return gradient natural log of the three body Jastrow function 
    (\f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$)
    */
  Array2D<double> * getGradientLnJastrow();
  
  /**
    Gets the laplacian of the natural log of the three body Jastrow function 
    with respect to the cartesian electronic coordinates for the last evaluated
    electronic configuration and parameter set.  
    \f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$

    @return laplacian natural log of the three body Jastrow function 
    (\f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$)
    */
  double getLaplacianLnJastrow();

 protected:
  double sum_U;
  Array2D<double> grad_sum_U;
  double laplacian_sum_U;
  QMCInput* Input;

  QMCWalkerData * wd;

  Array1D<double> p_a;
  Array1D< Array2D<double> > p2_xa;
  Array1D<double> p3_xxa;

  Array1D<QMCThreeBodyCorrelationFunctionParameters> * EupEdnNuclear;
  Array1D<QMCThreeBodyCorrelationFunctionParameters> * EupEupNuclear;
  Array1D<QMCThreeBodyCorrelationFunctionParameters> * EdnEdnNuclear;

private:    
  void packageDerivatives();

  void calculateDistances(Array2D<double> &X1, int x1particle, int x2particle, 
			  Array2D<double> &X3, int x3particle, 
			  Array1D<double> &position1, double &r1, 
			  Array1D<double> &position2, double &r2);

  void collectForPair(int Electron1, 
		      int Electron2,
		      int Nuclei,
		      QMCThreeBodyCorrelationFunctionParameters * paramset);
  
};
#endif