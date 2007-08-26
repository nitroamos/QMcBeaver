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
#include "QMCThreeBodyJastrow.h"
#include "Stopwatch.h"
#include "IeeeMath.h"
#include "QMCWalkerData.h"

#ifdef QMC_GPU
#include "GPUQMCJastrowElectronElectron.h"
#endif

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
  QMCJastrow();

  ~QMCJastrow();

  /**
    Initializes the class with the data controlling the calculation. 
    
    @param input input data for the calculation
  */
  void initialize(QMCInput * input);

  /**
    Evaluates the Jastrow function and its derivatives at X using the 
    QMCJastrowParameters stored in the QMCInput class.

    This function only evaluates on the CPU. It is set up to process
    the walkers that the GPU did not process.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix indexed by walker
    @param num how many walkers to process
    @param start where the GPU left off (is 0 when the GPU isn't used)
    */
  void evaluate(Array1D<QMCWalkerData *> &walkerData,
		Array1D<Array2D<double>*> &X,
		int num, int start);

  void evaluate(Array2D<double> & R);

#ifdef QMC_GPU

  /**
    This gets the GPU started in evaluating Jastrow functions.

    @param a the texture ID for alpha electrons
    @param b the texture ID for beta electrons
    @param num how many walkers are to be processed
  */
  void setUpGPU(GLuint aElectrons, GLuint bElectrons, int num);

  /**
    Calculates the ElectronNuclear Jastrow data and then merges
    that with the ElectronElectron Jastrow data.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix indexed by which walker.
    @param num how many walkers are to be processed
  */
  void gpuEvaluate(Array1D<Array2D<double>*> &X, int num);

#endif

  /**
    Evaluates the Jastrow function and its derivatives at X using a
    given set of QMCJastrowParameters.

    @param JP Jastrow parameters to use during the evaluation
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */
  void evaluate( QMCJastrowParameters & JP,
		 Array1D<QMCWalkerData *> &walkerData,
		 Array1D<Array2D<double>*> &X,
		 int num, int start);

  /**
    Gets the value of the Jastrow function for the last evaluated
    electronic configuration and parameter set.  
    \f$J=exp(\sum{u_{i,j}(r_{i,j})})\f$

    @return Jastrow function value (\f$J=exp(\sum{u_{i,j}(r_{i,j})})\f$).
  */
  double getJastrow(int which);

  /**
     The partial derivative of this function with repspect to parameter
     ai.
   */
  double get_p_a(int which, int ai);
  
  /**
    Gets the value of the natural log of the Jastrow function for the last 
    evaluated electronic configuration and parameter set.  
    \f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$

    @return natural log of the Jastrow function 
    (\f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$)
  */
  double getLnJastrow(int which);

  /**
     The first partial derivative of the ln of this function with repspect to parameter
     ai.
   */
  double get_p_a_ln(int which, int ai);
  
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
     The second partial derivative of the ln of this function with repspect to parameters
     x and ai.
   */
  Array2D<double> * get_p2_xa_ln(int which, int ai);

  /**
    Gets the laplacian of the natural log of the Jastrow function with
    respect to the cartesian electronic coordinates for the last 
    evaluated electronic configuration and parameter set.  
    \f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$

    @return gradient natural log of the Jastrow function 
    (\f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$)
  */
  double getLaplacianLnJastrow(int which);

  /**
     The third partial derivative of the ln of this function with repspect to parameters
     x, x, and ai.
   */
  double get_p3_xxa_ln(int which, int ai);
  
  void operator=(const QMCJastrow & rhs );

protected:
  Array1D<double> sum_U;
  Array1D< Array2D<double> > grad_sum_U;
  Array1D<double> laplacian_sum_U;

  Array2D<double>  p_a;
  Array2D< Array2D<double> > p2_xa;
  Array2D<double> p3_xxa;
  
private:
  QMCJastrowElectronNuclear JastrowElectronNuclear;

#ifdef QMC_GPU
  GPUQMCJastrowElectronElectron gpuJEE;
#endif

  QMCJastrowElectronElectron JastrowElectronElectron;
  QMCThreeBodyJastrow ThreeBodyJastrow;
  QMCInput* Input;
};

#endif
