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

#ifndef QMCThreeBodyCorrelationFunction_H
#define QMCThreeBodyCorrelationFunction_H

#include "Array1D.h"

/** 
  Interface for a parameterized function describing the interaction of a
  a nucleus and two electrons.
  The trial wavefunction for QMC is 
  \f$\Psi_{QMC}=\Psi_{Trial}J\f$ where \f$J=exp(\sum{u_{i,j}(r_{i,j})})\f$. 
  \f$u_{ij}(r_{ij})\f$ are the QMCCorrelationFunctions describing the 
  interactions of particles \f$i\f$ and \f$j\f$.
*/

class QMCThreeBodyCorrelationFunction
{
 public:
  /**
     Virtual destructor.
  */
  virtual ~QMCThreeBodyCorrelationFunction(){};

  /**
     Initializes the correlation function with a specified set of parameters.  
     This must be called every time the parameters are changed.
  */
  virtual void initializeParameters(int electron_nucleus, 
				    int electron_electron, 
				    Array1D<double> &Parameters,int power, 
				    double max_dist) = 0;

  virtual bool setElectron(bool first, Array1D<double> &xyz, double dist) = 0;

  /**
    Evaluates the correlation function and its first two derivatives at 
    \f$r\f$.
    */
  virtual void evaluate(Array1D<double> &xyz12, double r12) = 0;
  
  /**
    Gets the value of the correlation function for the last evaluated \f$r\f$.
    */
  virtual double getFunctionValue() = 0;

  /**
     Evaluate the function as fast as possible by skipping
     the evaluation of the derivatives.
  */
  virtual double getFunctionValue(double r12, double r1, double r2) = 0;

  /**
     Partial derivative of function with respect to parameter ai.
  */
  virtual double get_p_a(int ai) = 0;
  
  /**
     Gets the gradient of the correlation function for electron 1 at the last
     evaluated configuration.
  */
  virtual Array1D<double> * getElectron1Gradient() = 0;

  /**
     Gets the gradient of the correlation function for electron 2 at the last 
     evaluated configuration.
  */
  virtual Array1D<double> * getElectron2Gradient() = 0;

  /**
     Second Partial derivative of function with respect to parameters x and ai.
  */
  virtual double get_p2_xa(bool e1, int xyz, int ai) = 0;
  
  /**
     Gets the value of the Laplacian of the correlation function with respect
     to electrons one and two at the last evaluated configuration.
  */
  virtual double getLaplacianValue() = 0;

  /**
     Third Partial derivative of function with respect to parameters x, x, and 
     ai.
  */
  virtual double get_p3_xxa(int ai) = 0;

  /**
     Returns the cutoff for the electron-nucleus distance for this function.
  */
  virtual double getCutoffDist() = 0;

  /**
     Override this function if there's some Jastrow specific
     message you want to print.
     It will be called right after the Jastrow is initialized.
  */
  virtual void print(ostream& strm){}
};


#endif





