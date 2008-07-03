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

#ifndef QMCCorrelationFunction_H
#define QMCCorrelationFunction_H

#include "Array1D.h"
#include "Complex.h"

/** 
  Interface for a parameterized function describing the interaction of 
  two particles.
  The trial wavefunction for QMC is 
  \f$\Psi_{QMC}=\Psi_{Trial}J\f$ where \f$J=exp(\sum{u_{i,j}(r_{i,j})})\f$. 
  \f$u_{ij}(r_{ij})\f$ are the QMCCorrelationFunctions describing the 
  interactions of particles \f$i\f$ and \f$j\f$.
  */

class QMCCorrelationFunction
{
public:
  /**
    Virtual destructor.
    */
  virtual ~QMCCorrelationFunction(){};

  /**
    Initializes the correlation function with a specified set of parameters.  
    This must be called every time the parameters are changed.
   */

  virtual void initializeParameters(Array1D<int> & 
	      BeginningIndexOfParameterType, 
              Array1D<double> &Parameters,
	      Array1D<int> & BeginningIndexOfConstantType, 
	      Array1D<double> &Constants) = 0;

  /** 
    Returns \f$true\f$ if the correlation function has a singularity in the
    domain \f$r\geq0\f$, and false otherwise.
    */
  virtual bool isSingular() = 0;

  /**
     Returns all of the poles of the correlation function.
  */
  virtual Array1D<Complex> getPoles() = 0;

  /**
    Evaluates the correlation function and it's first two derivatives at 
    \f$r\f$.
  */
  virtual void evaluate( double r ) = 0;
  
  /**
     Gets the value of the correlation function for the last evaluated \f$r\f$.
  */
  virtual double getFunctionValue() = 0;

  /**
     Evaluate the function as fast as possible by skipping
     the evaluation of the derivatives.
  */
  virtual double getFunctionValue(double r) = 0;
  
  /**
     Partial derivative of function with respect to parameter ai.
  */
  virtual double get_p_a(int ai) = 0;
  
  /**
     Gets the value of the first derivative of the correlation function 
     for the last evaluated \f$r\f$.
  */
  virtual double getFirstDerivativeValue() = 0;
  
  /**
     Second Partial derivative of function with respect to parameters x and ai.
  */
  virtual double get_p2_xa(int ai) = 0;
  
  /**
     Gets the value of the second derivative of the correlation function 
     for the last evaluated \f$r\f$.
  */
  virtual double getSecondDerivativeValue() = 0;
  
  /**
     Third Partial derivative of function with respect to parameters x, x, and ai.
  */
  virtual double get_p3_xxa(int ai) = 0;
  
  /**
     Returns the coefficients for the numerator of the Jastrow's function
  */
  virtual Array1D<double> getNumeratorCoeffs() = 0;
  
  /**
     Returns the coefficients for the denominator of the Jastrow's function
  */
  virtual Array1D<double> getDenominatorCoeffs() = 0;
  
  /**
     Override this function if there's some Jastrow specific
     message you want to print.
     It will be called right after the Jastrow is initialized.
  */
  virtual void print(ostream& strm){}
};


#endif





