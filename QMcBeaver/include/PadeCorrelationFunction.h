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

#ifndef PadeCorrelationFunction_H
#define PadeCorrelationFunction_H

#include "QMCCorrelationFunction.h"
#include "QMCPolynomial.h"


/** 
  Correlation function which uses a Pade expansion to describe 
  particle-particle interactions.  All parameters will be adjusted during 
  an optimization.
  */

class PadeCorrelationFunction : public QMCCorrelationFunction
{
private:    
  Polynomial Numerator;
  QMCPolynomial Denominator;
  double FunctionValue;
  double dFunctionValue;
  double d2FunctionValue;
  
public:

  void initializeParameters(Array1D<int> & BeginningIndexOfParameterType, 
			    Array1D<double> &Parameters,
			    Array1D<int> & BeginningIndexOfConstantType, 
			    Array1D<double> & Constants);

  void evaluate( double r );

  bool isSingular();

  Array1D<Complex> getPoles();
  
  double getFunctionValue();
  
  double getFirstDerivativeValue();
  
  double getSecondDerivativeValue();

  Array1D<double> getNumeratorCoeffs();

  Array1D<double> getDenominatorCoeffs();
};

#endif
