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

#ifndef FixedCuspPadeCorrelationFunction_H
#define FixedCuspPadeCorrelationFunction_H

#include "QMCCorrelationFunction.h"
#include "QMCPolynomial.h"


/** 
  Correlation function which uses a Pade expansion to describe 
  particle-particle interactions.  The cusp condition is a fixed constant, 
  and all other parameters will be adjusted during an optimization.
  */

class FixedCuspPadeCorrelationFunction : public QMCCorrelationFunction
{
 protected:
  double FunctionValue;
  double dFunctionValue;
  double d2FunctionValue;
  
 private:    
  Polynomial Numerator;
  QMCPolynomial Denominator;
  
 public:

  void initializeParameters(Array1D<int> & BeginningIndexOfParameterType, 
			    Array1D<double> &Parameters,
			    Array1D<int> & BeginningIndexOfConstantType, 
			    Array1D<double> & Constants);

  void evaluate( double r );

  bool isSingular();

  Array1D<Complex> getPoles();

  double getFunctionValue();
  double get_p_a(int ai);
  
  double getFirstDerivativeValue();
  double get_p2_xa(int ai);
  
  double getSecondDerivativeValue();
  double get_p3_xxa(int ai);

  Array1D<double> getNumeratorCoeffs();

  Array1D<double> getDenominatorCoeffs();
};

#endif
