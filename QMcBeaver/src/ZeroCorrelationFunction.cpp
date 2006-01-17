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


#include "ZeroCorrelationFunction.h"

void ZeroCorrelationFunction::initializeParameters(
		  Array1D<int> & BeginningIndexOfParameterType, 
		  Array1D<double> &Parameters,			    
		  Array1D<int> & BeginningIndexOfConstantType, 
		  Array1D<double> & Constants)
{
}

void ZeroCorrelationFunction::evaluate( double r )
{
}

bool ZeroCorrelationFunction::isSingular()
{
  return false;
}

Array1D<Complex> ZeroCorrelationFunction::getPoles()
{
  Array1D<Complex> result;
  return result;
}

double ZeroCorrelationFunction::getFunctionValue()
{
  return 0.0;
}

double ZeroCorrelationFunction::getFirstDerivativeValue()
{
  return 0.0;
}

double ZeroCorrelationFunction::getSecondDerivativeValue()
{
  return 0.0;
}

Array1D<double> ZeroCorrelationFunction::getNumeratorCoeffs()
{
  Array1D<double> temp(1);
  temp = 0;
  return temp;
}

Array1D<double> ZeroCorrelationFunction::getDenominatorCoeffs()
{
  Array1D<double> temp(1);
  temp = 0;
  return temp;
}

