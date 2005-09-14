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

#include "FixedCuspPadeCorrelationFunction.h"

void FixedCuspPadeCorrelationFunction::initializeParameters(
  Array1D<int> & BeginningIndexOfParameterType,
  Array1D<double> &Parameters,
  Array1D<int> & BeginningIndexOfConstantType,
  Array1D<double> & Constants)
{
  if( Constants.dim1() < 1 )
    {
      cerr << "ERROR: Fixed Cusp Pade correlation function entered "
      << "without a constant for the cusp condition!" << endl;
      exit(0);
    }
  else if( Constants.dim1() > 1 )
    {
      cerr << "ERROR: Fixed Cusp Pade correlation function entered "
      << "with more than one constant!" << endl;
      exit(0);
    }

  if( BeginningIndexOfParameterType.dim1() != 2 )
    {
      cerr << "ERROR: Pade correlation function entered with an incorrect "
      << "number of parameter types!" << endl;
      exit(0);
    }

  Array1D<double> numeratorParameters(BeginningIndexOfParameterType(1)+2);

  numeratorParameters(0) = 0.0;
  numeratorParameters(1) = Constants(0);
  for(int i=0; i<numeratorParameters.dim1()-2; i++)
    {
      numeratorParameters(i+2) = Parameters(i);
    }

  Numerator.initialize(numeratorParameters);

  Array1D<double> denominatorParameters(Parameters.dim1()-
                                        numeratorParameters.dim1()+3);

  denominatorParameters(0) = 1.0;

  for(int i=0; i<denominatorParameters.dim1()-1; i++)
    {
      denominatorParameters(i+1) =
        Parameters(i+BeginningIndexOfParameterType(1));
    }

  Denominator.initialize(denominatorParameters);
}

bool FixedCuspPadeCorrelationFunction::isSingular()
{
  return Denominator.hasNonNegativeZeroes();
}

Array1D<Complex> FixedCuspPadeCorrelationFunction::getPoles()
{
  return Denominator.getRoots();
}

void FixedCuspPadeCorrelationFunction::evaluate( double r )
{
  // p = a/b
  // p' = a'/b - a b'/b^2
  // p'' = a"/b - 2 a' b'/b^2 + 2a (b')^2 /b^3 -a b"/b^2

  Numerator.evaluate(r);
  double a   = Numerator.getFunctionValue();
  double ap  = Numerator.getFirstDerivativeValue();
  double app = Numerator.getSecondDerivativeValue();

  Denominator.evaluate(r);
  double b   = Denominator.getFunctionValue();
  double bp  = Denominator.getFirstDerivativeValue();
  double bpp = Denominator.getSecondDerivativeValue();

  double aDIVb = a/b;
  double bpDIVb = bp/b;
  double apDIVb = ap/b;

  FunctionValue = aDIVb;
  dFunctionValue = apDIVb - bpDIVb*aDIVb;
  d2FunctionValue = (app - bpp*aDIVb)/b + 2*bpDIVb*(bpDIVb*aDIVb - apDIVb);
}

double FixedCuspPadeCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double FixedCuspPadeCorrelationFunction::getFirstDerivativeValue()
{
  return dFunctionValue;
}

double FixedCuspPadeCorrelationFunction::getSecondDerivativeValue()
{
  return d2FunctionValue;
}


