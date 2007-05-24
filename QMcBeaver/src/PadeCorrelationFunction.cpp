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

#include "PadeCorrelationFunction.h"

void PadeCorrelationFunction::initializeParameters(
		  Array1D<int> & BeginningIndexOfParameterType, 
		  Array1D<double> &Parameters,			    
		  Array1D<int> & BeginningIndexOfConstantType, 
		  Array1D<double> & Constants)
{
  if( BeginningIndexOfParameterType.dim1() != 2 )
    {
    cerr << "ERROR: Pade correlation function entered with an incorrect "
	 << "number of parameter types!" << endl;
    exit(0);
  }


  Array1D<double> numeratorParameters(BeginningIndexOfParameterType(1)+1);

  numeratorParameters(0) = 0.0;

  for(int i=0; i<numeratorParameters.dim1()-1; i++)
  {
    numeratorParameters(i+1) = Parameters(i);
  }

  Numerator.initialize(numeratorParameters);

  Array1D<double> denominatorParameters(Parameters.dim1()-
					numeratorParameters.dim1()+2);

  denominatorParameters(0) = 1.0;

  for(int i=0; i<denominatorParameters.dim1()-1; i++)
  {
    denominatorParameters(i+1) = 
      Parameters(i+BeginningIndexOfParameterType(1));
  }

  Denominator.initialize(denominatorParameters);
}

bool PadeCorrelationFunction::isSingular()
{
  return Denominator.hasNonNegativeZeroes();
}

Array1D<Complex> PadeCorrelationFunction::getPoles()
{
  return Denominator.getRoots();
}

void PadeCorrelationFunction::evaluate( double r )
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

  FunctionValue = a/b;
  dFunctionValue = ap/b - a*bp/b/b;
  d2FunctionValue = app/b - 2*ap*bp/b/b + 2*a*bp*bp/b/b/b - a*bpp/b/b;
}

double PadeCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double PadeCorrelationFunction::get_p_a(int ai)
{
  cout << "Parameter derivatives not implemented yet!!\n";
  return 0.0;
}

double PadeCorrelationFunction::getFirstDerivativeValue()
{
  return dFunctionValue;
}

double PadeCorrelationFunction::get_p2_xa(int ai)
{
  cout << "Parameter derivatives not implemented yet!!\n";
  return 0.0;
}

double PadeCorrelationFunction::getSecondDerivativeValue()
{
  return d2FunctionValue;
}

double PadeCorrelationFunction::get_p3_xxa(int ai)
{
  cout << "Parameter derivatives not implemented yet!!\n";
  return 0.0;
}

Array1D<double> PadeCorrelationFunction::getNumeratorCoeffs()
{
  return Numerator.getCoefficients();
}

Array1D<double> PadeCorrelationFunction::getDenominatorCoeffs()
{
  return Denominator.getCoefficients();
}






