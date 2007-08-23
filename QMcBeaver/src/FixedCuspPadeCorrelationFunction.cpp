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
  d2FunctionValue = (app - bpp*aDIVb)/b - 2*bpDIVb*dFunctionValue;
}

double FixedCuspPadeCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double FixedCuspPadeCorrelationFunction::get_p_a(int ai)
{
  if(ai < Numerator.getNumberCoefficients() - 2)
    {
      //The parameter is in the numerator. The first
      //two terms in the numerator are not optimizable, so the first
      //one is at 2.
      ai += 2;
      return Numerator.get_p_a(ai) / Denominator.getFunctionValue();
      return 0.0;
    }

  //The parameter is in the denominator, where first optimizable
  //parameter has index 1.
  ai -= (Numerator.getNumberCoefficients() - 2) -1;

  double t = -FunctionValue/Denominator.getFunctionValue();
  return t * Denominator.get_p_a(ai);
}

double FixedCuspPadeCorrelationFunction::getFirstDerivativeValue()
{
  return dFunctionValue;
}

double FixedCuspPadeCorrelationFunction::get_p2_xa(int ai)
{
  double div    = Denominator.getFunctionValue();
  double div2   = div*div;
  double  p_x   = Denominator.getFirstDerivativeValue();

  if(ai < Numerator.getNumberCoefficients() - 2)
    {
      ai += 2;
      double t = div * Numerator.get_p2_xa(ai);
      t -= Numerator.get_p_a(ai) * p_x;
      return t / div2;
    }

  ai -= (Numerator.getNumberCoefficients() - 2) -1; 

  double  p_a   = Denominator.get_p_a(ai);
  double p2_xa  = Denominator.get_p2_xa(ai);

  double t1 =  2.0 * p_a * p_x / div;
  double t2 = -1.0 * p2_xa;
  double sum = (t1 + t2) * Numerator.getFunctionValue();

  double t3 = -1.0 * p_a;
  sum += t3 * Numerator.getFirstDerivativeValue();

  return sum/div2;
}

double FixedCuspPadeCorrelationFunction::getSecondDerivativeValue()
{
  return d2FunctionValue;
}

double FixedCuspPadeCorrelationFunction::get_p3_xxa(int ai)
{
  double div    = Denominator.getFunctionValue();
  double div2   = div*div;
  double  p_x   = Denominator.getFirstDerivativeValue();
  double p2_xx  = Denominator.getSecondDerivativeValue();

  if(ai < Numerator.getNumberCoefficients() - 2)
    {
      ai += 2;
      double t = 2.0 * p_x * p_x * Numerator.get_p_a(ai);
      t -= div * p2_xx * Numerator.get_p_a(ai);
      t -= 2.0 * div * p_x * Numerator.get_p2_xa(ai);
      t += div2 * Numerator.get_p3_xxa(ai);
      return t/(div2 * div);
    }

  ai -= (Numerator.getNumberCoefficients() - 2) -1; 
  double  p_a   = Denominator.get_p_a(ai);
  double p2_xa  = Denominator.get_p2_xa(ai);
  double p3_xxa = Denominator.get_p3_xxa(ai);
  /*
    There are 7 terms. First, we'll handle the 4 terms where the numerator
    doesn't have any derivatives.
  */

  double t1 = -6.0 * p_a * p_x * p_x / div2;
  double t2 =  4.0 * p_x * p2_xa     / div;
  double t3 =  2.0 * p_a * p2_xx     / div;
  double t4 = -1.0 * p3_xxa;
  double sum = (t1 + t2 + t3 + t4) * Numerator.getFunctionValue();
  
  //Now the 2 terms with the first derivative of the numerator
  double t5 = -2.0 * p2_xa;
  double t6 =  4.0 * p_a * p_x       / div;
  sum += (t5 + t6) * Numerator.getFirstDerivativeValue();

  //Lastly, the term with the secton derivative of the numerator
  double t7 = -1.0 * p_a;
  sum += t7 * Numerator.getSecondDerivativeValue();

  return sum/div2;
}

Array1D<double> FixedCuspPadeCorrelationFunction::getNumeratorCoeffs()
{
  return Numerator.getCoefficients();
}

Array1D<double> FixedCuspPadeCorrelationFunction::getDenominatorCoeffs()
{
  return Denominator.getCoefficients();
}

void FixedCuspPadeCorrelationFunction::print(ostream& strm)
{
  strm << "Fixed Cusp Pade Jastrow parameters:" << endl;
  strm << "(";
  Numerator.print(strm);
  strm << ")/" << endl;
 
  strm << "(";
  Denominator.print(strm);
  strm << ")" << endl;
}


