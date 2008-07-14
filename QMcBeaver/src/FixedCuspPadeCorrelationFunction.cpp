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
  if( BeginningIndexOfParameterType.dim1() > 2 )
    {
      cerr << "Error: for FixedCuspPade Jastrows, we only permit up to 2 types of parameters, but you entered: "
	   << BeginningIndexOfParameterType << endl;
      exit(0);
    }
  if(BeginningIndexOfConstantType.dim1() > 2)
    {
      cerr << "Error: for FixedCuspPade Jastrows, we only permit up to 2 types of constants, but you entered: "
	   << BeginningIndexOfConstantType << endl;
      exit(0);
    }

  // Number of numerator constants and optimizable parameters
  if(BeginningIndexOfConstantType.dim1() == 1)
    numNC = Constants.dim1();
  else
    numNC = BeginningIndexOfConstantType(1);
  numDC = Constants.dim1() - numNC;

  // Number of denominator constants and optimizable parameters
  if(BeginningIndexOfParameterType.dim1() == 1)
    numNP = Parameters.dim1();
  else
    numNP = BeginningIndexOfParameterType(1);
  numDP = Parameters.dim1() - numNP;

  /*
    The numerator is filled with constants, then parameters.
    We assume the constant term is 0.0
  */
  Array1D<double> numeratorParameters(numNC+numNP+1);
  numeratorParameters(0) = 0.0;
  for(int i=0; i<numNC; i++)
    numeratorParameters(i+1) = Constants(i);
  for(int i=0; i<numNP; i++)
    numeratorParameters(i+numNC+1) = Parameters(i);
  Numerator.initialize(numeratorParameters);

  /*
    We assume constant term is 1.0
  */
  Array1D<double> denominatorParameters(numDC+numDP+1);
  denominatorParameters(0) = 1.0;
  for(int i=0; i<numDC; i++)
    denominatorParameters(i+1) = Constants(i+numNC);
  for(int i=0; i<numDP; i++)
    denominatorParameters(i+numDC+1) = Parameters(i+numNP);
  Denominator.initialize(denominatorParameters);
  
  /*
    Now that we're done reading arrays,
    set these to include the constant term
   */
  numDC++;
  numNC++;
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
  
  if(fabs(FunctionValue) > 300){
    print(cout);
    printf("r = %20.10e F=%20.10e dF=%20.10e d2F=%20.10e\n",r,FunctionValue,dFunctionValue,d2FunctionValue);
  }
}

double FixedCuspPadeCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double FixedCuspPadeCorrelationFunction::getFunctionValue(double r)
{
  return Numerator.getFunctionValue(r) / Denominator.getFunctionValue(r);
}

double FixedCuspPadeCorrelationFunction::get_p_a(int ai)
{
  if(ai < numNP)
    {
      //The parameter is in the numerator.
      return Numerator.get_p_a(ai+numNC) / Denominator.getFunctionValue();
    } else {      
      //The parameter is in the denominator
      ai -= numNP;      
      double t = -FunctionValue/Denominator.getFunctionValue();
      return t * Denominator.get_p_a(ai+numDC);
    }
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

  if(ai < numNP)
    {
      double t = div * Numerator.get_p2_xa(ai+numNC);
      t -= Numerator.get_p_a(ai+numNC) * p_x;
      return t / div2;
    }

  ai -= numNP;

  double  p_a   = Denominator.get_p_a(ai+numDC);
  double p2_xa  = Denominator.get_p2_xa(ai+numDC);

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

  if(ai < numNP)
    {
      double t = 2.0 * p_x * p_x * Numerator.get_p_a(ai+numNC);
      t -= div * p2_xx * Numerator.get_p_a(ai+numNC);
      t -= 2.0 * div * p_x * Numerator.get_p2_xa(ai+numNC);
      t += div2 * Numerator.get_p3_xxa(ai+numNC);
      return t/(div2 * div);
    }

  ai -= numNP;
  double  p_a   = Denominator.get_p_a(ai+numDC);
  double p2_xa  = Denominator.get_p2_xa(ai+numDC);
  double p3_xxa = Denominator.get_p3_xxa(ai+numDC);
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


