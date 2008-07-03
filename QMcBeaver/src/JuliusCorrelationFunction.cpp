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

// Added to preferentially sample large electron-electron distances.

#include "JuliusCorrelationFunction.h"

void JuliusCorrelationFunction::initializeParameters(
			          Array1D<int> & BeginningIndexOfParameterType,
				  Array1D<double> &Parameters,
				  Array1D<int> & BeginningIndexOfConstantType,
				  Array1D<double> & Constants)
{
  if( BeginningIndexOfParameterType.dim1() != 2 )
    {
      cerr << "ERROR: Julius correlation function entered with an incorrect "
	   << "number of parameter types!" << endl;
      exit(0);
    }

  // Julius correlation function.
  // First constant defines function type.
  // Other constants used as parameters.
  type = Constants(0);
  if (type == 0.0)
    {
      // define Jastrow of form
      // f(r) = a0 + a1 r + a2 r^2 etc.
      // parameters given in constants section only.
      Array1D<double> numeratorParameters(BeginningIndexOfParameterType(1) + 
					  Constants.dim1() - 1);
      for (int i = 1; i < Constants.dim1(); i++)
        numeratorParameters(i - 1) = Constants(i);
      Numerator.initialize(numeratorParameters);
    }
  else if (type == 1.0)
    {
      for (int i = 1; i < Constants.dim1(); i++)
	{
	  params[i - 1] = Constants(i);
	  if (i > 10)
	    {
	      printf("Exceeded max num of parameters in ");
	      printf("JuliusCorrelationFunction.\n"); 
	      exit(1);
	    }
	}
    }

  else
    {
      printf("Jastrow type %f unknown in JuliusCorrelationFunction, exiting.\n", Constants(0)); 
      exit(1);
    }
}

bool JuliusCorrelationFunction::isSingular()
{
  return false;
}

Array1D<Complex> JuliusCorrelationFunction::getPoles()
{
  Array1D<Complex> result;
  return result;
}

void JuliusCorrelationFunction::evaluate( double r )
{
  // p = a
  // p' = a'
  // p'' = a"
  
  double a, ap, app;
  if (type == 0.0)
    {
      // define Jastrow of form
      //  f(r) = a0 + a1 r + a2 r^2 etc.
      Numerator.evaluate(r);
      a   = Numerator.getFunctionValue();
      ap  = Numerator.getFirstDerivativeValue();
      app = Numerator.getSecondDerivativeValue();
    }
  else if (type == 1.0)
    {
      // define Jastrow of form 
      //  f(r) = (1 - (a0 - a1 r^2 - a2 r^3) Exp(-((1/a0) - 1) r))^1/2
      double a0 = params[0], a1 = params[1], a2 = params[2]; 
      double t1, t2;
      t1 = exp(-(1-a0)*r/a0);
      a = sqrt(1-(a0-r*r*(a1+a2*r))* t1);
      ap = -t1*(a0*a0+r*r*(a1+a2*r)-a0*(1+a1*r*(2+r)+a2*r*r*(3+r))) / 
	(2*a0*a);
      t2 = a0*a0+r*r*(a1+a2*r)-a0*(1+a1*r*(2+r)+a2*r*r*(3+r));
      app = (t1*t1*(-t2*t2-2*(-a0+1.0/t1+r*r*(a1+a2*r))*(a0*a0*a0-r*r*(a1+a2*r)
          + a0*(1+2*a1*r*(2+r)+2*a2*r*r*(3+r))
          - a0*a0*(2+a1*(2+r*(4+r))+a2*r*(6+r*(6+r)))))) / (4*a0*a0*a*a*a);
    }
  else
    {
      printf("Jastrow type %f unknown in JuliusCorrelationFunction, exiting.\n", type); 
      exit(1);
    }
  FunctionValue = a;
  dFunctionValue = ap;
  d2FunctionValue = app;
}

double JuliusCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double JuliusCorrelationFunction::getFunctionValue(double r)
{
  /**
     I'm not going to separate the evaluate function t
  */
  if (type == 0.0)
    {
      return Numerator.getFunctionValue(r);
    }
  else if (type == 1.0)
    {
      // define Jastrow of form 
      //  f(r) = (1 - (a0 - a1 r^2 - a2 r^3) Exp(-((1/a0) - 1) r))^1/2
      double a0 = params[0], a1 = params[1], a2 = params[2]; 
      double t1;
      t1 = exp(-(1-a0)*r/a0);
      return sqrt(1-(a0-r*r*(a1+a2*r))* t1);
    } 
  else
    {
      printf("Jastrow type %f unknown in JuliusCorrelationFunction, exiting.\n", type); 
      exit(1);
    }
}

double JuliusCorrelationFunction::get_p_a(int ai)
{
  cout << "Parameter derivatives not implemented!\n";
  return 0.0;
}

double JuliusCorrelationFunction::getFirstDerivativeValue()
{
  return dFunctionValue;
}

double JuliusCorrelationFunction::get_p2_xa(int ai)
{
  cout << "Parameter derivatives not implemented!\n";
  return 0.0;
}

double JuliusCorrelationFunction::getSecondDerivativeValue()
{
  return d2FunctionValue;
}

double JuliusCorrelationFunction::get_p3_xxa(int ai)
{
  cout << "Parameter derivatives not implemented!\n";
  return 0.0;
}

Array1D<double> JuliusCorrelationFunction::getNumeratorCoeffs()
{
  return Numerator.getCoefficients();
}

Array1D<double> JuliusCorrelationFunction::getDenominatorCoeffs()
{
  Array1D<double> temp(1);
  temp(0) = 1;
  return temp;
}




