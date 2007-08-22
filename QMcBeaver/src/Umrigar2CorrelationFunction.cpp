#include "Umrigar2CorrelationFunction.h"
#include <sstream>

void Umrigar2CorrelationFunction::initializeParameters(
  Array1D<int> & BeginningIndexOfParameterType,
  Array1D<double> &Parameters,
  Array1D<int> & BeginningIndexOfConstantType,
  Array1D<double> & Constants)
{
  if( Constants.dim1() != 1 )
  {
    cerr << "ERROR: Umrigar correlation function entered "
      << "without a constant for the cusp condition!" << endl;    
    exit(0);
  }

  if( BeginningIndexOfParameterType.dim1() != 2 )
  {
    cerr << "ERROR: Umrigar correlation function entered with an incorrect "
      << "number of parameter types!" << endl;
    exit(0);
  }

  //Cusp condition
  g = Parameters(0);
  //parameter in denominator
  b = Parameters(1);
  //length scaling
  k = Parameters(2);

  /**
     This form might be easier to optimize since
     it permits b to change sign.
    
     Additionally, the a parameter can help scale
     the term, which might also stabilize optimization.
  */
  a = 1.0;
  B  = a * b * b;
  dB = 2.0 * a * b;

  /*
    In some circumstances, we may want to optimize the cusp condition.
    The 96 Umrigar paper describes the cusp condition as partially fulfilled
    by the HF portion of the code, so the Jastrow portion can be optimizable.
  */
  optimizeG = false;
  if(Constants(0) == 1)
    optimizeG = true;
}

void Umrigar2CorrelationFunction::evaluate( double _r )
{
  r = _r;

  dU_dr  = exp(-k * r);
  dU_drr = -1.0 * dU_dr * k;
  U = (1.0 - dU_dr) / k;

  dU_dkr  = -1.0 * r * dU_dr;
  dU_dk   = -1.0 * (dU_dkr + U) / k;
  dU_dkrr = dU_dr * ( k * r - 1.0 );

  den   = (1.0 + B * U);
  iden  = 1.0 / den;
  iden2 = iden * iden;

  FunctionValue   = g * U * iden;
  dFunctionValue  = g * dU_dr * iden2;
  d2FunctionValue = dFunctionValue * ( dU_drr / dU_dr - 2.0 * B * dU_dr * iden);
}

double Umrigar2CorrelationFunction::get_p_a(int ai)
{
  double p_a = 0.0;
  switch(ai)
    {
      //derivative w.r.t. g
    case 0:
      if(optimizeG)
	p_a = U * iden;
      else
	p_a = 0.0;
      break;

      //derivative w.r.t. b
    case 1:
      p_a  = -1.0 * dB * U * U * iden2;
      p_a *= g;
      break;

      //derivative w.r.t. k
    case 2:
      p_a  = dU_dk * iden2;
      p_a *= g;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      clog << "Error: Umrigar Jastrow function doesn't use parameter"
	   << " index = " << ai << endl;
      break;
    }

  return p_a;
}

double Umrigar2CorrelationFunction::get_p2_xa(int ai)
{
  double p2_xa = 0.0;
  switch(ai)
    {
      //derivative w.r.t. g
    case 0:
      if(optimizeG)
	p2_xa = dU_dr * iden2;
      else
	p2_xa = 0.0;
      break;

      //derivative w.r.t. b
    case 1:
      p2_xa  = -2.0 * dB * U * dU_dr * iden2 * iden;
      p2_xa *= g;
      break;

      //derivative w.r.t. k
    case 2:
      p2_xa  = -2.0 * B * dU_dk * dU_dr;
      p2_xa += den * dU_dkr;
      p2_xa *= iden2 * iden;
      p2_xa *= g;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      clog << "Error: Umrigar Jastrow function doesn't use parameter"
	   << " index = " << ai << endl;
      break;
    }

  return p2_xa;
}

double Umrigar2CorrelationFunction::get_p3_xxa(int ai)
{
  double p3_xxa = 0.0;
  switch(ai)
    {
      //derivative w.r.t. g
    case 0:
      if(optimizeG)
	{
	  p3_xxa  = -2.0 * B * dU_dr * dU_dr;
	  p3_xxa += den * dU_drr;
	  p3_xxa *= iden2 * iden;
	}
      else
	p3_xxa = 0.0;
      break;

      //derivative w.r.t. b
    case 1:
      p3_xxa  = (1.0 - 2.0 * B * U) * dU_dr * dU_dr;
      p3_xxa += U * den * dU_drr;
      p3_xxa *= -2.0 * iden2 * iden2;
      p3_xxa *= g * dB;
      break;

      //derivative w.r.t. k
    case 2:
      /*
      //this derivative was found with mathematica, where I attempted to
      //minimize the number of operations necessary. 
      p3_xxa  = r * k + r * b * (1.0 + 5.0 * dU_dr) + 2.0;
      p3_xxa += 6.0 * b * b * dU_dk * dU_dr * iden;
      p3_xxa *= iden;
      p3_xxa = dFunctionValue * (p3_xxa - 3.0);
      */
      p3_xxa  = 3.0 * B * dU_dr * dU_dr - den * dU_drr;
      p3_xxa *= 2.0 * B * dU_dk;
      p3_xxa += den * ( den * dU_dkrr - 4.0 * B * dU_dr * dU_dkr);
      p3_xxa *= iden2 * iden2;
      p3_xxa *= g;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      clog << "Error: Umrigar Jastrow function doesn't use parameter"
	   << " index = " << ai << endl;
      break;
    }
  return p3_xxa;
}

void Umrigar2CorrelationFunction::print(ostream& strm)
{
  strm << "Umrigar Jastrow parameters:" << endl;
  strm << "   g = " << g;
  if(optimizeG)
    strm << "  (optimized)" << endl;
  else
    strm << "  (not optimized)" << endl;
  strm << "   b = " << b << endl;
  strm << "   B = " << B << endl;
  strm << "   k = " << k << endl;

  stringstream Uterm;
  Uterm << " (1.0 - Exp[ -" << k << " r]) / " << k;
  strm << "  (" << g << Uterm.str() << ")/(1.0 + " << B << Uterm.str() << ")" << endl;  
}

bool Umrigar2CorrelationFunction::isSingular()
{
  return k < 0.0;
}
