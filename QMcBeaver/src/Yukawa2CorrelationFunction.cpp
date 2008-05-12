#include "Yukawa2CorrelationFunction.h"
#include <sstream>

void Yukawa2CorrelationFunction::initializeParameters(
  Array1D<int> & BeginningIndexOfParameterType,
  Array1D<double> &Parameters,
  Array1D<int> & BeginningIndexOfConstantType,
  Array1D<double> & Constants)
{
  if( Constants.dim1() != 1 )
  {
    cerr << "ERROR: Yukawa correlation function entered "
	 << "without a constant for the cusp condition!" << endl;    
    exit(0);
  }

  if( BeginningIndexOfParameterType.dim1() != 1 )
  {
    cerr << "ERROR: Yukawa correlation function entered with an incorrect "
	 << "number of parameter types!" << endl;
    exit(0);
  }

  // Cusp condition
  g = Constants(0);
  // Where Williamson et al put an A, I'm putting A^2
  // Parameter related to Omega / (4 pi N),
  // where N/Omega is the number of electrons per unit volume.
  A = Parameters(0);

  if(g <= 0){
    cerr << "ERROR: Yukawa correlation function can not be used to set cusp condition "
	 << " g = " << g << endl;
    exit(0);
  }

  s2g = -1.0 * sqrt(2.0 * g);
  F   = s2g / A;
  A2  = A * A;
  A2F = A2 * F;
}

void Yukawa2CorrelationFunction::evaluate( double _r )
{
  r = _r;
  
  if(isSingular())
    {
      FunctionValue   = 0.0;
      dFunctionValue  = 0.0;
      d2FunctionValue = -1e10;
      return;
    }

  ir = 1.0 / r;
  t1 = exp( r * F );
  t2 = A2F * t1 * ir;
  
  FunctionValue   = - A2 * (1.0 - t1) * ir;
  dFunctionValue  = t2 - FunctionValue * ir;
  d2FunctionValue = t2 * F - 2.0 * dFunctionValue * ir;
  
  /*
  printf(" r = %20.10f F = %20.10f dF = %20.10f d2F = %20.10f\n",
  r,
  FunctionValue,
  dFunctionValue,
  d2FunctionValue);
  */
}

double Yukawa2CorrelationFunction::get_p_a(int ai)
{
  double p_a = 0.0;
  switch(ai)
    {
      //derivative w.r.t. A
    case 0:
      p_a = 2.0 * FunctionValue / A - t1 * s2g; 
      //printf(" r = %20.10f p = %20.10f",r,p_a);
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      clog << "Error: Yukawa Jastrow function doesn't use parameter"
	   << " index = " << ai << endl;
      break;
    }
  
  return p_a;
}

double Yukawa2CorrelationFunction::get_p2_xa(int ai)
{
  double p2_xa = 0.0;
  switch(ai)
    {
      //derivative w.r.t. g
    case 0:
      p2_xa = 2.0 * dFunctionValue / A - t1 * F * s2g;
      //printf(" p2 = %20.10f",p2_xa);
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      clog << "Error: Yukawa Jastrow function doesn't use parameter"
	   << " index = " << ai << endl;
      break;
    }
  
  return p2_xa;
}

double Yukawa2CorrelationFunction::get_p3_xxa(int ai)
{
  double p3_xxa = 0.0;
  switch(ai)
    {
      //derivative w.r.t. A
    case 0:
      p3_xxa = 2.0 * d2FunctionValue / A - t1 * F * F * s2g;
      //printf(" p3 = %20.10f\n",p3_xxa);
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      clog << "Error: Yukawa Jastrow function doesn't use parameter"
	   << " index = " << ai << endl;
      break;
    }
  return p3_xxa;
}

void Yukawa2CorrelationFunction::print(ostream& strm)
{
  strm << "Yukawa Jastrow parameters:" << endl;
  strm << "   g = " << g;
  strm << "   A = " << A << endl;
  strm << "   F = " << F << endl;

  strm << "-" << A2 << " (1.0 - Exp[ " << F << " r ])/r" << endl;
}

bool Yukawa2CorrelationFunction::isSingular()
{
  return g <= 0.0 || A <= 0.0;
}
