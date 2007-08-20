#include "Cambridge2CorrelationFunction.h"
#include <iomanip>

void Cambridge2CorrelationFunction::initializeParameters(
  Array1D<int> & BeginningIndexOfParameterType,
  Array1D<double> &Parameters,
  Array1D<int> & BeginningIndexOfConstantType,
  Array1D<double> & Constants)
{
  if( Constants.dim1() < 2 )
  {
    cerr << "ERROR: Cambridge 2 particle Jastrow needs two constants, gamma and C." << endl
	 << "  You requested " << Constants.dim1() << " constant parameters." << endl;
    exit(0);
  }

  g = Constants(0);
  C = (int)(Constants(1));

  if( C != 2 && C != 3 )
    {
      cerr << "ERROR: Cambridge 2 particle Jastrow requires C = 2 or C = 3, "
	   << "but you set C = " << C << endl;

    }

  if( BeginningIndexOfParameterType.dim1() != 2 || Parameters.dim1() < 2)
  {
    cerr << "ERROR: Cambridge 2 particle Jastrow needs two parameter types, L and alpha," << endl
	 << "  including at least one alpha." << endl
	 << "  You requested " << Parameters.dim1() << " parameters"
	 << " and " << BeginningIndexOfParameterType.dim1() << " parameter types." << endl;
    exit(0);
  }
  
  L = Parameters(0);
  optimizeL = false;

  if(Constants.dim1() >= 3)
    optimizeL = (int)(Constants(2)) == 1 ? true : false;
  if(Constants.dim1() >= 4)
    f = Constants(3);
  else
    {
      f = 1.0;
      //f = 1e-4;
    }

  //L = L / f;

  if( L <= 0.0 )
    {
      cerr << "Warning: Cambridge 2 particle Jastrow's L needs to be greater than 0, "
	   << "but you set L = " << L << endl;
      //If L was requested to be negative, then fabs will permit the calculation to continue
      //even if L really is bad.
      L = fabs(L);
      exit(0);
    }

  int N = Parameters.dim1() - 1;  
  Array1D<double> a(N+1);

  alpha_0 = Parameters(1);
  //We are using a modified version of their Jastrow.
  //alpha_1 is the coefficient of rij in the polynomial
  alpha_1 = g / pow( -1.0 , C ) + alpha_0 * C * f * L;
  //derivative of alpha_1 w.r.t. L
  d_a1_dL = alpha_0 * C * f;

  a(0) = alpha_0;
  a(1) = alpha_1;
  for(int ai=2; ai<=N; ai++)
    a(ai) = Parameters(ai);
  alpha.initialize(a);
}

bool Cambridge2CorrelationFunction::isSingular()
{
  return L <= 0.0;
}

Array1D<Complex> Cambridge2CorrelationFunction::getPoles()
{
  Array1D<Complex> temp;
  return temp;
}

void Cambridge2CorrelationFunction::evaluate( double _r )
{
  r = _r;
    FunctionValue = 0.0;
   dFunctionValue = 0.0;
  d2FunctionValue = 0.0;
  d   = r * f * L - 1.0;
  //This is the Heaviside function
  if( d > 0.0 )
    {
      //cout << "r = " << setw(15) << r << endl;
      return;
    }

  alpha.evaluate(r);
  double a   = alpha.getFunctionValue();
  double da  = alpha.getFirstDerivativeValue();
  double d2a = alpha.getSecondDerivativeValue();
  dpc = pow(d, C);

  //Derivatives of dpc with respect to L and r
  dpc_r  = C * f * L * dpc / d;
  dpc_rr = (C - 1.0) * f * L * dpc_r / d;
  dpc_L  = dpc_r * r / L;
  dpc_Lr = dpc_L / d * (C * f * L - 1.0 / r);  //is dividing by r here going to be a problem?
  dpc_Lrr = dpc_rr / d * (C * f * r - 2.0 / L);

  //FunctionValue = (r - L) ^ C * \sum_0^Nu a_i r^i
  FunctionValue = dpc * a;

  //dFunctionValue -> g as r -> 0
  dFunctionValue = dpc_r * a + dpc * da;

  d2FunctionValue =
    dpc_rr *   a       +
    dpc_r  *  da * 2.0 +
    dpc    * d2a;

  /*
  if(g < 0)
    {
      cout << "r = " << setw(15) << r
	   << " dpc = " << setw(15) << dpc
	   << "   a = " << setw(15) << a
	   << " F = " << setw(15) << FunctionValue
	   << " dF = " << setw(15) << dFunctionValue
	   << " d2F = " << setw(15) << d2FunctionValue << endl;
    }
  //*/
}

double Cambridge2CorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double Cambridge2CorrelationFunction::get_p_a(int ai)
{
  double p_a = 0.0;
  if( d > 0.0 ) return 0.0;

  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(optimizeL)
	{
	  p_a  = dpc * d_a1_dL * r;
	  p_a += dpc_L * alpha.getFunctionValue();
	} else {
	  p_a = 0.0;
	}
      break;
      
      //derivative w.r.t. alpha_0
    case 1:
      p_a = dpc * (1.0 + C * r * f * L);
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p_a = dpc * alpha.get_p_a(ai);
      break;
    }
  return p_a;
}

double Cambridge2CorrelationFunction::getFirstDerivativeValue()
{
  return dFunctionValue;
}

double Cambridge2CorrelationFunction::get_p2_xa(int ai)
{
  double p2_xa = 0.0;
  if( d > 0.0 ) return 0.0;

  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(optimizeL)
	{
	  p2_xa  = d_a1_dL * ( dpc + dpc_r * r);
	  p2_xa += dpc_Lr  * alpha.getFunctionValue();
	  p2_xa += dpc_L   * alpha.getFirstDerivativeValue();
	} else {
	  p2_xa = 0.0;
	}
      break;
      
      //derivative w.r.t. alpha_0
    case 1:
      p2_xa = dpc_L * (1.0 + C) * f * L * L;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p2_xa  = dpc_r * alpha.get_p_a(ai);
      p2_xa += dpc   * alpha.get_p2_xa(ai);
      break;
    }
  return p2_xa;
}

double Cambridge2CorrelationFunction::getSecondDerivativeValue()
{
  return d2FunctionValue;
}

double Cambridge2CorrelationFunction::get_p3_xxa(int ai)
{
  double p3_xxa = 0.0;
  if( d > 0.0 ) return 0.0;

  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(optimizeL)
	{
	  p3_xxa  = dpc_Lrr * alpha.getFunctionValue();
	  p3_xxa += dpc_Lr  * alpha.getFirstDerivativeValue() * 2.0;
	  p3_xxa += dpc_L   * alpha.getSecondDerivativeValue();
	  p3_xxa += (2.0 * dpc_r + r * dpc_rr) * d_a1_dL;
	} else {
	  p3_xxa = 0.0;
	}
      break;
      
      //derivative w.r.t. alpha_0
    case 1:
      p3_xxa  = dpc_Lr * (1.0 + C) * f * L * L;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p3_xxa  = dpc_rr * alpha.get_p_a(ai);
      p3_xxa += dpc_r  * alpha.get_p2_xa(ai) * 2.0;
      p3_xxa += dpc    * alpha.get_p3_xxa(ai);
      break;
    }
  return p3_xxa;
}

Array1D<double> Cambridge2CorrelationFunction::getNumeratorCoeffs()
{
  return alpha.getCoefficients();
}

Array1D<double> Cambridge2CorrelationFunction::getDenominatorCoeffs()
{
  return alpha.getCoefficients();
}

void Cambridge2CorrelationFunction::print(ostream& strm)
{
  /*
#ifdef QMC_DEBUG
  strm << *this;
#endif
  /*/
  strm << *this;
  //*/
}

ostream& operator <<(ostream& strm, Cambridge2CorrelationFunction &rhs)
{    
  strm.unsetf(ios::scientific);
  strm << "Cambridge 2 particle jastrow parameters:" << endl
       << "  g = " << rhs.g << endl
       << "  C = " << rhs.C << endl
       << "  f = " << rhs.f << endl
       << "  L = " << rhs.L;
  if(rhs.optimizeL)
    strm << " (optimized)" << endl;
  else
    strm << " (not optimized)" << endl;

  strm.setf(ios::scientific);
  strm << "  a0 = " << setw(15) << rhs.alpha_0;
  strm << "  a1 = " << setw(15) << rhs.alpha_1 << endl;
  Array1D<double> a = rhs.alpha.getCoefficients();
  for(int i=2; i<a.dim1(); i++)
    strm << "  a" << i << " = " << setw(15) << a(i);
  strm << endl;
  
  strm.unsetf(ios::scientific);
  strm << "(" << (rhs.f * rhs.L) <<  " x - 1)^" << rhs.C << " (";
  rhs.alpha.print(strm);
  strm << ")" << endl;
  return strm;
}

