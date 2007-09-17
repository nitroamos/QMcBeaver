#include "Cambridge2CorrelationFunction.h"
#include <iomanip>
#include "QMCInput.h"

void Cambridge2CorrelationFunction::initializeParameters(
  Array1D<int> & BeginningIndexOfParameterType,
  Array1D<double> &Parameters,
  Array1D<int> & BeginningIndexOfConstantType,
  Array1D<double> & Constants)
{
  active = true;
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
      /*
	Anything less than 3 will produce a discontinuity in the local energy.
      */
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
    }
  fL = f * L;

  if( L <= 0.0 || fL > 100.0)
    {
      //hopefully this isn't a guiding function...
      cerr << "Warning: bad value for L in Cambridge2CorrelationFunction." << endl
	   << "    you set L = " << L << endl
	   << "           fL = " << fL << endl
	   << "         1/fL = " << (1.0/fL) << endl;
      cerr << " This Jastrow will be inactivated." << endl;
      active = false;
    }

  int N = Parameters.dim1() - 1;  
  Array1D<double> a(N+1);

  alpha_0 = Parameters(1);
  //We are using a modified version of their Jastrow.
  //alpha_1 is the coefficient of rij in the polynomial
  alpha_1 =  g / pow( -1.0 , C ) / ( fL ) + alpha_0 * C;
  //derivative of alpha_1 w.r.t. L
  d_a1_dL = -g / pow( -1.0 , C ) / ( fL * L );

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
  if(!active) return;

  d   = fL * r - 1.0;

  //This is the Heaviside function
  if( d > 0.0 )
    return;

  //The polynomial and its derivatives
  alpha.evaluate(fL * r);
   P     =           alpha.getFunctionValue();
  dP_r   = fL *      alpha.getFirstDerivativeValue();
  dP_rr  = fL * fL * alpha.getSecondDerivativeValue();
  // "_d_ to the _p_ower of _c_" and its derivatives
  dpc    = pow(d, C);
  dpc_r  = C * fL * dpc / d;
  dpc_rr = (C - 1.0) * fL * dpc_r / d;

  if(globalInput.flags.calculate_Derivatives == 1)
    {
      dP_L   = f * r *   alpha.getFirstDerivativeValue();
      dP_Lr  = dP_L / r + f * fL * r * alpha.getSecondDerivativeValue();
      dP_Lrr = 2.0 * f *      fL * alpha.getSecondDerivativeValue() +
	         r * f * fL * fL * alpha.getThirdDerivativeValue();
      
      dpc_L  = dpc_r * r / L;
      dpc_Lr = dpc_L / d * (C * fL - 1.0 / r);  //is dividing by r here going to be a problem?
      dpc_Lrr = dpc_rr / d * (C * f * r - 2.0 / L);
    }

  //The Jastrow, and its derivatives w.r.t. r
  FunctionValue = dpc * P;

  //dFunctionValue -> g as r -> 0
  dFunctionValue = dpc_r * P + dpc * dP_r;

  d2FunctionValue =
    dpc_rr *  P         +
    dpc_r  * dP_r * 2.0 +
    dpc    * dP_rr;
}

double Cambridge2CorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double Cambridge2CorrelationFunction::get_p_a(int ai)
{
  double p_a = 0.0;
  if( d > 0.0 || !active) return 0.0;

  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(optimizeL)
	{
	  p_a  = dpc * d_a1_dL * fL * r;
	  p_a += dpc   * dP_L;
	  p_a += dpc_L *  P;
	} else {
	  p_a = 0.0;
	}
      break;
      
      //derivative w.r.t. alpha_0
    case 1:
      p_a = dpc * (1.0 + C * fL * r);
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
  if( d > 0.0  || !active) return 0.0;

  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(optimizeL)
	{
	  p2_xa  = fL * d_a1_dL * ( dpc + dpc_r * r);
	  p2_xa += dpc     * dP_Lr;
	  p2_xa += dpc_L   * dP_r;
	  p2_xa += dpc_r   * dP_L;
	  p2_xa += dpc_Lr  *  P;
	} else {
	  p2_xa = 0.0;
	}
      break;
      
      //derivative w.r.t. alpha_0
    case 1:
      p2_xa  = dpc_r * (1.0 + C * fL * r);
      p2_xa += dpc   *        C * fL;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p2_xa  = dpc_r * alpha.get_p_a(ai);
      p2_xa += dpc   * alpha.get_p2_xa(ai) * fL;
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
  if( d > 0.0  || !active ) return 0.0;

  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(optimizeL)
	{
	  p3_xxa  = dpc_Lrr * P;
	  p3_xxa += dpc_rr  * dP_L;
	  p3_xxa += dpc_Lr  * dP_r  * 2.0;
	  p3_xxa += dpc_r   * dP_Lr * 2.0;
	  p3_xxa += dpc_L   * dP_rr;
	  p3_xxa += dpc     * dP_Lrr;
	  p3_xxa += (2.0 * dpc_r + r * dpc_rr) * d_a1_dL * fL;
	} else {
	  p3_xxa = 0.0;
	}
      break;
      
      //derivative w.r.t. alpha_0
    case 1:
      p3_xxa  = dpc_rr * (1.0 + C * fL * r);
      p3_xxa += dpc_r  *  2.0 * C * fL;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p3_xxa  = dpc_rr * alpha.get_p_a(ai);
      p3_xxa += dpc_r  * alpha.get_p2_xa(ai) * 2.0 * fL;
      p3_xxa += dpc    * alpha.get_p3_xxa(ai) * fL * fL;
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
  strm << "x = r / " << (1.0 / rhs.fL) << endl;
  strm << "(x - 1)^" << rhs.C << " (";
  rhs.alpha.print(strm);
  strm << ")" << endl;
  return strm;
}

