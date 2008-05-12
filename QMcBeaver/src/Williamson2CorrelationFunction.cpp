#include "Williamson2CorrelationFunction.h"

#include "QMCInput.h"
#include <sstream>
#include "StringManipulation.h"
#include "IeeeMath.h"

//If c<2, this won't satisfy the cusp conditions
//If c<3, the 2nd derivative won't be continuous
//Hopefully the compiler will optimize some stuff
#define WJ_C_TYPE 3
#define WJ_SCALE  16
#define WJ_OPT_A

//#define SYMM_CHBY
//#define All_CHBY
#define EXP_CHBY

#define SCALE_CHBY 1

void Williamson2CorrelationFunction::initializeParameters(
  Array1D<int> & BeginningIndexOfParameterType,
  Array1D<double> &Parameters,
  Array1D<int> & BeginningIndexOfConstantType,
  Array1D<double> & Constants)
{
  if( Constants.dim1() != 1 )
  {
    cerr << "ERROR: Williamson correlation function entered "
	 << "without a constant for the cusp condition!" << endl;    
    exit(0);
  }

  if( BeginningIndexOfParameterType.dim1() != 2 )
  {
    cerr << "ERROR: Williamson correlation function entered with an incorrect "
	 << "number of parameter types!" << endl;
    exit(0);
  }

  // Cusp condition
  g = Constants(0);

  L = Parameters(0);

  // Where Williamson et al put an A, I'm putting A^2
  // Parameter related to Omega / (4 pi N),
  // where N/Omega is the number of electrons per unit volume.
  A = Parameters(1);
  B = Parameters(2);

  if(g <= 0){
    //cerr << "ERROR: Williamson correlation function can not be used to set cusp condition "
    //<< " g = " << g << endl;
    //exit(0);
    s2g = 0.0;
    F   = 0.0;
    A2  = 0.0;
    A2F = 0.0;
  } else {
    s2g = -1.0 * sqrt(2.0 * g / L);
    F   = s2g / A;
    A2  = A * A;
    A2F = A2 * F;
  }

  n = Parameters.dim1() - 3;

  Tn.allocate(n);
  dTn.allocate(n);
  d2Tn.allocate(n);

  if(n >= 0)
    {
      Tn[0]   = 1.0;
      dTn[0]  = 0.0;
      d2Tn[0] = 0.0;
    }

  if(n >= 1)
    {
      dTn[1]  = 1.0;
      d2Tn[1] = 0.0;
    }

  co.allocate(n);

  for(int i=0; i<n; i++)
    co[i] = Parameters(i+3);
}

void Williamson2CorrelationFunction::evaluate( double _r )
{
  r = _r;
  double x = r * L;

  if(isSingular())
    {
      FunctionValue   = 0.0;
      dFunctionValue  = 0.0;
      d2FunctionValue = -1e10;
      return;
    }

  // ************ This is the Yukawa type contribution 
  ir = 1.0 / x;
  t1 = exp( x * F );
  t2 = A2F * t1 * ir;

  /*
    The cusp is only accurate to r ~ 1e-5 because t1 goes to 1 too fast
    
    There are 2 choices here. The code that runs fastest tends to send
    dFunctionValue -> inf as x -> 0. So instead, I choose this arrangement
    so that the numerator goes to zero first, and the derivative -> 0.
  */
  yuk   = - A2 * (1.0 - t1) * ir;
  dyuk  = (A2 + t1*(A2F * x - A2)) * ir * ir;
  d2yuk = (-2.0 * A2 + t1 * (2.0 * A2 + A2F*x*(F * x - 2.0))) * ir * ir * ir; 

#ifdef WJ_OPT_A
  dyuk_a  = 2.0 *   yuk / A - t1 * s2g;
  d2yuk_a = 2.0 *  dyuk / A - t1 * F * s2g;
  d3yuk_a = 2.0 * d2yuk / A - t1 * F * F * s2g;
#endif

  temper   = exp( - WJ_SCALE*x*x );
  dtemper  = -2.0*WJ_SCALE*x*temper;
  d2temper = 2.0*WJ_SCALE*temper*(2.0*WJ_SCALE*x*x - 1.0);

  FunctionValue   =   yuk*temper;
  dFunctionValue  =  dyuk*temper +      yuk*dtemper;
  d2FunctionValue = d2yuk*temper + 2.0*dyuk*dtemper + yuk*d2temper;

  double dpc0, dpc1, dpc2;
  if(x > 1.0)
    {
      dG_a   = 0.0;
      dG_xa  = 0.0;
      dG_xxa = 0.0;
      dpc0   = 0.0;
      dpc1   = 0.0;
      dpc2   = 0.0;
    } else {
      dpc0 = pow(x-1,WJ_C_TYPE);
      dpc1 = pow(x-1,WJ_C_TYPE-1);
      dpc2 = pow(x-1,WJ_C_TYPE-2);

      // ************ This is the Gaussian type contribution
      
      /*
	dG_a   = 0.0;
	dG_xa  = 0.0;
	dG_xxa = 0.0;
	/*/
      dG_a   = dpc0 * (1.0 / WJ_C_TYPE + x);
      dG_xa  = dpc1 * (1.0 + WJ_C_TYPE) * x;
      dG_xxa = dpc2 * (1.0 + WJ_C_TYPE) * (WJ_C_TYPE * x - 1.0);
      
      FunctionValue   += B * dG_a;
      dFunctionValue  += B * dG_xa;
      d2FunctionValue += B * dG_xxa;
      //*/
    }

  // ************ This is the Chebyshev type contribution
  double s, sd, s2d;

#ifdef SYMM_CHBY
  //This just selects the even numbered Chebys, range [0,1]->[-1,1]
  xbar   = 2.0*x*x - 1.0;
  dxbar  = 4.0*x;
  d2xbar = 4.0;
#endif
#ifdef All_CHBY
  //This is the regular Cheby, in range [0,1]->[-1,1]
  xbar   = 2.0 * x - 1.0;
  dxbar  = 2.0;
  d2xbar = 0.0;
#endif
#ifdef EXP_CHBY
  //This is a different transformation, range [0,inf]->[1,0]
  xbar   = exp( -SCALE_CHBY*x*x );
  dxbar  = -2.0*SCALE_CHBY*x*xbar;
  d2xbar = (-2.0*SCALE_CHBY + 4.0*x*x*SCALE_CHBY*SCALE_CHBY)*xbar;
#endif

  ChebyshevT(n,xbar,s,sd,s2d);
  s2d = d2xbar*sd + dxbar*dxbar*s2d;
  sd *= dxbar;

  //This is the term that is multiplied into the Chebys to enforce
  //contraints
#ifdef EXP_CHBY
  pre   = 1.0;
  dpre  = 0.0;
  d2pre = 0.0;
#else
  pre   = dpc0 * x * x;
  dpre  = dpc1 * x * (-2.0 + (2.0+WJ_C_TYPE)*x);
  d2pre = dpc2 *(2.0 + (1.0 + WJ_C_TYPE)*x*(-4.0 + (2.0 + WJ_C_TYPE)*x));
#endif

  FunctionValue   += pre*s;
  dFunctionValue  += pre*sd + dpre*s;
  d2FunctionValue += pre*s2d + 2.0*dpre*sd + d2pre*s;

  //Convert from x back to r
  dFunctionValue  *= L;
  d2FunctionValue *= L*L;

  if(IeeeMath::isNaN(d2FunctionValue) || IeeeMath::isNaN(FunctionValue))
    {
      cout << "Bad Jastrow " << isSingular() << endl;
      print(cout);
      printf(" x = %20.10f r = %20.10f F = %20.10f dF = %20.10e d2F = %20.10e\n",
	     x,r,
	     FunctionValue,
	     dFunctionValue,
	     d2FunctionValue);
    }
  /*
  printf(" x = %20.10f r = %20.10f F = %20.10f dF = %20.10e d2F = %20.10e\n",
	 x,r,
	 FunctionValue,
	 dFunctionValue,
	 d2FunctionValue);
  //exit(0);
  //*/
}

double Williamson2CorrelationFunction::get_p_a(int ai)
{
  double p_a = 0.0;
  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(globalInput.flags.optimize_L)
	{
	  p_a = dFunctionValue * r / L;
	} else {
	  p_a = 0.0;
	}
      break;

      //derivative w.r.t. A
    case 1:
#ifdef WJ_OPT_A
      p_a = temper*dyuk_a;
#else
      p_a = 0.0;
#endif
      break;

      //derivative w.r.t. B
    case 2:
      p_a = dG_a;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p_a = pre * Tn[ai-3];
      break;
    }
  return p_a;
}

double Williamson2CorrelationFunction::get_p2_xa(int ai)
{
  double p2_xa = 0.0;
  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(globalInput.flags.optimize_L)
	{
	  p2_xa = dFunctionValue + d2FunctionValue * r / L;
	} else {
	  p2_xa = 0.0;
	}
      break;

      //derivative w.r.t. A
    case 1:
#ifdef WJ_OPT_A
      p2_xa = d2yuk_a*temper + dyuk_a*dtemper;
#else
      p2_xa = 0.0;
#endif
      break;
      
      //derivative w.r.t. B
    case 2:
      p2_xa = dG_xa;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p2_xa = dxbar*pre*dTn[ai-3] + dpre*Tn[ai-3];
      break;
    }
  return p2_xa * L;
}

double Williamson2CorrelationFunction::get_p3_xxa(int ai)
{
  double p3_xxa = 0.0;
  switch(ai)
    {
      //derivative w.r.t. L
    case 0:
      if(globalInput.flags.optimize_L)
	{
	  p3_xxa = 2.0 * d2FunctionValue / L;// + d3FunctionValue * r / L;
	} else {
	  p3_xxa = 0.0;	  
	}
      break;

      //derivative w.r.t. A
    case 1:
#ifdef WJ_OPT_A
      p3_xxa = d3yuk_a*temper + 2.0*d2yuk_a*dtemper + dyuk_a*d2temper;
#else
      p3_xxa = 0.0;
#endif
      break;

      //derivative w.r.t. B
    case 2:
      p3_xxa = dG_xxa;
      break;
      
      //derivative w.r.t. one of the terms in the summation
    default:
      p3_xxa  = d2pre*Tn[ai-3];
      p3_xxa += (2.0*dpre*dxbar + pre*d2xbar)*dTn[ai-3];
      p3_xxa += pre*dxbar*dxbar*d2Tn[ai-3];
      break;
    }
  return p3_xxa * L * L;
}

void Williamson2CorrelationFunction::print(ostream& strm)
{
  strm << "Williamson Jastrow parameters ";
#ifdef WJ_OPT_A
  strm << " (A optimized):" << endl;
#else
  strm << " (A not optimized):" << endl;
#endif

  int prec = 10;
  int width = 18;
  strm.precision(prec);
  strm << "   g = " << g;
  strm << "   A = " << A << endl;
  strm << "   B = " << B << endl;
  strm << "   L = " << L << endl;
  strm << "   F = " << F << endl;

  strm << "x = r / " << (1.0 / L) << endl;
  strm << "U[x_] := " << StringManipulation::fancyDoubleToString(prec,0,-A2)
       << " (1.0 - Exp[ " << StringManipulation::fancyDoubleToString(prec,0,F) << " x ])/x";
  strm << " Exp[ -" << WJ_SCALE << " x^2 ] ";
  strm << StringManipulation::fancyDoubleToString(prec,0,B)
       << " (1/" << WJ_C_TYPE << " + x)(x-1)^" << WJ_C_TYPE;

#ifdef EXP_CHBY
  strm << " +(";
#else
  strm << " + x^2 (x-1)^" << WJ_C_TYPE << " *(";
#endif

  for(int i=0; i<n; i++)
    {
      if(i%3 == 0) strm << endl;
      strm << StringManipulation::fancyDoubleToString(prec,width,co[i]);
#ifdef SYMM_CHBY
      strm << " ChebyshevT[" << setw(2) << 2*i << ",x]";
#endif
#ifdef All_CHBY
      strm << " ChebyshevT[" << i << ",2x-1]";
#endif
#ifdef EXP_CHBY
      strm << " ChebyshevT[" << i << ",Exp[-" << SCALE_CHBY << " x^2]]";
#endif
    }
  strm << ")" << endl;
}

bool Williamson2CorrelationFunction::isSingular()
{
  return (g > 0.0 && A < 0.0) || L <= 0.0;
}

void Williamson2CorrelationFunction::ChebyshevT(int n, double x, double &s, double &sd, double &s2d)
{
  s2d = 0.0;
  
  if (n == 0){
    s  = co[0];
    sd = 0.0;
    return;
  }

  Tn[1]   = x;

  s  = co[0] + co[1]*x;
  sd = co[1];

  if (n == 1) return;
  
  for (int k = 2; k < n; k++) {
      Tn[k] = x * 2.0 *   Tn[k-1]                  -   Tn[k-2];
     dTn[k] = x * 2.0 *  dTn[k-1] + 2.0 *  Tn[k-1] -  dTn[k-2];
    d2Tn[k] = x * 2.0 * d2Tn[k-1] + 4.0 * dTn[k-1] - d2Tn[k-2];

    s      += co[k]*Tn[k];
    sd     += co[k]*dTn[k];
    s2d    += co[k]*d2Tn[k];
  }
}

#undef WJ_C_TYPE
