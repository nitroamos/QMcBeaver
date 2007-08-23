//            QMcBeaver
//
//         Constructed by
//
//     Michael Todd Feldmann
//              and
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   
 
Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/


#include "Polynomial.h"

Polynomial::Polynomial()
{
  initialize();
}

Polynomial::Polynomial(Array1D<double> & coeffs)
{
  initialize();
  initialize(coeffs);
}

void Polynomial::operator=(Polynomial rhs)
{
  initialize();
  initialize(rhs.coefficients);
}

void Polynomial::initialize()
{
  evaluatedF   = false;
  evaluatedDF  = false;
  evaluatedD2F = false;

  // Assign some non-random values so that bugs can be traced faster
  f   = 1234.5;
  df  = 2345.6;
  d2f = 3456.7;
  x   = 4567.8;
}

void Polynomial::initialize(Array1D<double> & coeffs)
{
  // Polynomial's coefficients
  coefficients = coeffs;

  // Polynomial's first derivative coefficients
  firstDerivativeCoefficients.allocate(coefficients.dim1()-1);

  for(int i=0; i<firstDerivativeCoefficients.dim1(); i++)
    {
      firstDerivativeCoefficients(i) = (i+1)*coefficients(i+1);
    }

  // Polynomial's second derivative coefficients
  secondDerivativeCoefficients.allocate(coefficients.dim1()-2);

  for(int i=0; i<secondDerivativeCoefficients.dim1(); i++)
    {
      secondDerivativeCoefficients(i) = (i+2)*(i+1)*coefficients(i+2);
    }
}

void Polynomial::evaluate(double x)
{
  this->x = x;
  evaluatedF   = false;
  evaluatedDF  = false;
  evaluatedD2F = false;
}

double Polynomial::evaluate(double x, Array1D<double> &coeffs)
{
  this->x = x;
  int n = coeffs.dim1()-1;

  if( n < 0 )
    {
      return 0;
    }

  double val = coeffs(n);
  for(int i=n-1; i >= 0; i--)
    {
      val = val*x + coeffs(i);
    }

  return val;
}

void Polynomial::evaluateAll(double x, Array1D<double> & coeffs)
{
  this->x = x;
  //This method was sort of adapted from 5.3 of Numerical Recipies.
  evaluatedF = true;
  evaluatedDF = true;
  evaluatedD2F = true;

  int n = coeffs.dim1()-1;

  if( n < 0 ) return;

  f = coeffs(n);
  df = 0.0;
  d2f = 0.0;
  for(int i=n-1; i >= 0; i--)
    {
      d2f = d2f*x + df;
      df  = df*x  + f;
      f   = f*x   + coeffs(i);
    }
  d2f *= 2;
}

double Polynomial::getFunctionValue()
{
  if( evaluatedF ) return f;

  evaluateAll(x,coefficients);
  //f = evaluate(x,coefficients);
  //evaluatedF = true;

  return f;
}

double Polynomial::get_p_a(int ai)
{
  return pow(x,ai);
}

double Polynomial::getFirstDerivativeValue()
{
  if( evaluatedDF ) return df;

  evaluateAll(x,coefficients);
  //df = evaluate(x,firstDerivativeCoefficients);
  //evaluatedDF = true;

  return df;
}

double Polynomial::get_p2_xa(int ai)
{
  //only one term has a_i in it.
  return ai * pow(x,ai-1);
}

double Polynomial::getSecondDerivativeValue()
{
  if( evaluatedD2F ) return d2f;

  evaluateAll(x,coefficients);
  //d2f = evaluate(x,secondDerivativeCoefficients);
  //evaluatedD2F = true;

  return d2f;
}

double Polynomial::get_p3_xxa(int ai)
{
  return (ai-1) * ai * pow(x,ai-2);
}

int Polynomial::getNumberCoefficients()
{
  return coefficients.dim1();
}

double Polynomial::getCoefficient(int i)
{
  return coefficients(i);
}

Array1D<double> Polynomial::getCoefficients()
{
  return coefficients;
}

void Polynomial::print(ostream & strm)
{
  strm.unsetf(ios::fixed);
  strm.unsetf(ios::scientific);
  for(int i=0; i<coefficients.dim1(); i++)
    {
      strm << coefficients(i);
      if(i > 0)
	strm << " x^" << i;
      if(i < coefficients.dim1() - 1 &&
	 coefficients(i+1) >= 0)
	strm << " +";
      else
	strm << " ";
    }
}

Array1D<Complex> Polynomial::getRoots()
{
  int numRoots = coefficients.dim1();

  /**
     If the highest order terms have 0 for
     a coefficient, then we have fewer roots.
  */
  for(int i=coefficients.dim1()-1; i>=0; i--)
    if( fabs(coefficients(i)) <  1.0e-50)
      numRoots--;
    else
      break;

  Array1D<Complex> complexCoeffs(numRoots);

  for(int i=0; i<numRoots; i++)
    {
      complexCoeffs(i) = coefficients(i);
    }

  bool polish = true;
  bool calcOK = true;

  Array1D<Complex> roots = zroots(complexCoeffs, polish, &calcOK);

  if( !calcOK )
    {
      throw Exception("ERROR: Problems during the calculation of Polynomial roots!");
    }

  return roots;
}

/*
  EPSS is the fractional estimated roundoff error.  1.0e-7 was the default
  for a float point implementation.  1.0e-14 is probably safe for this double 
  implementation.
 
  Try to break (rare) limit cycles with MR different fractional values once
  every MT steps for MAXIT total allowed iterations.
  */

#define EPSS 1.0e-14 //1.0e-7 default for float implementation
#define MR 8         //don't change this unless you know how to change frac
#define MT 10
#define MAXIT (MT*MR)

void Polynomial::laguer(Array1D<Complex> &a, int m, Complex &x, int *its,
                        bool *calcOK)
{
  Complex cZero(0.0,0.0);

  double abx,abp,abm,err;
  Complex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
  static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

  *calcOK = *calcOK && true;

  for(int iter=1; iter<=MAXIT; iter++)
    {
      *its = iter;
      b = a(m);
      err = b.abs();
      d = cZero;
      f = cZero;
      abx = x.abs();

      for(int j=m-1;j>=0;j--)
        {
          f = (x*f)+d;
          d = (x*d)+b;
          b = (x*b)+a(j);
          err = b.abs() + abx*err;
        }

      err *= EPSS;

      if( b.abs() <= err) return;

      g = d/b;
      g2 = g*g;
      h = g2-((f/b)*2.0);
      sq = (((h*(double)m)-g2)*(double)(m-1)).squareroot();
      gp = g+sq;
      gm = g-sq;
      abp = gp.abs();
      abm = gm.abs();

      if (abp < abm) gp=gm;

      dx=( (abp>abm?abp:abm) > 0.0 ? (Complex(m,0.0)/gp)
           : (Complex(cos((double)iter),sin((double)iter)) * (1+abx)) );

      x1 = x-dx;

      if (x.real() == x1.real() && x.imaginary() == x1.imaginary()) return;

      if (iter % MT) x = x1;
      else x = x-(dx*frac[iter/MT]);
    }

  *calcOK = false;
  return;
}

#undef EPSS
#undef MR
#undef MT
#undef MAXIT



#define EPS 2.0e-13  // 2.0e-6 default in float implementation

Array1D<Complex> Polynomial::zroots(Array1D<Complex> & a, bool polish,
                                    bool *calcOK)
{
  *calcOK = true;

  int m = a.dim1() - 1;

  Array1D<Complex> roots(m);

  int its;
  Complex x,b,c;
  Array1D<Complex> ad(m+1);

  for(int j=0;j<=m;j++)
    {
      ad(j) = a(j);
    }

  for(int j=m;j>=1;j--)
    {
      x = 0.0;

      laguer(ad,j,x,&its,calcOK);

      if (fabs(x.imaginary()) <= 2.0*EPS*fabs(x.real()))
        {
          x = Complex(x.real(),0.0);
        }

      roots(j-1) = x;
      b = ad(j);

      for (int jj=j-1;jj>=0;jj--)
        {
          c = ad(jj);
          ad(jj) = b;
          b = (x*b)+c;
        }
    }
  if(polish)
    {
      for(int j=1;j<=m;j++)
        {
          laguer(a,m,roots(j-1),&its,calcOK);
        }
    }

  for (int j=2;j<=m;j++)
    {
      x = roots(j-1);

      int i;
      for (i=j-1;i>=1;i--)
        {
          if ( roots(i-1).real() <= x.real() ) break;
          roots(i) = roots(i-1);
        }
      roots(i) = x;
    }

  return roots;
}
#undef EPS





