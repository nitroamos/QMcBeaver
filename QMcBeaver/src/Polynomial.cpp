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

double Polynomial::getFunctionValue()
{
  if( evaluatedF ) return f;

  f = evaluate(x,coefficients);
  evaluatedF = true;

  return f;
}

double Polynomial::getFirstDerivativeValue()
{
  if( evaluatedDF ) return df;

  df = evaluate(x,firstDerivativeCoefficients);
  evaluatedDF = true;

  return df;
}

double Polynomial::getSecondDerivativeValue()
{
  if( evaluatedD2F ) return d2f;

  d2f = evaluate(x,secondDerivativeCoefficients);
  evaluatedD2F = true;

  return d2f;
}

int Polynomial::getNumberCoefficients()
{
  return coefficients.dim1();
}
  
double Polynomial::getCoefficient(int i)
{
  return coefficients(i);
}

Array1D<Complex> Polynomial::getRoots()
{
  Array1D<Complex> complexCoeffs(coefficients.dim1());

  for(int i=0; i<coefficients.dim1(); i++)
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
        }
    }

  return roots;
}
#undef EPS





