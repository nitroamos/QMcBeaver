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

#include "CambridgeThreeBodyCorrelationFunction.h"

void CambridgeThreeBodyCorrelationFunction::initializeParameters(
  int electron_nucleus, int electron_electron, Array1D<double> &Parameters, 
  int power, double max_dist)
{
  cutoff = max_dist;
  C = power;

  Nen = electron_nucleus;
  Nee = electron_electron;

  coeffs.allocate(Nen,Nen,Nee);

  for (int l=0; l<Nen; l++)
    for (int m=0; m<Nen; m++)
      for (int n=0; n<Nee; n++)
	coeffs(l,m,n) = Parameters(l*Nen*Nee+m*Nee+n);
}

// The code for evaluating the function and its derivatives is in this 
// function.
// We also have to figure out how to make this work with Amos's parameter
// derivatives.

void CambridgeThreeBodyCorrelationFunction::evaluate(Array1D<double> &xyz1,
						     double dist1, 
						     Array1D<double> &xyz2, 
						     double dist2)
{
  grad1.allocate(3);
  grad1 = 0.0;

  grad2.allocate(3);
  grad2 = 0.0;

  FunctionValue = 0.0;
  LaplacianValue = 0.0;

  if (dist1 > cutoff || dist2 > cutoff)
    {
      FunctionValue = 0.0;
      grad1 = 0.0;
      grad2 = 0.0;
      LaplacianValue = 0.0;
    }

  else
    {
      double d2pre1 = pow(dist1-cutoff,C-2);
      double d1pre1 = d2pre1*(dist1-cutoff);
      double pre1   = d1pre1*(dist1-cutoff);

      d1pre1 *= C;
      d2pre1 *= C*(C-1);

      double d2pre2 = pow(dist2-cutoff,C-2);
      double d1pre2 = d2pre2*(dist2-cutoff);
      double pre2   = d1pre2*(dist2-cutoff);

      d1pre2 *= C;
      d2pre2 *= C*(C-1);

      Array1D<double> xyz12(3);
      xyz12 = 0.0;

      double r12 = 0.0;
      for (int i=0; i<3; i++)
	{
	  xyz12(i) = xyz1(i)-xyz2(i);
	  r12 += xyz12(i)*xyz12(i);
	}

      r12 = sqrt(r12);
      
      Array1D<double> d1pow(Nen);
      d1pow = 1.0;
      Array1D<double> d2pow(Nen);
      d2pow = 1.0;
      Array1D<double> r12pow(Nee);
      r12pow = 1.0;

      for (int i=1; i<Nen; i++)
	{
	  d1pow(i) = d1pow(i-1)*dist1;
	  d2pow(i) = d2pow(i-1)*dist2;
	}

      for (int i=1; i<Nee; i++)
	r12pow(i) = r12pow(i-1)*r12;

      double polynomial_sum   = 0.0;

      double lpolynomial_sum  = 0.0;
      double l2polynomial_sum = 0.0;
      double lnpolynomial_sum = 0.0;

      double mpolynomial_sum  = 0.0;
      double m2polynomial_sum = 0.0;
      double mnpolynomial_sum = 0.0;

      double npolynomial_sum  = 0.0;
      double n2polynomial_sum = 0.0;

      double coeff = 0.0;

      for (int l=0; l<Nen; l++)
	for (int m=0; m<Nen; m++)
	  for (int n=0; n<Nee; n++)
	    {
	      coeff = coeffs(l,m,n);

	      polynomial_sum += 
		coeff*d1pow(l)*d2pow(m)*r12pow(n);

	      if (l>0)
		{
		  lpolynomial_sum += 
		    coeff*l*d1pow(l-1)*d2pow(m)*r12pow(n);
		  if (n>0)
		    lnpolynomial_sum +=
		      coeff*l*d1pow(l-1)*d2pow(m)*n*r12pow(n-1);
		}		    
	      if (l>1)
		l2polynomial_sum += 
		  coeff*l*(l-1)*d1pow(l-2)*d2pow(m)*r12pow(n);

	      if (m>0)
		{
		  mpolynomial_sum +=
		    coeff*d1pow(l)*m*d2pow(m-1)*r12pow(n);
		  if (n>0)
		    mnpolynomial_sum += 
		      coeff*d1pow(l)*m*d2pow(m-1)*n*r12pow(n-1);
		}
	      if (m>1)
		m2polynomial_sum += 
		  coeff*d1pow(l)*m*(m-1)*d2pow(m-2)*r12pow(n);

	      if (n>0)
		npolynomial_sum += 
		  coeff*d1pow(l)*d2pow(m)*n*r12pow(n-1);
	      if (n>1)
		n2polynomial_sum += 
		  coeff*d1pow(l)*d2pow(m)*n*(n-1)*r12pow(n-2);
	    }

      FunctionValue = pre1*pre2*polynomial_sum;

      double grad1a = d1pre1*pre2*polynomial_sum/dist1;
      double grad1b = lpolynomial_sum/dist1;
      double grad1c = npolynomial_sum/r12;

      double grad2a = pre1*d1pre2*polynomial_sum/dist2;
      double grad2b = mpolynomial_sum/dist2;
      double grad2c = grad1c;

      for (int i=0; i<3; i++)
	{
	  grad1(i) = grad1a*xyz1(i)+pre1*pre2*(grad1b*xyz1(i)+grad1c*xyz12(i));
	  grad2(i) = grad2a*xyz2(i)+pre1*pre2*(grad2b*xyz2(i)-grad2c*xyz12(i));
	}

      double dr1f = d1pre1*pre2*polynomial_sum + 
	pre1*pre2*(lpolynomial_sum + npolynomial_sum);
      
      double d2r1f = d2pre1*pre2*polynomial_sum +
	2*d1pre1*pre2*(lpolynomial_sum + npolynomial_sum) + 
	pre1*pre2*(l2polynomial_sum + 2*lnpolynomial_sum + n2polynomial_sum);

      double dr2f = pre1*d1pre2*polynomial_sum +
	pre1*pre2*(mpolynomial_sum - npolynomial_sum);

      double d2r2f = pre1*d2pre2*polynomial_sum +
	2*pre1*d1pre2*(mpolynomial_sum - npolynomial_sum) +
	pre1*pre2*(m2polynomial_sum - 2*mnpolynomial_sum + n2polynomial_sum);

      LaplacianValue = 2*dr1f/dist1 + d2r1f + 2*dr2f/dist2 + d2r2f;
    }
}

double CambridgeThreeBodyCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double CambridgeThreeBodyCorrelationFunction::get_p_a(int ai)
{
  return 0.0;
}

Array1D<double>* CambridgeThreeBodyCorrelationFunction::getElectron1Gradient()
{
  return &grad1;
}

Array1D<double>* CambridgeThreeBodyCorrelationFunction::getElectron2Gradient()
{
  return &grad2;
}

double CambridgeThreeBodyCorrelationFunction::get_p2_xa(int ai)
{
  return 0.0;
}

double CambridgeThreeBodyCorrelationFunction::getLaplacianValue()
{
  return LaplacianValue;
}

double CambridgeThreeBodyCorrelationFunction::get_p3_xxa(int ai)
{
  return 0.0;
}

double CambridgeThreeBodyCorrelationFunction::getCutoffDist()
{
  return cutoff;
}

