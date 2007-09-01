#include "CambridgeThreeBodyCorrelationFunction.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "QMCInput.h"

void CambridgeThreeBodyCorrelationFunction::initializeParameters(
  int electron_nucleus, int electron_electron, Array1D<double> &Parameters, 
  int power, double max_dist)
{
  cutoff = max_dist;
  C = power;

  Nen = electron_nucleus;
  Nee = electron_electron;

  coeffs.allocate(Nen,Nen,Nee);

  if(globalInput.flags.calculate_Derivatives == 1)
    {
      int numParams = Parameters.dim1();
      p_a.allocate(numParams);
      p2_x1a.allocate(3,numParams);
      p2_x2a.allocate(3,numParams);
      p3_xxa.allocate(numParams);
    }

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
						     double dist2,
						     Array1D<double> &xyz12,
						     double dist12)
{
  r1 = dist1;
  r2 = dist2;
  r12 = dist12;

  r1v = xyz1;
  r2v = xyz2;
  r12v = xyz12;

  grad1.allocate(3);
  grad1 = 0.0;

  grad2.allocate(3);
  grad2 = 0.0;

  FunctionValue = 0.0;
  LaplacianValue = 0.0;

  if (r1 > cutoff || r2 > cutoff)
    return;

  double L = 1.0 / cutoff;

  /*
    I switched this from (r - L)^C to (L' r - 1)^C
    because it's easier to change L and not have all
    the other parameters go crazy.
  */
  double d1 = L * r1 - 1.0;
  pre1 = pow(d1,C);
  d1pre1 = pre1 * C * L / d1;
  d2pre1 = (C - 1.0) * L * d1pre1 / d1;

  double d2 = L * r2 - 1.0;
  pre2 = pow(d2,C);
  d1pre2 = pre2 * C * L / d2;
  d2pre2 = (C - 1.0) * L * d1pre2 / d2;

  //calculate powers of r1, r2, and r12
  Array1D<double> d1pow(Nen);
  d1pow = 1.0;
  Array1D<double> d2pow(Nen);
  d2pow = 1.0;
  Array1D<double> r12pow(Nee);
  r12pow = 1.0;
  
  for (int i=1; i<Nen; i++)
    {
      d1pow(i) = d1pow(i-1)*r1;
      d2pow(i) = d2pow(i-1)*r2;
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

  int ai = 0;
  double term;
  double _l, _m, _n;
  double ln, mn;
  double l2, m2, n2;

  for (int l=0; l<Nen; l++)
    for (int m=0; m<Nen; m++)
      for (int n=0; n<Nee; n++)
	{
	  coeff = coeffs(l,m,n);
	  
	  //the term in the polynomial
	  term = d1pow(l)*d2pow(m)*r12pow(n);
	  
	  //terms in the first derivative of the polynomial
	  _l = l == 0 ? 0 : l*d1pow(l-1)*d2pow(m)*r12pow(n);
	  _m = m == 0 ? 0 : d1pow(l)*m*d2pow(m-1)*r12pow(n);
	  _n = n == 0 ? 0 : d1pow(l)*d2pow(m)*n*r12pow(n-1);

	  //cross terms in the second derivative of the polynomial
	  ln = l == 0 || n == 0 ? 0 : l*d1pow(l-1)*d2pow(m)*r12pow(n-1)*n;
	  mn = m == 0 || n == 0 ? 0 : d1pow(l)*m*d2pow(m-1)*n*r12pow(n-1);

	  //diagonal terms in the second derivative of the polynomial
	  l2 = l < 2 ? 0 : l*(l-1)*d1pow(l-2)*d2pow(m)*r12pow(n);
	  m2 = m < 2 ? 0 : d1pow(l)*m*(m-1)*d2pow(m-2)*r12pow(n);	  
	  n2 = n < 2 ? 0 : d1pow(l)*d2pow(m)*n*(n-1)*r12pow(n-2);

	  if(fabs(coeff) > 1.0e-100)
	    {
	        polynomial_sum += coeff * term;
	       lpolynomial_sum += coeff * _l;
	       mpolynomial_sum += coeff * _m;
	       npolynomial_sum += coeff * _n;
	      lnpolynomial_sum += coeff * ln;
	      mnpolynomial_sum += coeff * mn;
	      l2polynomial_sum += coeff * l2;
	      m2polynomial_sum += coeff * m2;
	      n2polynomial_sum += coeff * n2;
	    }

	  if(globalInput.flags.calculate_Derivatives == 1)
	    {
	      p_a(ai) = pre1*pre2*term;
	      
	      dU_dr12 = pre1*pre2*_n;
	      dU_dr1  = (d1pre1*term + pre1*_l)*pre2;
	      dU_dr2  = (d1pre2*term + pre2*_m)*pre1;
	      
	      for (int i=0; i<3; i++)
		{
		  p2_x1a(i,ai) = dU_dr1*r1v(i) + dU_dr12*r12v(i);
		  p2_x2a(i,ai) = dU_dr2*r2v(i) - dU_dr12*r12v(i);
		}
	      
	      p3_xxa(ai) = getLapPoly(term, _l, _m, _n,
				      l2, m2, n2, ln, mn);
	  
	      ai++;
	    }
	}
  
  FunctionValue = pre1*pre2*polynomial_sum;

  dU_dr12 = pre1*pre2*npolynomial_sum;
  dU_dr1  = (d1pre1*polynomial_sum + pre1*lpolynomial_sum)*pre2;
  dU_dr2  = (d1pre2*polynomial_sum + pre2*mpolynomial_sum)*pre1;

  for (int i=0; i<3; i++)
    {
      grad1(i) = dU_dr1*r1v(i) + dU_dr12*r12v(i);
      grad2(i) = dU_dr2*r2v(i) - dU_dr12*r12v(i);
    }

  LaplacianValue = getLapPoly(polynomial_sum,
			      lpolynomial_sum, mpolynomial_sum, npolynomial_sum,
			      l2polynomial_sum, m2polynomial_sum, n2polynomial_sum,
			      lnpolynomial_sum, mnpolynomial_sum);
}

double CambridgeThreeBodyCorrelationFunction::getLapPoly(double term,
							 double lterm, double mterm, double nterm,
							 double l2term, double m2term, double n2term,
							 double lnterm, double mnterm)
{
  /*
    First, I took the derivatives of U ( = term) w.r.t r1, r2, and r12. Then, the
    laplacian is:

    \nabla^2_{r_1} U + \nabla^2_{r_2} U &=& 
     \frac{4}{r_{12}} \frac{\partial U}{\partial r_{12}} + 2 \frac{\partial^2 U}{\partial r_{12}^2}\\
    &+& \frac{2}{r_1} \frac{\partial U}{\partial r_1} +      \frac{\partial^2 U}{\partial r_1^2}
     +  \frac{2}{r_2} \frac{\partial U}{\partial r_2} +      \frac{\partial^2 U}{\partial r_2^2}\\
    &+& 2 \vec{r'}_{12} \cdot
    \left(
    \vec{r'}_1 \frac{\partial^2 U}{\partial r_1 \partial r_{12}} -
    \vec{r'}_2 \frac{\partial^2 U}{\partial r_2 \partial r_{12}}
    \right)

    where \vec{r'}_1 is the unit vector of r_1, etc
  */
  double lap = 2.0*(2.0 * dU_dr12 / r12 +
		    dU_dr1  / r1 +
		    dU_dr2  / r2);

  double d2U_d2r1  = pre2*(d2pre1 * term + 2.0 * d1pre1 * lterm + pre1 * l2term);
  double d2U_d2r2  = pre1*(d2pre2 * term + 2.0 * d1pre2 * mterm + pre2 * m2term);
  double d2U_d2r12 = pre1*pre2*n2term;
  lap += 2.0 * d2U_d2r12 + d2U_d2r1 + d2U_d2r2;

  double d2U_dr1r12 = 2.0 * pre2 * (d1pre1 * nterm + pre1 * lnterm);
  double d2U_dr2r12 = 2.0 * pre1 * (d1pre2 * nterm + pre2 * mnterm);
  for(int i=0; i<3; i++)
    lap += r12v(i)*( r1v(i)*d2U_dr1r12 - r2v(i)*d2U_dr2r12 );

  return lap;
}

void CambridgeThreeBodyCorrelationFunction::print(ostream & strm)
{
  strm.unsetf(ios::scientific);
  strm << "Cambridge 3 particle jastrow parameters:" << endl
       << "Nen = " << Nen << endl
       << "Nee = " << Nee << endl
       << "  C = " << C << endl
       << "  L = " << cutoff;
  if(false)
    strm << " (optimized)" << endl;
  else
    strm << " (not optimized)" << endl;

  strm.unsetf(ios::scientific);
  strm << "(" << (1.0 / cutoff) << " r1 - 1)^" << C;
  strm << " (" << (1.0 / cutoff) << " r2 - 1)^" << C;
  strm << " (" << endl;
  int counter = 0;
  bool symmetric = true;
  for (int n=0; n<Nee; n++)
    for (int l=0; l<Nen; l++)
      for (int m=0; m<=l; m++)
	{
	  double coeff = coeffs(l,m,n);
	  if(fabs(coeff) < 1e-50 && true)
	    continue;

	  if(counter % 4 == 0 && counter != 0)
	    strm << endl;
	  counter++;

	  if(fabs(coeff - coeffs(m,l,n)) > 1e-20)
	    symmetric = false;

	  stringstream co;
	  co.precision(7);
	  co.width(10);
	  co.setf(ios::showpos);
	  co << left << coeff;

	  stringstream temp;
	  stringstream lm;
	  if(l > 0) lm << "r1";
	  if(l > 1) lm << "^" << l;
	  if(m > 0) lm << " r2";
	  if(m > 1) lm << "^" << m;

	  if(l != m)
	    {
	      lm << " + ";
	      if(l > 0) lm << "r2";
	      if(l > 1) lm << "^" << l;
	      if(m > 0) lm << " r1";
	      if(m > 1) lm << "^" << m;
	      temp << " (" << lm.str() << ")";
	    } else {
	      if(l > 0 && m > 0)
		temp << " " << lm.str();
	    }

	  if(n > 0) temp << " r12";
	  if(n > 1) temp << "^" << n;	  

	  int extrawidth = co.str().length() - 10;
	  strm << setw(10) << co.str() << setw(31-extrawidth) << left << temp.str();
	}
  strm << ")";

  if(!symmetric)
    strm << " (not symmetric!!)";

  strm << endl;
  strm << right;
}

double CambridgeThreeBodyCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double CambridgeThreeBodyCorrelationFunction::get_p_a(int ai)
{
  if (r1 > cutoff || r2 > cutoff)
    return 0.0;

  return p_a(ai);
}

Array1D<double>* CambridgeThreeBodyCorrelationFunction::getElectron1Gradient()
{
  return &grad1;
}

Array1D<double>* CambridgeThreeBodyCorrelationFunction::getElectron2Gradient()
{
  return &grad2;
}

double CambridgeThreeBodyCorrelationFunction::get_p2_xa(bool e1, int xyz, int ai)
{
  if (r1 > cutoff || r2 > cutoff)
    return 0.0;

  if(e1)
    return p2_x1a(xyz,ai);
  return p2_x2a(xyz,ai);
}

double CambridgeThreeBodyCorrelationFunction::getLaplacianValue()
{
  return LaplacianValue;
}

double CambridgeThreeBodyCorrelationFunction::get_p3_xxa(int ai)
{
  if (r1 > cutoff || r2 > cutoff)
    return 0.0;

  return p3_xxa(ai);
}

double CambridgeThreeBodyCorrelationFunction::getCutoffDist()
{
  return cutoff;
}

