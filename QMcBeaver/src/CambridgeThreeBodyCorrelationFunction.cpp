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
  L = cutoff;
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
      p2_x1L.allocate(3);
      p2_x2L.allocate(3);
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
  r1 = dist1 * L;
  r2 = dist2 * L;
  r12 = dist12 * L;

  r1v = xyz1;
  r2v = xyz2;
  r12v = xyz12;

  grad1.allocate(3);
  grad1 = 0.0;

  grad2.allocate(3);
  grad2 = 0.0;

  FunctionValue  = 0.0;
  LaplacianValue = 0.0;

  d1 = r1 - 1.0;
  d2 = r2 - 1.0; 
  if(d1 > 0.0 || d2 > 0.0) return;

  /*
    I switched this from (r - L)^C to (L' r - 1)^C
    because it's easier to change L and not have all
    the other parameters go crazy.
  */
  pre1 = pow(d1,C);
  d1pre1 = pre1 * C / d1;
  d2pre1 = (C - 1.0) * d1pre1 / d1;
  
  pre2 = pow(d2,C);
  d1pre2 = pre2 * C / d2;
  d2pre2 = (C - 1.0) * d1pre2 / d2;

  if(globalInput.flags.calculate_Derivatives == 1)
    {
      pre1_L = d1pre1 * r1 / L;
      pre2_L = d1pre2 * r2 / L;

      d1pre1_L = pre1_L / d1 * (C - 1.0 / r1);
      d1pre2_L = pre2_L / d2 * (C - 1.0 / r2);

      d2pre1_L = d2pre1 / d1 * (C * r1 - 2.0) / L;
      d2pre2_L = d2pre2 / d2 * (C * r2 - 2.0) / L;

      p2_x1L = 0.0;
      p2_x2L = 0.0;
      dU_L = 0.0;
      dU_Lrr = 0.0;
    }

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
	      
	      double dL = coeff * (l + m + n) / L;
	      
	      for (int i=0; i<3; i++)
		{
		  p2_x1a(i,ai) = L*(dU_dr1*r1v(i) + dU_dr12*r12v(i));
		  p2_x2a(i,ai) = L*(dU_dr2*r2v(i) - dU_dr12*r12v(i));
		  p2_x1L(i) += dL * p2_x1a(i,ai);
		  p2_x2L(i) += dL * p2_x2a(i,ai);
		}
	      
	      p3_xxa(ai) = getLapPoly(term, _l, _m, _n,
				      l2, m2, n2, ln, mn,
				      pre1, pre2, d1pre1,
				      d1pre2, d2pre1, d2pre2);

	      dU_L += dL * p_a(ai);
	      dU_Lrr += dL * p3_xxa(ai);

	      ai++;
	    }
	}

  FunctionValue = pre1*pre2*polynomial_sum;

  dU_dr12 = pre1*pre2*npolynomial_sum;
  dU_dr1  = (d1pre1*polynomial_sum + pre1*lpolynomial_sum)*pre2;
  dU_dr2  = (d1pre2*polynomial_sum + pre2*mpolynomial_sum)*pre1;

  for (int i=0; i<3; i++)
    {
      grad1(i) = L*(dU_dr1*r1v(i) + dU_dr12*r12v(i));
      grad2(i) = L*(dU_dr2*r2v(i) - dU_dr12*r12v(i));
    }

  LaplacianValue = getLapPoly(polynomial_sum,
			      lpolynomial_sum, mpolynomial_sum, npolynomial_sum,
			      l2polynomial_sum, m2polynomial_sum, n2polynomial_sum,
			      lnpolynomial_sum, mnpolynomial_sum,
			      pre1, pre2, d1pre1, d1pre2, d2pre1, d2pre2);

  if(globalInput.flags.calculate_Derivatives == 1 &&
     globalInput.flags.optimize_L == 1)
    {
      dU_L += pre1_L * pre2 * polynomial_sum;
      dU_L += pre1 * pre2_L * polynomial_sum;

      dU_dr12 = pre1_L*pre2*npolynomial_sum;
      dU_dr1  = (d1pre1_L*polynomial_sum + pre1_L*lpolynomial_sum)*pre2;
      dU_dr2  = (d1pre2*polynomial_sum + pre2*mpolynomial_sum)*pre1_L;

      for (int i=0; i<3; i++)
	{
	  p2_x1L(i) += L*(dU_dr1*r1v(i) + dU_dr12*r12v(i));
	  p2_x2L(i) += L*(dU_dr2*r2v(i) - dU_dr12*r12v(i));
	}

      dU_Lrr += getLapPoly(polynomial_sum,
			   lpolynomial_sum, mpolynomial_sum, npolynomial_sum,
			   l2polynomial_sum, m2polynomial_sum, n2polynomial_sum,
			   lnpolynomial_sum, mnpolynomial_sum,
			   pre1_L, pre2, d1pre1_L, d1pre2, d2pre1_L, d2pre2);
      
      dU_dr12 = pre1*pre2_L*npolynomial_sum;
      dU_dr1  = (d1pre1*polynomial_sum + pre1*lpolynomial_sum)*pre2_L;
      dU_dr2  = (d1pre2_L*polynomial_sum + pre2_L*mpolynomial_sum)*pre1;

      for (int i=0; i<3; i++)
	{
	  p2_x1L(i) += L*(dU_dr1*r1v(i) + dU_dr12*r12v(i));
	  p2_x2L(i) += L*(dU_dr2*r2v(i) - dU_dr12*r12v(i));
	}

      dU_Lrr += getLapPoly(polynomial_sum,
			   lpolynomial_sum, mpolynomial_sum, npolynomial_sum,
			   l2polynomial_sum, m2polynomial_sum, n2polynomial_sum,
			   lnpolynomial_sum, mnpolynomial_sum,
			   pre1, pre2_L, d1pre1, d1pre2_L, d2pre1, d2pre2_L);
    }
}

double CambridgeThreeBodyCorrelationFunction::getLapPoly(double term,
							 double lterm, double mterm, double nterm,
							 double l2term, double m2term, double n2term,
							 double lnterm, double mnterm,
							 double   p1, double   p2,
							 double d1p1, double d1p2,
							 double d2p1, double d2p2)
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

  double d2U_d2r1  = p2*(d2p1 * term + 2.0 * d1p1 * lterm + p1 * l2term);
  double d2U_d2r2  = p1*(d2p2 * term + 2.0 * d1p2 * mterm + p2 * m2term);
  double d2U_d2r12 = p1*p2*n2term;
  lap += 2.0 * d2U_d2r12 + d2U_d2r1 + d2U_d2r2;

  double d2U_dr1r12 = 2.0 * p2 * (d1p1 * nterm + p1 * lnterm);
  double d2U_dr2r12 = 2.0 * p1 * (d1p2 * nterm + p2 * mnterm);
  for(int i=0; i<3; i++)
    lap += r12v(i)*( r1v(i)*d2U_dr1r12 - r2v(i)*d2U_dr2r12 );

  return lap * L * L;
}

double CambridgeThreeBodyCorrelationFunction::getFunctionValue()
{
  return FunctionValue;
}

double CambridgeThreeBodyCorrelationFunction::get_p_a(int ai)
{
  if(d1 > 0.0 || d2 > 0.0) return 0.0;

  if(globalInput.flags.optimize_L == 1)
    {
      if(ai == 0) return dU_L;
    } else {
      if(ai == 0) return 0.0;
    }
  return p_a(ai-1);
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
  if(d1 > 0.0 || d2 > 0.0) return 0.0;

  if(globalInput.flags.optimize_L == 1)
    {
      if(ai == 0)
	{
	  if(e1) return p2_x1L(xyz);
	  else   return p2_x2L(xyz);
	}
    } else {
      if(ai == 0) return 0.0;
    }
  if(e1) return p2_x1a(xyz,ai-1);
  return p2_x2a(xyz,ai-1);
}

double CambridgeThreeBodyCorrelationFunction::getLaplacianValue()
{
  return LaplacianValue;
}

double CambridgeThreeBodyCorrelationFunction::get_p3_xxa(int ai)
{
  if(d1 > 0.0 || d2 > 0.0) return 0.0;

  if(globalInput.flags.optimize_L == 1)
    {
      if(ai == 0) return dU_Lrr;
    } else {
      if(ai == 0) return 0.0;
    }
  return p3_xxa(ai-1);
}

double CambridgeThreeBodyCorrelationFunction::getCutoffDist()
{
  return cutoff;
}

void CambridgeThreeBodyCorrelationFunction::print(ostream & strm)
{
  strm.unsetf(ios::scientific);
  strm << "Cambridge 3 particle jastrow parameters:" << endl
       << "Nen = " << Nen << endl
       << "Nee = " << Nee << endl
       << "  C = " << C << endl
       << "  L = " << cutoff;
  if(globalInput.flags.optimize_L == 1)
    strm << " (optimized)" << endl;
  else
    strm << " (not optimized)" << endl;

  strm.unsetf(ios::scientific);
  strm << "r1  = x1  / " << (1.0 / L) << endl;
  strm << "r2  = x2  / " << (1.0 / L) << endl;
  strm << "r12 = x12 / " << (1.0 / L) << endl;
  strm << "(r1 - 1)^" << C;
  strm << " (r2 - 1)^" << C;
  strm << " (" << endl;

  //**************************
  bool symmetric = true;
  bool extraPrec = false;
  bool printZero = false;
  //**************************

  int coWidth, coPrec;
  if(extraPrec)
    {
      coPrec  = 15;
      coWidth = 20;
    } else {
      coPrec  = 7;
      coWidth = 10;
    }

  int counter = 0;
  for(int p=0; p<max(Nen,Nee); p++)
    for (int n=0; n<Nee; n++)
      for (int l=0; l<Nen; l++)
	for (int m=0; m<=l; m++)
	  if( (l == p || m == p || n == p) &&
	      (l <= p && m <= p && n <= p))	    
	    {
	      double coeff = coeffs(l,m,n);
	      if(fabs(coeff) < 1e-50 && !printZero)
		continue;
	      
	      if(counter % 4 == 0 && counter != 0)
		strm << endl;
	      counter++;
	      
	      if(fabs(coeff - coeffs(m,l,n)) > 1e-20)
		symmetric = false;
	      
	      stringstream co;
	      co.precision(coPrec);
	      co.width(coWidth);
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
	      
	      int extrawidth = co.str().length() - coWidth;
	      strm << setw(coWidth) << co.str() << setw(31-extrawidth) << left << temp.str();
	    }
  strm << ")";

  if(!symmetric)
    strm << " (not symmetric!!)";

  strm << endl;
  strm << right;
}
