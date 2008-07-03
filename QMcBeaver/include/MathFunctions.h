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

#ifndef MathFunctions_H
#define MathFunctions_H

#include <stdio.h>
#include <stdlib.h>
#include "IeeeMath.h"
#include <math.h>
#include "Array2D.h"

/**
   A set of basic mathematical functions.
*/

class MathFunctions
{
 public:

  /**
     Legendre polynomial order l at x.
  */
  static double legendre(int l, double x);

  /**
     Returns the error function of a double.
     @param x A double value.
     @return The error function of x.
  */
  static double erf(double x);

  /**
     Returns the complementary error function of a double.
     @param x A double value.
     @return The complementary error function of x.
  */
  static double erfc(double x);

  static double F_gamma(double v, double t);

  static double gamma_inc(double a, double x);

  static double gamma_log(double x2);

  static double gamma_series(double a, double x);

  static double gamma_cf(double a, double x);

  /**
     Used for Umrigar's spherical metropolis sampling.

     Look in A.4 of 
     "Variational Monte Carlo Basics and Applications to Atoms and Molecules"
     by C. J. Umrigar
     NATO ASI Series, Series C, Mathematical and Physical Sciences, Vol. C-525,
     (Kluwer Academic Publishers, Boston, 1999)
  */
  static double accel_G(double n, double x, double y,
			double a, double zeta);
  static void fitU(double Z, double Frk, double r, double delta,
		   double & zeta, double & a, double & I);

  /**
     Calculates distance between electrons i and j
     based on the coordinates found in R
  */
  static double rij(Array2D<double> &R, int i, int j);
  
  /**
     Calculates the distance between particles i and j,
     where i is found in positioni and
     j is found in positionj.
  */
  static double rij(Array2D<double> &positioni,
		    Array2D<double> &positionj, 
		    int i, int j);

  /**
     This function was copied from:
     On a Fast, Compact Approximation of the Exponential Function
     by Gavin C Cawley in Neural Computation 12, 2009-2012 (2000)

     which is itself based on:
     A Fast, Compact Approximation of the Exponential Function
     by Nicol N Schraudolph in Neural Computation 11, 853-862 (1999)

     Essentially, it's taking advantage of the fact that the exponent term
     in a floating point double is 2^x and is only different from exp(x) by a change of base
     formula. The function is essentially free, but you get ~1% error.
  */
  static double fast_exp(double x);

  static const double pi = 3.141592653589793;
 private:
  /**
     Evaluates a Chebyschev series
  */
  static double csevl(double x, double * coef, int lengthCoef);

  /*
    Returns the value of x with the sign of y.
  */
  static double sign(double x, double y);
};

#endif
