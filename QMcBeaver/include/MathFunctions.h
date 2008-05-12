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


/**
   A set of basic mathematical functions.
*/

class MathFunctions
{
 public:

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
