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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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

 private:
  /**
     Evaluates a Chebyschev series
  */
  static double csevl(double x, double coef[], int lengthCoef);

  /*
    Returns the value of x with the sign of y.
  */
  static double sign(double x, double y);
};

#endif
