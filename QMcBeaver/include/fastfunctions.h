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

/**
  @file fastfunctions.h
  This is a fast function library originally
  intended to speed up QMcBeaver a Quantum
  Monte Carlo program. 
  */

#ifndef fastfunctions_H
#define fastfunctions_H

#include <math.h>

/**
  Fast power function for use when the exponent is a small integer.

  @param x base 
  @param n exponent
  @return \f$x^n\f$
  */
double fastPower(double x, int n); 

/**
   A safe implementation of pythagoran's
   theorem.

   It shouldn't under or over flow.
*/
double pythag(double a, double b);

#endif
