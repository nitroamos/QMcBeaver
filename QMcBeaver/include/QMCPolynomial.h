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

#ifndef QMCPolynomial_H
#define QMCPolynomial_H

#include "Polynomial.h"

/**
  An extension of Polynomial which adds QMC specific functionality.
  */

class QMCPolynomial : public Polynomial
{

public:
  /**
    Constructs an uninitialized instance of this class.
    */
  QMCPolynomial();

  /**
    Constructs and initializes an intance of this class.

    @param coeffs set of polynomial coefficients to use for the polynomial.
    */
  QMCPolynomial(Array1D<double> & coeffs);



  /**
    Determines if this polynomial has any non-negative real zeroes.

    @return <code>true</code> if the polynomial has a non-negative real zeros 
    and <code>false</code> otherwise.

    @throws Exception if problems were encounted during the calculation.
    */
  bool hasNonNegativeZeroes();
};

#endif
