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

#include "QMCPolynomial.h"

QMCPolynomial::QMCPolynomial():Polynomial()
{
}

QMCPolynomial::QMCPolynomial(Array1D<double> & coeffs):Polynomial(coeffs)
{
}

bool QMCPolynomial::hasNonNegativeZeroes()
{
  const double tol = 1e-8;

  Array1D<Complex> roots = getRoots();

  for(int i=0; i<roots.dim1(); i++)
    {
      double re = roots(i).real();
      double im = roots(i).imaginary();
      
      if( re >= 0 && fabs(im) < tol ) return true; 
    }

  return false;
}




