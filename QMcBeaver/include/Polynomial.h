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

#ifndef Polynomial_H
#define Polynomial_H

#include "Array1D.h"
#include "FunctionR1toR1.h"
#include "Complex.h"
#include "Exception.h"

/**
  A one dimensional real polynomial.
  \f[
  P(x) = \sum^{n}_{i=0}c_{i}x^{i}
  \f]
  */

class Polynomial : public FunctionR1toR1
{
private:
  /**
    Coefficients of the polynomial.  They are arranged from smallest to
    largest powers of x.
    */
  Array1D<double> coefficients;

  /**
    Coefficients of the polynomial's first derivative.  They are arranged
    from smallest to largest powers of x.
    */
  Array1D<double> firstDerivativeCoefficients;

  /**
    Coefficients of the polynomial's second derivative.  They are arranged
    from smallest to largest powers of x.
    */
  Array1D<double> secondDerivativeCoefficients;

  /**
    Last calculated function value.
    */
  double f;

  /**
    Last calculated function first derivative value.
    */
  double df;

  /**
    Last calculated function second derivative value.
    */
  double d2f;

  /**
    Last evaluated x value.
    */
  double x;

  /**
    Has the function value been calculated for the last evaluated
    x value.
    */
  bool evaluatedF;

  /**
    Has the function first derivative value been calculated for the 
    last evaluated x value.
    */
  bool evaluatedDF;

  /**
    Has the function second derivative value been calculated for the 
    last evaluated x value.
    */
  bool evaluatedD2F;

public:
  /**
    Constructs an uninitialized instance of this class.
    */
  Polynomial();

  /**
    Constructs and initializes an intance of this class.

    @param coeffs set of polynomial coefficients to use for the polynomial.
    */
  Polynomial(Array1D<double> & coeffs);

  /**
    Initializes this object.

    @param coeffs set of polynomial coefficients to use for the polynomial.
    */
  void initialize(Array1D<double> & coeffs);


  void evaluate(double x);
  double getFunctionValue();
  double getFirstDerivativeValue();
  double getSecondDerivativeValue();

  /**
    Gets the roots of the polynomial.

    @return roots of the polynomial.

    @throws Exception if problems were encounted during the root calculation.
    */
  Array1D<Complex> getRoots();

protected:
  /**
    Gets the number of coefficients in the polynomial.  This is one larger 
    than the order of the polynomial.

    @return number of coefficients in the polynomial.
    */
  int getNumberCoefficients();
  
  /**
    Gets the ith coefficient of the polynomial.  Where the polynomial is 
    defined such that
    \f[
    P(x) = \sum^{n}_{i=0}c_{i}x^{i}
    \f]
    where \f$n\f$ is the order of the polynomial and \f$c_{i}\f$ is the
    ith coefficient.

    @param i index of the coefficient to return.

    @return ith coefficient of the polynomial.
    */
  double getCoefficient(int i);


private:
  /**
    Initialize all of the member variables of this object.
    */
  void initialize();

  /**
    Evaluates a polynomial at \f$x\f$ defined by the entered coefficients.

    @param x point to evaluate the polynomial at.
    @param coeffs coefficients describing the polynomial being evaluated.
    These are ordered from smallest to largest powers of \f$x\f$.

    @return value of the polynomial defined by <code>coeffs</code> at 
    <code>x</code>.
    */
  double evaluate(double x, Array1D<double> & coeffs);


  /**
    Evaluates a polynomial at \f$x\f$ defined by the entered coefficients.
    In addition, it evaluates the 1st two derivatives of that polynomial.
    It stores the results in f, df, and d2f from this class.

    @param x point to evaluate the polynomial at.
    @param coeffs coefficients describing the polynomial being evaluated.
    These are ordered from smallest to largest powers of \f$x\f$.
    */
  void evaluateAll(double x, Array1D<double> & coeffs);

  /**
    Given a complex point \f$x\f$ and a polynomial with complex coefficients,
    converges \f$x\f$ to the root of the polynomial, within achievable 
    roundoff limits, using Laguerre's method.  This method can be used to
    refine the roots of a polynomial with either real or complex coefficients.
    This was adapted from "Numerical Recipes in C".

    @param a complex coefficients of the polynomial ordered from lowest
    to highest powers of x.
    @param m degree of the polynomial.  This is one less than the number of 
    coefficients.
    @param x complex point to refine to the polynomial's root.
    @param its number of iterations taken to refine the root.
    @param calcOK <code>true</code> if there were no problems during the 
    calculation and <code>false</code> if too many iterations were
    performed.  When this happens try a different starting guess for the root.
    */
  void laguer(Array1D<Complex> &a, int m, Complex &x, int *its, bool *calcOK);

  /**
    Finds the roots of a polynomial with complex coefficients using
    Laguerre's method.  This was adapted from "Numerical Recipes in C".

    @param a complex coefficients of the polynomial ordered from lowest
    to highest powers of x.
    @param polish <code>true</code> if the roots are going to be polished
    using Laguerre's method and <code>false</code> if these roots will be
    given to another routine to polish.
    @param calcOK <code>true</code> if there were no problems during the 
    calculation and <code>false</code> if too many iterations were
    performed while calculating one of the roots.  

    @return roots of the polynomial described by a.
    */
  Array1D<Complex> zroots(Array1D<Complex> & a, bool polish, bool *calcOK);
};

#endif



