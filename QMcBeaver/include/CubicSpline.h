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

#ifndef CubicSpline_H
#define CubicSpline_H

#include <iostream>
#include <string>
#include <math.h>

#include "Array1D.h"
#include "FunctionR1toR1.h"

using namespace std;

/**
  A 1-dimensional (\f$\mathbf{R}^{1} \rightarrow \mathbf{R}^{1}\f$) cubic
  spline interpolation.
*/

class CubicSpline : public FunctionR1toR1
{
public:
  /**
     Creates an instance of this class.
  */

  CubicSpline();


  /**
     Sets two CubicSpline objects equal.

    @param rhs object to set this object equal to
  */

  void operator=( const CubicSpline & rhs );


  /**
     Initializes the spline with the function values at given points plus
     the derivative values at the end points.

     @param xInput x values of the given points.
     @param yInput y values of the given points.
     @param yPrimeFirst derivative value at the first point.
     @param yPrimeLast derivative value at the last point.
  */

  void initializeWithFunctionValues(const Array1D<double> &xInput, 
				    const Array1D<double> &yInput, 
				    double yPrimeFirst, double yPrimeLast);


  /**
     Initializes the spline with the derivative values at given points plus
     the function value at the first point.

     @param xInput x values of the given points.
     @param yPrimeInput derivative values of the given points.
     @param yFirst function value at the first point.
  */

  void initializeWithDerivativeValues(const Array1D<double> &xInput, 
				      const Array1D<double> &yPrimeInput, 
				      double yFirst);

  void evaluate(double x);
  double getFunctionValue();
  double getFirstDerivativeValue();
  double getSecondDerivativeValue();

  /**
     This function will integrate the spline between two
     indices.

     @param indexStart
     @param indexStop
    */
  
  double integrate(int indexStart, int indexStop);

  /**
    Writes the state of this object to an XML stream.

    @param strm XML stream
    */

  void toXML(ostream& strm);

  /**
    Writes the state of this object to an XML stream.
    */

  friend ostream& operator <<(ostream& strm, CubicSpline & rhs);

protected:
  /**
    Evaluate the function at \f$x\f$ when the index of the box of the domain
    containing \f$x\f$ is known.

    @param x point to evaluate the function.
    @param index index of the box of the domain containing <code>x</code>.
    */
  void evaluate(double x, int index);

private:
  Array1D<double> x_list;      //domain values
  Array1D<double> y_list;      //function values for each x
  Array1D<double> yp_list;     //1st derv of y at each point
  Array1D<double> a0_list;     //poly coeff of f at each point
  Array1D<double> a1_list;     //poly coeff of f at each point
  Array1D<double> a2_list;     //poly coeff of f at each point
  Array1D<double> a3_list;     //poly coeff of f at each point
     
  Array1D<double> row_right;
  Array1D<double> row_middle;
  Array1D<double> row_left;
  Array1D<double> col_rhs;
                  
  double y0;                   //y for smallest x
  double yp0;                  //y prime for smallest x
  double ypend;                //y prime for largest x
  int n;                       //number of input data points;

  double f;
  double dfdx;
  double ddfddx;

  void solve_tridiagonal_system();
  void solve_tridiagonal_system2();
};


#endif




