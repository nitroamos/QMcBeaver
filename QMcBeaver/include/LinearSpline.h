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

#ifndef LINEAR_SPLINE_H
#define LINEAR_SPLINE_H

#include <iostream>
#include <string>
#include <math.h>

#include "Array1D.h"
#include "FunctionR1toR1.h"

using namespace std;

/**
   Very crude linear spline class.
*/

class LinearSpline : FunctionR1toR1
{
public:
  LinearSpline(){;}                  //constructs a Spline_Linear

  void operator=( const LinearSpline & rhs );

  /** 
  Initialize the spline with the function values at given points.
  */
  void initializeWithFunctionValues(Array1D<double> &x_input, 
				    Array1D<double> &y_input);

  void evaluate(double x);

  double getFunctionValue();  

  double getFirstDerivativeValue();

  double getSecondDerivativeValue();


private:
  Array1D<double> x_list;      //domain values
  Array1D<double> y_list;      //function values for each x
  Array1D<double> a0_list;     //poly coeff of f at each point
  Array1D<double> a1_list;     //poly coeff of f at each point
     
  int n;                       //number of input data points;

  double f;
  double df;
};


#endif




