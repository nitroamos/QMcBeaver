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

#ifndef FunctionR1toR1_H
#define FunctionR1toR1_H

#include "Array1D.h"

using namespace std;

/**
  An interface for a function from \f$\mathbf{R}^{1} \rightarrow 
  \mathbf{R}^{1}\f$.
*/

class FunctionR1toR1
{
public:
  /**
    Virtual destructor.
    */
  virtual ~FunctionR1toR1(){};

  /**
     Evaluates the function at \f$x\f$.

     @param x point to evaluate the function.
  */

  virtual void evaluate(double x) = 0;


  /**
     Gets the function value at the last evaluated point.

     @return function value.
  */

  virtual double getFunctionValue() = 0;


  /**
     Gets the function's first deriviate at the last evaluated point.

     @return function's deriviative value.
  */

  virtual double getFirstDerivativeValue() = 0;


  /** 
      Gets the function's second deriviative at the last evaluated point.

      @return function's second derivative value.
  */
  virtual double getSecondDerivativeValue() = 0;
};


#endif




