//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2002.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef QMCOptimizationAlgorithm_H
#define QMCOptimizationAlgorithm_H

#include "Array1D.h"
#include "Array2D.h"

/**
  Interface for numerical optimization algorithms.
  */

class QMCOptimizationAlgorithm
{
public:
  /**
    Virtual destructor.
    */
  virtual ~QMCOptimizationAlgorithm(){};

  /**
    Optimize the function starting from the provided initial guess
    parameters.

    @param initialGuess initial guess parameters for the optimization.
    @return optimized parameters.
    */
  virtual Array1D<double> optimize(Array1D<double> & initialGuess,
				   double value,
				   Array1D<double> & gradient,
				   Array2D<double> & hessian,
				   double) = 0;
};

#endif
