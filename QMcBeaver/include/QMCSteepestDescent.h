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

#ifndef QMCSteepestDescent_H
#define QMCSteepestDescent_H

#include "QMCLineSearch.h"

/**
  Steepest descent line search numerical optimization algorithm.
  As is standard in the field, the optimization is a minimization.
  */

class QMCSteepestDescent : public QMCLineSearch
{
public:
  /**
    Constructs and initializes an instance of this class.

    @param function objective function to optimize.
    @param stepAlg algorithm to use in determining the line search step length.
    @param maxSteps maximum number of steps to be performed during the line 
    search.
    @param tol tolerance to converge the solution to.  Calculation is converged
    when \f$\left| 1-\frac{f(x_{i+1})}{f(x_{i})} \right| < tol \f$.
    */
  QMCSteepestDescent(QMCObjectiveFunction * function, 
		     QMCLineSearchStepLengthSelectionAlgorithm * stepAlg, 
		     int maxSteps, double tol);


private:

  /**
    Calculates the search direction at x. 
    x' = x + StepLength * SearchDirection.
    Where the search direction is the negative of the gradient.
    */
  Array1D<double> searchDirection(Array1D<double> & x,Array1D<double> & g);
};

#endif
