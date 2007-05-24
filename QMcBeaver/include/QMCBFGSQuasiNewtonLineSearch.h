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

#ifndef QMCBFGSQuasiNewtonLineSearch_H
#define QMCBFGSQuasiNewtonLineSearch_H

#include "QMCLineSearch.h"

/**
  BFGS Quasi-Newton line search numerical optimization algorithm.
  As is standard in the field, the optimization is a minimization.
  */

class QMCBFGSQuasiNewtonLineSearch : public QMCLineSearch
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
  QMCBFGSQuasiNewtonLineSearch(QMCObjectiveFunction * function, 
		     QMCLineSearchStepLengthSelectionAlgorithm * stepAlg, 
		     int maxSteps, double tol);


private:

  /**
    Calculates the search direction at x. 
    x' = x + StepLength * SearchDirection.
    Where the search direction is the approximately the Newton step 
    direction.  The step is determined using the BFGS Quasi-Newton algorithm.
    */
  void calculateHessian();
};

#endif
