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

#ifndef QMCEigenSearch_H
#define QMCEigenSearch_H

#include <math.h>

#include <vector>
#include "QMCObjectiveFunction.h"
#include "QMCOptimizationAlgorithm.h"
#include "QMCLineSearchStepLengthSelectionAlgorithm.h"

/**
  Abstract implementation of a line search numerical optimization algorithm.
  As is standard in the field, the optimization is a minimization.
  */

class QMCEigenSearch : public QMCOptimizationAlgorithm
{
public:
  /**
    Constructs and initializes an instance of this class.

    @param function objective function to optimize.
    @param stepAlg algorithm to use claculate the step length.
    @param maxSteps maximum number of steps to be performed during the line 
    search.
    @param tol tolerance to converge the solution to.  Calculation is converged
    when \f$\left| 1-\frac{f(x_{i+1})}{f(x_{i})} \right| < tol \f$.
    */
  QMCEigenSearch(QMCObjectiveFunction * function, 
		      QMCLineSearchStepLengthSelectionAlgorithm * stepAlg, 
		      int maxSteps, double tol);
  
  /**
    Virtual destructor.
    */
  virtual ~QMCEigenSearch(){};

  Array1D<double> optimize(Array1D<double> & initialGuess,
			   QMCDerivativeProperties & dp,
			   double a_diag_factor,
			   int optStep);

protected:

  /**
    Gets the objective function for the calculation.
    */
  QMCObjectiveFunction * getObjectiveFunction();

  int dim;
  vector<double> f;
  vector< Array1D<double> > x;

private:  
  /**
     Objective function to optimize.
  */
  QMCObjectiveFunction *OF;
  
  /**
     Step length determining algorithm to use.
  */
  QMCLineSearchStepLengthSelectionAlgorithm * stepLengthAlg;
  
  /**
     Tolerance to converge the solution to.
  */
  double epsilon;
  
  /**
     Maximum number of steps allowed for the search.
  */
  int maximumSteps;
};

#endif


