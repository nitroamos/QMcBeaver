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

#ifndef QMCLineSearchStepLengthSelectionAlgorithm_H
#define QMCLineSearchStepLengthSelectionAlgorithm_H

#include "QMCObjectiveFunction.h"

/**
  Interface to algorithms which determine the proper step length to use
  during a line search optimization (QMCLineSearch).
  */

class QMCLineSearchStepLengthSelectionAlgorithm
{
public:
  /**
    Virtual destructor.
    */
  virtual ~QMCLineSearchStepLengthSelectionAlgorithm(){}

  /**
    Calculates the step length to use when performing a line search 
    optimization.

    @param function objective function being optimized.
    @param position current location of the optimization.
    @param searchDirection direction to optimize along.
    @param gradient current gradient value.
    @param functionValue current function value.
    */
  virtual double stepLength(QMCObjectiveFunction *function, 
			    Array1D<double> & position,
			    Array1D<double> & searchDirection,
			    Array1D<double> & gradient,
			    double functionValue) = 0;

};

#endif
