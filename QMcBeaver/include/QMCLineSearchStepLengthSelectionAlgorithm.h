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
    
     You will have to look at the individual algorithms to
     see how the input parameters are defined.
  */
  virtual double stepLength(QMCObjectiveFunction *function, 
			    Array1D<double> & array1,
			    Array1D<double> & array2,
			    Array1D<double> & array3,
			    Array2D<double> & matrix1,
			    double scalar1) = 0;

};

#endif
