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

#include "QMCSteepestDescent.h"

QMCSteepestDescent::QMCSteepestDescent(QMCObjectiveFunction * function, 
	 QMCLineSearchStepLengthSelectionAlgorithm * stepAlg,
	 int maxSteps, double tol):QMCLineSearch(function,stepAlg,maxSteps,tol)
{
}


Array1D<double> QMCSteepestDescent::searchDirection(Array1D<double> & x,
						    Array1D<double> & g)
{
  return g*-1.0;
}
