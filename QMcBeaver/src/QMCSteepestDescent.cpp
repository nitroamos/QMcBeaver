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

void QMCSteepestDescent::calculateHessian()
{
  /*
    Steepest_Descent is merely the Conjugate Gradient
    method where the identity matrix is used as the
    inverse Hessian.

    The inverseHessian was initialized
    to be the identity matrix.

    Therefore, there's nothing for us to do here.
  */
}
