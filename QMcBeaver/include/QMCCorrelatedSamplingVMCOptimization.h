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

#ifndef QMCCorrelatedSamplingVMCOptimization_H
#define QMCCorrelatedSamplingVMCOptimization_H

#include "QMCInput.h"
#include "QMCObjectiveFunction.h"
#include "QMCOptimizationFactory.h"
#include "QMCReadAndEvaluateConfigs.h"
#include "QMCProperties.h"
#include "QMCDerivativeProperties.h"


/**
  Optimize the parameters in a variational QMC (VMC) calculation using
  the correlated sampling method.
*/

class QMCCorrelatedSamplingVMCOptimization
{
 public:
  /**
     Optimizes the parameters in a variational QMC (VMC) calculation using
     the correlated sampling method.

     @param input data input to control the calculation.
  */
  static void optimize(QMCInput * input,
		       QMCProperties & lastRun,
		       QMCFutureWalkingProperties & fwLastRun,
		       int configsToSkip);
};

#endif

