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

#ifndef QMCOptimizationFactory_H
#define QMCOptimizationFactory_H

#include "QMCOptimizationAlgorithm.h"
#include "QMCObjectiveFunction.h"
#include "QMCInput.h"

#include "QMCLineSearchStepLengthSelectionFactory.h"
#include "QMCSteepestDescent.h"
#include "QMCBFGSQuasiNewtonLineSearch.h"
#include "CKGeneticAlgorithm1.h"

/**
  Object factory which returns the correct QMCOptimizationAlgorithm
  specified in the calculation input data.  Optimization assumed to mean
  minimization, as is standard in the field.
  */

class QMCOptimizationFactory
{
  public:
  /** 
    Returns the correct QMCOptimizationAlgorithm specified in the calculation
    input data.
  
    @param objFunc object function to optimize.
    @param input input data to control the calculation.
    */

  static QMCOptimizationAlgorithm * 
    optimizationAlgorithmFactory(QMCObjectiveFunction &objFunc,
				 QMCInput * input);

};

#endif
