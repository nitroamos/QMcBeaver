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

#ifndef QMCObjectiveFunction_H
#define QMCObjectiveFunction_H

#include <iostream>
#include <string>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCObjectiveFunctionResult.h"
#include "QMCProperties.h"
#include "QMCReadAndEvaluateConfigs.h"

using namespace std;

/**
  Objective function optimized during a variational QMC (VMC) calculation
  to find the optimal wavefunction parameters.  As is standard in the
  field of numerical optimization, optimization means minimization.  The
  particular form of the objective function is determined by parameters in
  the input file.
  */

class QMCObjectiveFunction
{
 public:  
  /**
    Initializes this object.  This must be called before any other functions
    in this object are called.

    @param input input data for the calculation
    */
  void initialize(QMCInput *input, int configsToSkip);

  /**
    Evaluates and returns the result of the objective function 
    evaluated with a single set of parameters.

    @param params set of parameters to evaluate the objective function with.

    @return result of the objective function evaluation.
    */
  QMCObjectiveFunctionResult evaluate(Array1D<double> &params);

  /**
    Evaluates and returns the result of the objective function 
    evaluated with multiple single sets of parameters.

    @param params sets of parameters to evaluate the objective function with.

    @return results of the objective function evaluations.  The index of the 
    input parameters corresponds to the index of the returned values.
    */
  Array1D< QMCObjectiveFunctionResult > evaluate(Array1D< Array1D<double> > &
						 params);
  
  /**
    Evaluates and returns the gradient of the objective function
    for one set of parameters.

    @param params sets of parameters to evaluate the gradient with.

    @return gradient of the objective function.
    */
  Array1D<double> grad(Array1D<double> & params);

  /**
    Evaluates and returns the gradient of the objective function 
    for multiple sets of parameters.

    @param params sets of parameters to evaluate the gradient with.

    @return gradients of the objective function.  The index of the 
    input parameters corresponds to the index of the returned values.
    */
  Array1D< Array1D<double> > grad(Array1D< Array1D<double> > & params);
 
 private:
  QMCInput *Input;

  QMCReadAndEvaluateConfigs RAEC;
  
  void numerical_grad(Array1D<double> & Params, Array1D<double> & GRAD);
};

#endif
