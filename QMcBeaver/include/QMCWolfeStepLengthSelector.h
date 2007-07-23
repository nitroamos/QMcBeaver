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

#ifndef QMCWolfeStepLengthSelector_H
#define QMCWolfeStepLengthSelector_H

#include <iostream>

#include "QMCLineSearchStepLengthSelectionAlgorithm.h"

using namespace std;

/**
  Algorithm to determine the step length for a line search optimization
  so that the step length satisfies the Wolfe conditions.  This algorithm
  is described in Nocedal and Wright.
  */

class QMCWolfeStepLengthSelector : 
public QMCLineSearchStepLengthSelectionAlgorithm
{
public:
  double stepLength(QMCObjectiveFunction *function, 
		    Array1D<double> & position,
		    Array1D<double> & searchDirection,
		    Array1D<double> & gradient,
		    Array2D<double> & unused,
		    double functionValue);

private:
  /**
     The objective function used in calculating the step length.
  */
  QMCObjectiveFunction *OF;

    /*
    Calculates the value of the 1-D objective function from the line search.
    This function is 
    \f[
    \phi(\alpha) = r^2(p+\alpha d)
    \f]
    where \f$p\f$ is the current set of parameters, \f$d\f$ is the search
    direction, and \f$\alpha\f$ is the step length parameter.
 
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param stepLength parameter determining how long of a step to take in the
    line search.

    @return value of the 1-D objective function for the given step length 
    parameter.
    */
  double calculateLineSearchObjectiveFunction(Array1D<double> & params,
                                  Array1D<double> & searchDirection,
						     double stepLength);

  /*
    Calculates the derivative of the 1-D objective function from the line 
    search.  This function is 
    \f[
    \phi(\alpha) = r^2(p+\alpha d)
    \f]
    where \f$p\f$ is the current set of parameters, \f$d\f$ is the search
    direction, and \f$\alpha\f$ is the step length parameter.
 
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param stepLength parameter determining how long of a step to take in the
    line search.

    @return derivate of the 1-D objective function for the given step length 
    parameter.
    */
  double calculateLineSearchObjectiveFunctionDerivative(
                                  Array1D<double> & params,
                                  Array1D<double> & searchDirection,
                                  double stepLength);
  /**
    Perform a cubic interpolation to choose a step length in the 
    interval [a_lo,a_hi]. See Nocedal and Wright p 56-58.

    @param alphaLo lower bound on the step length.
    @param alphaHi upper bound on the step length.
    @param phi_0 objective function value at the current parameters
    @param phi_alphaLo 1-D objective function value at the smallest step 
    length.
    @param phi_alphaHi 1-D objective function value at the largest step length.
    @param phiPrime_0 derivative of the 1-D objective function at the
    current parameters.

    @return calculated step length 
    */
  double cubicInterpolateStep(double alphaLo, double alphaHi, 
                                     double phi_0, 
                                     double phi_alphaLo, double phi_alphaHi, 
                                     double phiPrime_0);

  /**
    Finds a point in the interval [a_lo, a_hi] that satisfy the strong Wolfe
    conditions.  See Nocedal and Wright p 60.

    @param alphaLo lower bound on the step length.
    @param alphaHi upper bound on the step length.
    @param phi_0 objective function value at the current parameters
    @param phi_alphaLo 1-D objective function value at the smallest step 
    length.
    @param phi_alphaHi 1-D objective function value at the largest step length.
    @param phiPrime_0 derivative of the 1-D objective function at the
    current parameters.
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.

    @return step length satisfying the strong Wolfe conditions.
    */
  double zoom(double alphaLo, double alphaHi, double phi_0, 
                     double phi_alphaLo, double phi_alphaHi, 
                     double phiPrime_0, Array1D<double> & params,
                     Array1D<double> & searchDirection);

  /**
    Finds a point between 0 and the maximum step length which satisfies the
    Wolfe conditions.  See Nocedal and Wright p 59.

    @param alphaGuess guess at a step length satisfying the Wolfe conditions.
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param gradient gradient of the objective function.
    @param functionValue value of the objective function at the current 
    set of parameters.

    @return step length satisfying the strong Wolfe conditions.
    */
  double wolfeStepLength(double alphaGuess, Array1D<double> & params,
                                Array1D<double> & searchDirection,
                                Array1D<double> & gradient,
                                double functionValue);
 

};

#endif
