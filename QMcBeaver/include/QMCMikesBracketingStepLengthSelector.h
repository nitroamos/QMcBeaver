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

#ifndef QMCMikesBracketingStepLengthSelector_H
#define QMCMikesBracketingStepLengthSelector_H

#include <iostream>

#include "QMCLineSearchStepLengthSelectionAlgorithm.h"

using namespace std;

/**
  Algorithm to determine the step length for a line search optimization
  developed by Michael Todd Feldmann.  This algorithm is purely huristic
  and does not insure the Wolfe conditions or other such properties.  
  Again, much work could be done to do this part of a line search better.
  */

class QMCMikesBracketingStepLengthSelector : 
public QMCLineSearchStepLengthSelectionAlgorithm
{
public:
  double stepLength(QMCObjectiveFunction *function, 
		    Array1D<double> & position,
		    Array1D<double> & searchDirection,
		    Array1D<double> & gradient,
		    double functionValue);

private:
  void bracket(QMCObjectiveFunction *OF,Array1D<double> & position,
	       Array1D<double> & searchDirection,
	       double alpha_zero,double scale,int recursion_depth);

  void quadratic_rebracketer(QMCObjectiveFunction * OF,
			     Array1D<double> & position,
			     Array1D<double> & searchDirection);

  void print_interval();
  
  //The next six variables make up an interval to do a parabolic fit to.
  //Think of the alphas as "x" and the scores as "y" to find the best "x"
  //which will be the new alpha_middle
  double alpha_left;      
  double alpha_middle;
  double alpha_right;
  double score_left;
  double score_middle;
  double score_right;


};

#endif
