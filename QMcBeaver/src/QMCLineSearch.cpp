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

#include "QMCLineSearch.h"

QMCLineSearch::QMCLineSearch(QMCObjectiveFunction *function, 
		         QMCLineSearchStepLengthSelectionAlgorithm * stepAlg,
			     int MaxSteps, 
			     double tol)
{
  OF            = function;
  stepLengthAlg = stepAlg;
  maximumSteps  = MaxSteps;
  epsilon       = tol;
}

double QMCLineSearch::stepLength(Array1D<double> & x,Array1D<double> & p,
				 Array1D<double> & g, double f)
{
  return stepLengthAlg->stepLength(OF,x,p,g,f); 
}

Array1D<double> QMCLineSearch::optimize(Array1D<double> & InitialGuess)
{
  cout << endl;
  cout << "Beginning Line Search Optimization ... " << endl;

  // Allocate the point used in the search and initialize its value.
  Array1D<double> x = InitialGuess;
  
  // Begin the line search
  double f_old = 0;

  for(int i=0; i<maximumSteps; i++)
    {
      // Calculate the function value, step length, and search direction.
      double f = OF->evaluate(x).getScore();
      
      cout << "\tIteration: " << i << endl;
      cout << "\t\tFunction Value: " << f << endl;
      cout << "\t\tParameters:     " << x << endl;
      
      // if converged quit
      if( i>0 && fabs(1.0-f/f_old) < epsilon ) 
	{
	  cout << "Line Search Optimization Has Converged in " 
	       << i << " Iterations... " << endl; 
	  break;
	}

      Array1D<double> grad_k = getObjectiveFunction()->grad(x);
      Array1D<double> p_k   = searchDirection(x,grad_k);
      double alpha_k = stepLength(x,p_k,grad_k,f);
	
      // Calculate the next step
      for(int j=0; j<x.dim1(); j++)
	{
	  x(j) += alpha_k * p_k(j);
	}
	
      // keep track of the previous function value to know when to finish 
      // the calc
      f_old = f;
    }

  cout << "Ending Line Search Optimization ... " << endl;
  cout << endl;

  return x;
}

QMCObjectiveFunction * QMCLineSearch::getObjectiveFunction()
{
  return OF;
}
