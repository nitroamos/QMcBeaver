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

#include "QMCBFGSQuasiNewtonLineSearch.h"

QMCBFGSQuasiNewtonLineSearch::QMCBFGSQuasiNewtonLineSearch(
         QMCObjectiveFunction * function, 
	 QMCLineSearchStepLengthSelectionAlgorithm * stepAlg,
	 int maxSteps, double tol):QMCLineSearch(function,stepAlg,maxSteps,tol)
{
}


Array1D<double> QMCBFGSQuasiNewtonLineSearch::searchDirection(
						    Array1D<double> & x,
						    Array1D<double> & gradient)
{
  const int dimension = x.dim1();

  // update the approximate inverse Hessian

  if( approximateInverseHessian.dim1() != dimension || 
      approximateInverseHessian.dim2() != dimension )
    {
      // The inverse hessian needs to be initialized

      approximateInverseHessian.allocate(dimension,dimension);

      // as an initial guess, use the identity matrix

      approximateInverseHessian = 0.0;

      for(int i=0; i<dimension; i++)
	{
	  approximateInverseHessian(i,i) = 1.0;
	}
    }
  else
    {
      // perform a normal update to the hessian

      Array1D<double> s = x - xOld;
      Array1D<double> y = gradient - gradientOld;
      double rho = 1.0/(s*y);

      Array2D<double> tempMatrix1(dimension,dimension);

      for(int i=0; i<dimension; i++)
	{
	  for(int j=0; j<dimension; j++)
	    {
	      tempMatrix1(i,j) = -rho*s(i)*y(j);
	    }
	  tempMatrix1(i,i) += 1.0;
	}
      
      Array2D<double> tempMatrix2 = tempMatrix1 * approximateInverseHessian;

      for(int i=0; i<dimension; i++)
	{
	  for(int j=0; j<dimension; j++)
	    {
	      tempMatrix1(i,j) = -rho*s(j)*y(i);
	    }
	  tempMatrix1(i,i) += 1.0;
	}

      approximateInverseHessian = tempMatrix2 * tempMatrix1;

      for(int i=0; i<dimension; i++)
	{
	  for(int j=0; j<dimension; j++)
	    {
	      approximateInverseHessian(i,j) += rho*s(i)*s(j);
	    }
	}
    }

  // calculate the search direction, p_k

  Array1D<double> p_k(dimension);

  for(int i=0; i<dimension; i++)
    {
      p_k(i) = 0.0;
      for(int j=0; j<dimension; j++)
	{
	  p_k(i) -= approximateInverseHessian(i,j) * gradient(j);
	}
    }

  // set the old values of the gradient and position equal to the current
  // values

  xOld = x;
  gradientOld = gradient;

  // return the search direction

  return p_k;
}
