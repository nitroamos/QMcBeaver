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


void QMCBFGSQuasiNewtonLineSearch::calculateHessian()
{

  int last = gradient.size() - 1;

  /*
    This method will only provide an inverse hessian matrix
    in the second round.

    If we simply return right here, then the method is just
    steepest descent.
  */
  if(last <= 0)
    return;

  Array1D<double> s = x[last]        - x[last-1];
  Array1D<double> y = gradient[last] - gradient[last-1];

  /*
    approximateInverseHessian is just an alias for the inverseHessian
    that we need to update
  */
  Array2D<double> & approximateInverseHessian = inverseHessian[last];

  double rho = 1.0/(s*y);
  
  Array2D<double> tempMatrix1(dim,dim);
  
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
	{
	  tempMatrix1(i,j) = -rho*s(i)*y(j);
	}
      tempMatrix1(i,i) += 1.0;
    }
  
  Array2D<double> tempMatrix2 = tempMatrix1 * approximateInverseHessian;
  
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
	{
	  tempMatrix1(i,j) = -rho*s(j)*y(i);
	}
      tempMatrix1(i,i) += 1.0;
    }
  
  approximateInverseHessian = tempMatrix2 * tempMatrix1;
  
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
	{
	  approximateInverseHessian(i,j) += rho*s(i)*s(j);
	}
    }
}
