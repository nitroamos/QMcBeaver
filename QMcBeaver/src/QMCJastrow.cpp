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

#include "QMCJastrow.h"

void QMCJastrow::initialize(QMCInput * input)
{
  Input = input;

  JastrowElectronNuclear.initialize(Input);
  JastrowElectronElectron.initialize(Input);
}

double QMCJastrow::getJastrow()
{
  return exp(sum_U);
}

double QMCJastrow::getLnJastrow()
{
  return sum_U;
}

Array2D<double> * QMCJastrow::getGradientLnJastrow()
{
  return &grad_sum_U;
}

double QMCJastrow::getLaplacianLnJastrow()
{
  return laplacian_sum_U;
}

void QMCJastrow::evaluate(Array2D < double > & X)
{
  evaluate(Input->JP,X);
}

void QMCJastrow::evaluate(QMCJastrowParameters & JP, Array2D<double> & X)
{
  JastrowElectronNuclear.evaluate(JP,X);
  JastrowElectronElectron.evaluate(JP,X);

  sum_U = JastrowElectronNuclear.getLnJastrow() + 
    JastrowElectronElectron.getLnJastrow();

  laplacian_sum_U = JastrowElectronNuclear.getLaplacianLnJastrow() +
    JastrowElectronElectron.getLaplacianLnJastrow();

  Array2D<double> * grad_JEN = JastrowElectronNuclear.getGradientLnJastrow();
  Array2D<double> * grad_JEE = JastrowElectronElectron.getGradientLnJastrow();

  grad_sum_U.allocate(X.dim1(),3);

  for(int i=0; i<grad_JEE->dim1(); i++)
    {
      for(int j=0; j<grad_JEE->dim2(); j++)
	{
	  grad_sum_U(i,j) = (*grad_JEE)(i,j) + 
	    (*grad_JEN)(i,j);
	}
    }
}
