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
  int walkersPerPass = Input->flags.walkers_per_pass;
  sum_U.allocate(walkersPerPass);
  grad_sum_U.allocate(walkersPerPass);
  laplacian_sum_U.allocate(walkersPerPass);

  JastrowElectronNuclear.initialize(Input);
  JastrowElectronElectron.initialize(Input);
}

double QMCJastrow::getJastrow(int which)
{
  return exp(sum_U(which));
}

double QMCJastrow::getLnJastrow(int which)
{
  return sum_U(which);
}

Array2D<double> * QMCJastrow::getGradientLnJastrow(int which)
{
  return &grad_sum_U(which);
}

double QMCJastrow::getLaplacianLnJastrow(int which)
{
  return laplacian_sum_U(which);
}

void QMCJastrow::evaluate(Array1D<Array2D<double>*> &X, int num)
{
  evaluate(Input->JP,X,num);
}

void QMCJastrow::evaluate(QMCJastrowParameters & JP, Array1D<Array2D<double>*> &X, int num)
{
  Array2D<double> * grad_JEN;
  Array2D<double> * grad_JEE;
  for(int walker = 0; walker < num; walker++)
    {
      JastrowElectronNuclear.evaluate(JP,*X(walker));
      JastrowElectronElectron.evaluate(JP,*X(walker));

      sum_U(walker) = JastrowElectronNuclear.getLnJastrow() +
                      JastrowElectronElectron.getLnJastrow();

      laplacian_sum_U(walker) = JastrowElectronNuclear.getLaplacianLnJastrow() 
	+ JastrowElectronElectron.getLaplacianLnJastrow();

      grad_JEN = JastrowElectronNuclear.getGradientLnJastrow();
      grad_JEE = JastrowElectronElectron.getGradientLnJastrow();

      grad_sum_U(walker).allocate(X(walker)->dim1(),3);

      for(int i=0; i<grad_JEE->dim1(); i++)
	for(int j=0; j<grad_JEE->dim2(); j++)
	  (grad_sum_U(walker))(i,j) = (*grad_JEE)(i,j) + (*grad_JEN)(i,j);
    }
}
