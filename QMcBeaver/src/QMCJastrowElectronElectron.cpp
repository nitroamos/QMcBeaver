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

#include "QMCJastrowElectronElectron.h"

void QMCJastrowElectronElectron::initialize(QMCInput * input)
{
  Input = input;
}

/**
* Find the unit vector and distance between X1 and X2.  The unit vector is in
 * the direction of X1-X2.
 */

void QMCJastrowElectronElectron::calculateDistanceAndUnitVector(
  Array2D<double> & X1, int x1particle, Array2D<double> &X2,
  int x2particle, double & r, Array1D<double> & UnitVector)
{
  double r_sq = 0;

  for(int i=0; i<3; i++)
    {
      UnitVector(i) = X1(x1particle,i) - X2(x2particle,i);
      r_sq += UnitVector(i) * UnitVector(i);
    }

  r = sqrt( r_sq );

  UnitVector *= 1.0/r;
}


double QMCJastrowElectronElectron::getLaplacianLnJastrow()
{
  return laplacian_sum_U;
}

Array2D<double> * QMCJastrowElectronElectron::getGradientLnJastrow()
{
  return &grad_sum_U;
}

double QMCJastrowElectronElectron::getLnJastrow()
{
  return sum_U;
}

void QMCJastrowElectronElectron::evaluate(QMCJastrowParameters & JP,
    Array2D<double> & X)
{
  // initialize the results

  sum_U = 0.0;
  laplacian_sum_U = 0.0;
  grad_sum_U.allocate(X.dim1(),3);
  grad_sum_U = 0.0;

  // Get values from JP that will be needed during the calc

  QMCCorrelationFunctionParameters * EupEdn = 0;

  if( Input->WF.getNumberAlphaElectrons() > 0 &&
      Input->WF.getNumberBetaElectrons() > 0 )
    {
      EupEdn = JP.getElectronUpElectronDownParameters();
    }

  QMCCorrelationFunctionParameters * EupEup = 0;

  if( Input->WF.getNumberAlphaElectrons() > 1 )
    {
      EupEup = JP.getElectronUpElectronUpParameters();
    }

  QMCCorrelationFunctionParameters * EdnEdn = 0;

  if( Input->WF.getNumberBetaElectrons() > 1 )
    {
      EdnEdn = JP.getElectronDownElectronDownParameters();
    }

  // Loop over each electron calculating the e-e jastrow function
  //Array1D<double> UnitVector;
  int numAlpha = Input->WF.getNumberAlphaElectrons();

  // Get the correct correlation function to use and evaluate it
  //I separated the collectForPair so that the inner loop didn't have
  //a giant if statement.
  QMCCorrelationFunction *U_Function = 0;
  if(EupEup != 0)
    U_Function = EupEup->getCorrelationFunction();
  for(int Electron1=0; Electron1<numAlpha; Electron1++)
    for(int Electron2=0; Electron2<Electron1; Electron2++)
      collectForPair(Electron1,Electron2,U_Function,X);

  if(EdnEdn != 0)
    U_Function = EdnEdn->getCorrelationFunction();
  for(int Electron1=numAlpha; Electron1<X.dim1(); Electron1++)
    for(int Electron2=numAlpha; Electron2<Electron1; Electron2++)
      collectForPair(Electron1,Electron2,U_Function,X);

  if(EupEdn != 0)
    U_Function = EupEdn->getCorrelationFunction();
  for(int Electron1=0; Electron1<numAlpha; Electron1++)
    for(int Electron2=numAlpha; Electron2<X.dim1(); Electron2++)
      collectForPair(Electron1,Electron2,U_Function,X);

}

inline void QMCJastrowElectronElectron::collectForPair(int Electron1, int Electron2,
    QMCCorrelationFunction *U_Function, Array2D<double> & X)
{
  // Find the unit vector between electron1 and electron2 and their
  // distance apart
  double * unitVector = new double[3];
  double r = 0;
  double firstDeriv;
  double temp;

  for(int i=0; i<3; i++)
    {
      unitVector[i] = X(Electron1,i) - X(Electron2,i);
      r += unitVector[i]*unitVector[i];
    }

  r = sqrt( r );

  U_Function->evaluate(r);

  // Update the values being calculated ...

  sum_U +=  U_Function->getFunctionValue();
  firstDeriv = U_Function->getFirstDerivativeValue();
  laplacian_sum_U += 2.0*(2.0/r *
                          firstDeriv +
                          U_Function->getSecondDerivativeValue());

  for(int i=0; i<3; i++)
    {
      unitVector[i] /= r;
      temp = firstDeriv * unitVector[i];

      grad_sum_U(Electron1,i) += temp;
      grad_sum_U(Electron2,i) -= temp;
    }

  delete [] unitVector;
  unitVector = 0;
}
