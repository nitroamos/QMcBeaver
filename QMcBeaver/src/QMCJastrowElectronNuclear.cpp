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

#include "QMCJastrowElectronNuclear.h"

void QMCJastrowElectronNuclear::initialize(QMCInput * input)
{
  Input = input;
}

/**
* Find the unit vector and distance between X1 and X2.  The unit vector is in
 * the direction of X1-X2.
 */

void QMCJastrowElectronNuclear::calculateDistanceAndUnitVector(
  Array2D<double> & X1, int x1particle, Array2D<double> &X2,
  int x2particle, double & r, Array1D<double> & UnitVector)
{
  double r_sq = 0;

  for(int i=0; i<3; i++)
    {
      UnitVector(i) =
        X1(x1particle,i) - X2(x2particle,i);

      r_sq += UnitVector(i) * UnitVector(i);
    }

  r = sqrt( r_sq );

  UnitVector *= 1.0/r;
}


double QMCJastrowElectronNuclear::getLaplacianLnJastrow()
{
  return laplacian_sum_U;
}

Array2D<double> * QMCJastrowElectronNuclear::getGradientLnJastrow()
{
  return &grad_sum_U;
}

double QMCJastrowElectronNuclear::getLnJastrow()
{
  return sum_U;
}

void QMCJastrowElectronNuclear::evaluate(QMCJastrowParameters & JP,
    Array2D<double> & X)
{
  // initialize the results

  sum_U = 0.0;
  laplacian_sum_U = 0.0;
  grad_sum_U.allocate(X.dim1(),3);
  grad_sum_U = 0.0;
  double firstDeriv;

  // Get values from JP that will be needed during the calc

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  Array1D<QMCCorrelationFunctionParameters> * EupNuclear =
    JP.getElectronUpNuclearParameters();

  Array1D<QMCCorrelationFunctionParameters> * EdnNuclear =
    JP.getElectronDownNuclearParameters();

  // Loop over each atom calculating the e-n jastrow function
  double r;
  double * unitVector = new double[3];
  for(int Nuclei=0; Nuclei<Input->Molecule.getNumberAtoms(); Nuclei++)
    {
      // Find the number of the current nucleus in the nuclei list
      int CurrentNucleiNumber = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
        {
          if( Input->Molecule.Atom_Labels(Nuclei) ==
              (*NucleiTypes)(i) )
            {
              CurrentNucleiNumber = i;
            }
        }

      for(int Electron=0; Electron<X.dim1(); Electron++)
        {
          // Find the unit vector between the nucleus and the electron and
          // their distance apart

          r = 0;
          for(int i=0; i<3; i++)
            {
              unitVector[i] = X(Electron,i) - Input->Molecule.Atom_Positions(Nuclei,i);
              r += unitVector[i] * unitVector[i];
            }
          r = sqrt( r );

          // Get the correct correlation function to use and evaluate it

          QMCCorrelationFunction *U_Function = 0;

          if( Electron < Input->WF.getNumberAlphaElectrons() )
            {
              U_Function =
                (*EupNuclear)(CurrentNucleiNumber).getCorrelationFunction();
            }
          else
            {
              U_Function =
                (*EdnNuclear)(CurrentNucleiNumber).getCorrelationFunction();
            }

          U_Function->evaluate(r);
          firstDeriv = U_Function->getFirstDerivativeValue();
          // Update the values being calculated ...

          sum_U +=  U_Function->getFunctionValue();

          laplacian_sum_U += 2.0/r * firstDeriv +
                             U_Function->getSecondDerivativeValue();

          for(int i=0; i<3; i++)
            {
              unitVector[i] /= r;
              grad_sum_U(Electron,i) +=
                firstDeriv * unitVector[i];
            }
        }
    }

  delete [] unitVector;
  unitVector = 0;
}
