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

QMCJastrowElectronElectron::QMCJastrowElectronElectron()
{

}

QMCJastrowElectronElectron::~QMCJastrowElectronElectron()
{
  grad_sum_U.deallocate();

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).deallocate();

  p_a.deallocate();
  p2_xa.deallocate();
  p3_xxa.deallocate();
}

void QMCJastrowElectronElectron::initialize(QMCInput * input)
{
  Input = input;

  grad_sum_U.allocate(Input->WF.getNumberElectrons(),3);

  int numEE = Input->JP.getNumberEEParameters();
  p_a.allocate(numEE);
  p2_xa.allocate(numEE);
  p3_xxa.allocate(numEE);

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).allocate(Input->WF.getNumberElectrons(),3);
}

/**
   Find the unit vector and distance between X1 and X2.  The unit vector is in
   the direction of X1-X2.
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

double QMCJastrowElectronElectron::get_p3_xxa_ln(int ai)
{
  return p3_xxa(ai);
}

Array2D<double> * QMCJastrowElectronElectron::getGradientLnJastrow()
{
  return &grad_sum_U;
}

Array2D<double> * QMCJastrowElectronElectron::get_p2_xa_ln(int ai)
{
  return &p2_xa(ai);
}

double QMCJastrowElectronElectron::getLnJastrow()
{
  return sum_U;
}

double QMCJastrowElectronElectron::get_p_a_ln(int ai)
{
  return p_a(ai);
}

void QMCJastrowElectronElectron::evaluate(QMCJastrowParameters & JP,
    Array2D<double> & X)
{
  // initialize the results
  sum_U            = 0.0;
  laplacian_sum_U  = 0.0;
  grad_sum_U       = 0.0;

  p_a              = 0.0;
  p3_xxa           = 0.0;
  for(int ai=0; ai<p2_xa.dim1(); ai++)
    p2_xa(ai)      = 0.0;

  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();

  // Get values from JP that will be needed during the calc

  QMCCorrelationFunctionParameters * EupEdn = 0;

  if(nalpha > 0 && nbeta > 0)
      EupEdn = JP.getElectronUpElectronDownParameters();

  QMCCorrelationFunctionParameters * EupEup = 0;

  if(nalpha > 1)
      EupEup = JP.getElectronUpElectronUpParameters();

  QMCCorrelationFunctionParameters * EdnEdn = 0;

  if(nbeta > 1 )
      EdnEdn = JP.getElectronDownElectronDownParameters();

  // Loop over each electron calculating the e-e jastrow function
  //Array1D<double> UnitVector;

  // Get the correct correlation function to use and evaluate it
  //I separated the collectForPair so that the inner loop didn't have
  //a giant if statement.

  QMCCorrelationFunction *U_Function = 0;

  int index = 0;
  int numP  = 0;

  if(EupEdn != 0)
    {
      U_Function = EupEdn->getCorrelationFunction();
      numP = EupEdn->getTotalNumberOfParameters();
    }
  for(int Electron1=0; Electron1<nalpha; Electron1++)
    for(int Electron2=nalpha; Electron2<X.dim1(); Electron2++)
	collectForPair(Electron1,Electron2,U_Function,X,index,numP);

  index += numP;

  if(EupEup != 0)
    {
      numP  = EupEup->getTotalNumberOfParameters();
      U_Function = EupEup->getCorrelationFunction();
    } else {
      numP = 0;
    }

  for(int Electron1=0; Electron1<nalpha; Electron1++)
    for(int Electron2=0; Electron2<Electron1; Electron2++)
	collectForPair(Electron1,Electron2,U_Function,X,index,numP);

  if(Input->flags.link_Jastrow_parameters == 0)
    index += numP;


  if(EdnEdn != 0)
    {
      numP  = EdnEdn->getTotalNumberOfParameters();
      U_Function = EdnEdn->getCorrelationFunction();
    } else {
      numP  = 0;
    }
  for(int Electron1=nalpha; Electron1<X.dim1(); Electron1++)
    for(int Electron2=nalpha; Electron2<Electron1; Electron2++)
      collectForPair(Electron1,Electron2,U_Function,X,index,numP);
}

inline void QMCJastrowElectronElectron::collectForPair(int Electron1, 
						       int Electron2,
						       QMCCorrelationFunction *U_Function,
						       Array2D<double> & X,
						       int index, int numP)
{
  // Find the unit vector between electron1 and electron2 and their
  // distance apart
  static double * unitVector = new double[3];
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
  laplacian_sum_U += 2.0*(2.0/r * firstDeriv +
                          U_Function->getSecondDerivativeValue());

  for(int i=0; i<3; i++)
    {
      unitVector[i] /= r;
      temp = firstDeriv * unitVector[i];

      grad_sum_U(Electron1,i) += temp;
      grad_sum_U(Electron2,i) -= temp;
    }

  /*
    These are for calculating the derivative of the energy
    with respect to parameter ai.
  */
  if(Input->flags.calculate_Derivatives == 1)
    {
      for(int ai=0; ai<numP; ai++)
	{
	  
	  p_a(ai+index)    += U_Function->get_p_a(ai);
	  firstDeriv        = U_Function->get_p2_xa(ai);
	  p3_xxa(ai+index) += 2.0*(2.0/r * firstDeriv +
				   U_Function->get_p3_xxa(ai));
	  
	  for(int i=0; i<3; i++)
	    {
	      temp = firstDeriv * unitVector[i];
	      
	      (p2_xa(ai+index))(Electron1,i) += temp;
	      (p2_xa(ai+index))(Electron2,i) -= temp;
	    }
	}
    }
}
