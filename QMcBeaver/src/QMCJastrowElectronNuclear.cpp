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

QMCJastrowElectronNuclear::QMCJastrowElectronNuclear()
{

}

QMCJastrowElectronNuclear::~QMCJastrowElectronNuclear()
{
  grad_sum_U.deallocate();

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).deallocate();

  p_a.deallocate();
  p2_xa.deallocate();
  p3_xxa.deallocate();
}

void QMCJastrowElectronNuclear::initialize(QMCInput * input)
{
  Input = input;

  grad_sum_U.allocate(Input->WF.getNumberElectrons(),3);

  int numNE = Input->JP.getNumberNEParameters();
  p_a.allocate(numNE);
  p2_xa.allocate(numNE);
  p3_xxa.allocate(numNE);

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).allocate(Input->WF.getNumberElectrons(),3);
}

/**
   Find the unit vector and distance between X1 and X2.  The unit vector is in
   the direction of X1-X2.
*/

void QMCJastrowElectronNuclear::calculateDistanceAndUnitVector(
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


double QMCJastrowElectronNuclear::getLaplacianLnJastrow()
{
  return laplacian_sum_U;
}

double QMCJastrowElectronNuclear::get_p3_xxa_ln(int ai)
{
  return p3_xxa(ai);
}

Array2D<double> * QMCJastrowElectronNuclear::getGradientLnJastrow()
{
  return &grad_sum_U;
}

Array2D<double> * QMCJastrowElectronNuclear::get_p2_xa_ln(int ai)
{
  return &p2_xa(ai);
}

double QMCJastrowElectronNuclear::getLnJastrow()
{
  return sum_U;
}

double QMCJastrowElectronNuclear::get_p_a_ln(int ai)
{
  return p_a(ai);
}

void QMCJastrowElectronNuclear::evaluate(QMCJastrowParameters & JP,
					 Array2D<double> & X)
{
  // initialize the results
  sum_U = 0.0;
  laplacian_sum_U = 0.0;
  grad_sum_U = 0.0;
  double firstDeriv;

  p_a              = 0.0;
  p3_xxa           = 0.0;
  for(int ai=0; ai<p2_xa.dim1(); ai++)
    p2_xa(ai)      = 0.0;

  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();

  // Get values from JP that will be needed during the calc

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  Array1D<QMCCorrelationFunctionParameters> * EupNuclear = 0;
  if (nalpha > 0)
    EupNuclear = JP.getElectronUpNuclearParameters();

  Array1D<QMCCorrelationFunctionParameters> * EdnNuclear = 0;
  if (nbeta > 0)
    EdnNuclear = JP.getElectronDownNuclearParameters();

  // Loop over each atom calculating the e-n jastrow function
  double r;
  static double * unitVector = new double[3];
  for(int Nuclei=0; Nuclei<Input->Molecule.getNumberAtoms(); Nuclei++)
    {

      /*
	For labeling the ai parameters, we need to make sure that the order
	we place the ai derivatives in our vectors lines up with how the ai
	are placed in the parameter vector in QMCJastrowParameters.

	In particular, the loops are inverted relative to what is done here.
	The parameters need to have the nuclei indices as the "fast loop".
	Eup comes first, then if the Jastrows are not linked, we have the Edn
	parameters.

	nuc_E??_shift will bring us to the right nucleus
	ai_shift will bring us to the right Jastrow.
	ai indexes the parameter within the Jastrow.
      */
      int nuc_Eup_shift = 0;
      int nuc_Edn_shift = 0;
      int numP;

      // Find the number of the current nucleus in the nuclei list
      int CurrentNucleiNumber = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( Input->Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    CurrentNucleiNumber = i;
	    break;
	  } else {
	    nuc_Eup_shift += (*EupNuclear)(i).getTotalNumberOfParameters();
	    nuc_Edn_shift += (*EdnNuclear)(i).getTotalNumberOfParameters();
	  }

      for(int Electron=0; Electron<X.dim1(); Electron++)
        {
          // Find the unit vector between the nucleus and the electron and
          // their distance apart

          r = 0;
          for(int i=0; i<3; i++)
            {
              unitVector[i] = X(Electron,i) - 
		Input->Molecule.Atom_Positions(Nuclei,i);
              r += unitVector[i] * unitVector[i];
            }
          r = sqrt( r );

          // Get the correct correlation function to use and evaluate it

          QMCCorrelationFunction *U_Function = 0;

	  int ai_shift;
          if( Electron < nalpha )
	    {
	      U_Function = 
		(*EupNuclear)(CurrentNucleiNumber).getCorrelationFunction();

	      ai_shift = nuc_Eup_shift;
	      numP = (*EupNuclear)(CurrentNucleiNumber).getTotalNumberOfParameters();
	    }
          else
	    {
	      U_Function =
		(*EdnNuclear)(CurrentNucleiNumber).getCorrelationFunction();

	      ai_shift = nuc_Edn_shift;
	      numP = (*EdnNuclear)(CurrentNucleiNumber).getTotalNumberOfParameters();
	      if(Input->flags.link_Jastrow_parameters == 0)
		ai_shift += Input->JP.getNumberNEupParameters();
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
              grad_sum_U(Electron,i) += firstDeriv * unitVector[i];
            }

	  if(globalInput.flags.calculate_Derivatives == 1 &&
	     globalInput.flags.optimize_EN_Jastrows == 1)
	    {
	      for(int ai=0; ai<numP; ai++)
		{
		  p_a(ai+ai_shift)    += U_Function->get_p_a(ai);
		  firstDeriv           = U_Function->get_p2_xa(ai);
		  p3_xxa(ai+ai_shift) += 2.0/r * firstDeriv + U_Function->get_p3_xxa(ai);
		  
		  for(int i=0; i<3; i++)
		    (p2_xa(ai+ai_shift))(Electron,i) += firstDeriv * unitVector[i];
		}
	    }
	}
    }
}
