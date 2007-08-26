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

#include "QMCThreeBodyJastrow.h"

QMCThreeBodyJastrow::QMCThreeBodyJastrow()
{
}

QMCThreeBodyJastrow::~QMCThreeBodyJastrow()
{
  grad_sum_U.deallocate();

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).deallocate();

  p_a.deallocate();
  p2_xa.deallocate();
  p3_xxa.deallocate();
}

void QMCThreeBodyJastrow::initialize(QMCInput * input)
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
   Finds the relative positions of two electrons with respect to a nucleus.
*/

void QMCThreeBodyJastrow::calculateDistances(Array2D<double> &X1, 
					     int x1particle, int x2particle,
					     Array2D<double> &X3,
					     int x3particle,
					     Array1D<double> &position1,
					     double &r1,
					     Array1D<double> &position2,
					     double &r2)

{
  r1 = 0.0;
  r2 = 0.0;

  for (int i=0; i<3; i++)
    {
      position1(i) = X1(x1particle,i) - X3(x3particle,i);
      r1 += position1(i)*position1(i);

      position2(i) = X1(x2particle,i) - X3(x3particle,i);
      r2 += position2(i)*position2(i);
    }

  r1 = sqrt(r1);
  r2 = sqrt(r2);
}


double QMCThreeBodyJastrow::getLaplacianLnJastrow()
{
  return laplacian_sum_U;
}

double QMCThreeBodyJastrow::get_p3_xxa_ln(int ai)
{
  return p3_xxa(ai);
}

Array2D<double> * QMCThreeBodyJastrow::getGradientLnJastrow()
{
  return &grad_sum_U;
}

Array2D<double> * QMCThreeBodyJastrow::get_p2_xa_ln(int ai)
{
  return &p2_xa(ai);
}

double QMCThreeBodyJastrow::getLnJastrow()
{
  return sum_U;
}

double QMCThreeBodyJastrow::get_p_a_ln(int ai)
{
  return p_a(ai);
}

void QMCThreeBodyJastrow::evaluate(QMCJastrowParameters & JP, 
				   QMCWalkerData * wd, Array2D<double> & X)
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

  Array1D<QMCThreeBodyCorrelationFunctionParameters> * EupEdnNuclear = 0;
  if (nalpha > 0 && nbeta > 0)
    EupEdnNuclear = JP.getElectronUpElectronDownNuclearParameters();

  Array1D<QMCThreeBodyCorrelationFunctionParameters> * EupEupNuclear = 0;
  if (nalpha > 1)
    EupEupNuclear = JP.getElectronUpElectronUpNuclearParameters();

  Array1D<QMCThreeBodyCorrelationFunctionParameters> * EdnEdnNuclear = 0;
  if (nbeta > 1)
    EdnEdnNuclear = JP.getElectronDownElectronDownNuclearParameters();

  // Loop over each atom calculating the three body jastrow function

  double dist1 = 0.0;
  double dist2 = 0.0;
  double cutoff = 0.0;

  Array1D<double> xyz1(3);
  xyz1 = 0.0;

  Array1D<double> xyz2(3);
  xyz2 = 0.0;

  Array1D<double> *grad1;
  Array1D<double> *grad2;

  QMCThreeBodyCorrelationFunction *U_function = 0;
  int ai_shift = 0;

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
      int nucEupEdnShift = 0;
      int nucEupEupShift = 0;
      int nucEdnEdnShift = 0;
      int numP;

      // Find the number of the current nucleus in the nuclei list
      int CurrentNucleiNumber = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( Input->Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    CurrentNucleiNumber = i;
	    break;
	  } 
	else 
	  {
	    nucEupEdnShift += (*EupEdnNuclear)(i).getNumberOfFreeParameters();
	    nucEupEupShift += (*EupEupNuclear)(i).getNumberOfFreeParameters();
	    nucEdnEdnShift += (*EdnEdnNuclear)(i).getNumberOfFreeParameters();
	  }

      // First we do all pairs of two alphas with this nucleus

      if (nalpha > 1)
	{
	  U_function =
       (*EupEupNuclear)(CurrentNucleiNumber).getThreeBodyCorrelationFunction();

	  ai_shift = nucEupEupShift;
	  numP = 
	    (*EupEupNuclear)(CurrentNucleiNumber).getNumberOfFreeParameters();

	  cutoff = (*EupEupNuclear)(CurrentNucleiNumber).getCutoffDist();

	  for(int Electron1=0; Electron1<nalpha-1; Electron1++)
	    for (int Electron2=Electron1+1; Electron2<nalpha; Electron2++)
	      {
		// Find the distances between the two electrons and the nucleus

		xyz1 = 0.0;
		dist1 = 0.0;

		xyz2 = 0.0;
		dist2 = 0.0;
	  
		calculateDistances(X, Electron1, Electron2,
				   Input->Molecule.Atom_Positions, Nuclei,
				   xyz1,dist1, xyz2, dist2);

		if (dist1 < cutoff && dist2 < cutoff)
		  {
		    U_function->evaluate(xyz1,dist1,xyz2,dist2);

		    sum_U += U_function->getFunctionValue();

		    grad1 = U_function->getElectron1Gradient();
		    grad2 = U_function->getElectron2Gradient();
	      
		    for (int i=0; i<3; i++)
		      {
			grad_sum_U(Electron1,i) += (*grad1)(i);
			grad_sum_U(Electron2,i) += (*grad2)(i);
		      }

		    laplacian_sum_U += U_function->getLaplacianValue();
		  }

	  /*
	    This part has to do with the parameter derivatives, and I haven't
	    figured out how this works yet.
	  if(globalInput.flags.calculate_Derivatives == 1 &&
	     globalInput.flags.optimize_NEE_Jastrows == 1)
	    for(int ai=0; ai<numP; ai++)
	      {
		p_a(ai+ai_shift)    += U_Function->get_p_a(ai);
		firstDeriv           = U_Function->get_p2_xa(ai);
		p3_xxa(ai+ai_shift) += 2.0/dist1 * firstDeriv + U_Function->get_p3_xxa(ai);
		  
		for(int i=0; i<3; i++)
		  (p2_xa(ai+ai_shift))(Electron1,i) += firstDeriv; 
	      }
	  */
	      }
	}
      // Now we do all opposite spin pairs with this nucleus

	 
      if (nalpha > 0 && nbeta > 0)
	{
	  U_function = 
       (*EupEdnNuclear)(CurrentNucleiNumber).getThreeBodyCorrelationFunction();

	  ai_shift = nucEupEdnShift;
	  numP = 
	    (*EupEdnNuclear)(CurrentNucleiNumber).getNumberOfFreeParameters();

	  cutoff = (*EupEdnNuclear)(CurrentNucleiNumber).getCutoffDist();

	  for(int Electron1=0; Electron1<nalpha; Electron1++)
	    for (int Electron2=nalpha; Electron2<nalpha+nbeta; Electron2++)
	      {
		// Find the distances between the two electrons and the nucleus

		xyz1 = 0.0;
		dist1 = 0.0;
		
		xyz2 = 0.0;
		dist2 = 0.0;
	  
		calculateDistances(X, Electron1, Electron2,
				   Input->Molecule.Atom_Positions, Nuclei,
				   xyz1,dist1, xyz2, dist2);

		if (dist1 < cutoff && dist2 < cutoff)
		  {
		    U_function->evaluate(xyz1,dist1,xyz2,dist2);

		    sum_U += U_function->getFunctionValue();

		    grad1 = U_function->getElectron1Gradient();
		    grad2 = U_function->getElectron2Gradient();
	      
		    for (int i=0; i<3; i++)
		      {
			grad_sum_U(Electron1,i) += (*grad1)(i);
			grad_sum_U(Electron2,i) += (*grad2)(i);
		      }

		    laplacian_sum_U += U_function->getLaplacianValue();
		  }

	  /*
	    This part has to do with the parameter derivatives.

	  if(globalInput.flags.calculate_Derivatives == 1 &&
	     globalInput.flags.optimize_NEE_Jastrows == 1)
	    for(int ai=0; ai<numP; ai++)
	      {
		p_a(ai+ai_shift)    += U_Function->get_p_a(ai);
		firstDeriv           = U_Function->get_p2_xa(ai);
		p3_xxa(ai+ai_shift) += 2.0/r * firstDeriv + U_Function->get_p3_xxa(ai);
		  
		for(int i=0; i<3; i++)
		  (p2_xa(ai+ai_shift))(Electron,i) += firstDeriv * unitVector[i];
	      }
	  */
	      }
	}

      // Finally we do all beta spin pairs with this nucleus
      
      if (nbeta > 1)
	{
	  U_function =
       (*EdnEdnNuclear)(CurrentNucleiNumber).getThreeBodyCorrelationFunction();

	  ai_shift = nucEdnEdnShift;
	  numP = 
	    (*EdnEdnNuclear)(CurrentNucleiNumber).getNumberOfFreeParameters();

	  cutoff = (*EdnEdnNuclear)(CurrentNucleiNumber).getCutoffDist();

	  for(int Electron1=nalpha; Electron1<nalpha+nbeta-1; Electron1++)
	    for (int Electron2=Electron1+1; Electron2<nalpha+nbeta; Electron2++)
	      {
		// Find the distances between the two electrons and the nucleus

		xyz1 = 0.0;
		dist1 = 0.0;

		xyz2 = 0.0;
		dist2 = 0.0;
	  
		calculateDistances(X, Electron1, Electron2,
				   Input->Molecule.Atom_Positions, Nuclei,
				   xyz1,dist1, xyz2, dist2);

		if (dist1 < cutoff && dist2 < cutoff)
		  {
		    U_function->evaluate(xyz1,dist1,xyz2,dist2);
		    
		    sum_U += U_function->getFunctionValue();

		    grad1 = U_function->getElectron1Gradient();
		    grad2 = U_function->getElectron2Gradient();
	      
		    for (int i=0; i<3; i++)
		      {
			grad_sum_U(Electron1,i) += (*grad1)(i);
			grad_sum_U(Electron2,i) += (*grad2)(i);
		      }

		    laplacian_sum_U += U_function->getLaplacianValue();
		  }
		
	  /*
	    This part has to do with the parameter derivatives.
	  
	  if(globalInput.flags.calculate_Derivatives == 1 &&
	     globalInput.flags.optimize_EEN_Jastrows == 1)
	    for(int ai=0; ai<numP; ai++)
	      {
		p_a(ai+ai_shift)    += U_Function->get_p_a(ai);
		firstDeriv           = U_Function->get_p2_xa(ai);
		p3_xxa(ai+ai_shift) += 2.0/r * firstDeriv + U_Function->get_p3_xxa(ai);
		  
		for(int i=0; i<3; i++)
		  (p2_xa(ai+ai_shift))(Electron,i) += firstDeriv * unitVector[i];
	      }
	  */
	      }
	}
    }
}







