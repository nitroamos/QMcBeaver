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
}

void QMCThreeBodyJastrow::initialize(QMCInput * input)
{
  Input = input;

  grad_sum_U.allocate(Input->WF.getNumberElectrons(),3);
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
				   QMCWalkerData * wData, Array2D<double> & X)
{
  wd = wData;
  // initialize the results
  sum_U = 0.0;
  laplacian_sum_U = 0.0;
  grad_sum_U = 0.0;

  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();

  // Get values from JP that will be needed during the calc

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  if (nalpha > 0 && nbeta > 0)
    {
      EupEdnNuclear = JP.getElectronUpElectronDownNuclearParameters();
      for(int nuc=0; nuc<EupEdnNuclear->dim1(); nuc++)
	(*EupEdnNuclear)(nuc).zeroOutDerivatives();
    }
  else
    EupEdnNuclear = 0;    

  if (nalpha > 1)
    {
      EupEupNuclear = JP.getElectronUpElectronUpNuclearParameters();
      for(int nuc=0; nuc<EupEupNuclear->dim1(); nuc++)
	(*EupEupNuclear)(nuc).zeroOutDerivatives();
    }
  else
    EupEupNuclear = 0;    
  
  if (nbeta > 1)
    {
      EdnEdnNuclear = JP.getElectronDownElectronDownNuclearParameters();
      for(int nuc=0; nuc<EdnEdnNuclear->dim1(); nuc++)
	(*EdnEdnNuclear)(nuc).zeroOutDerivatives();
    }
  else
    EdnEdnNuclear = 0;

  for(int Nuclei=0; Nuclei<Input->Molecule.getNumberAtoms(); Nuclei++)
    {
      // Find the number of the current nucleus in the nuclei list
      int NucleiType = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( Input->Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    NucleiType = i;
	    break;
	  } 

      // Now we do all opposite spin pairs with this nucleus	 
      if (nalpha > 0 && nbeta > 0)
	for(int Electron1=0; Electron1<nalpha; Electron1++)
	  for (int Electron2=nalpha; Electron2<nalpha+nbeta; Electron2++)
	    collectForPair(Electron2, Electron1, Nuclei,
			   & (*EupEdnNuclear)(NucleiType));

      // First we do all pairs of two alphas with this nucleus
      if (nalpha > 1)
	for(int Electron1=0; Electron1<nalpha; Electron1++)
	  for(int Electron2=0; Electron2<Electron1; Electron2++)	      
	    collectForPair(Electron1, Electron2, Nuclei,
			   & (*EupEupNuclear)(NucleiType));
      
      // Finally we do all beta spin pairs with this nucleus
      if (nbeta > 1)
	for(int Electron1=nalpha; Electron1<X.dim1(); Electron1++)
	  for(int Electron2=nalpha; Electron2<Electron1; Electron2++)
	    collectForPair(Electron1, Electron2, Nuclei,
			   & (*EdnEdnNuclear)(NucleiType));
    }

  packageDerivatives();
  wd = 0;
}

inline void QMCThreeBodyJastrow::collectForPair(int Electron1, 
						int Electron2,
						int Nuclei,
						QMCThreeBodyCorrelationFunctionParameters * paramset)
{
  QMCThreeBodyCorrelationFunction *U_Function = paramset->getThreeBodyCorrelationFunction();

  // Find the distances between the two electrons and the nucleus
  double dist1 = wd->riI(Electron1,Nuclei);
  double dist2 = wd->riI(Electron2,Nuclei);
  double r12   = wd->rij(Electron1,Electron2);

  Array1D<double> xyz1(3);
  Array1D<double> xyz2(3);
  Array1D<double> xyz12(3);

  for(int i=0; i<3; i++)
    {
      xyz1(i)  = wd->riI_uvec(Electron1,Nuclei,i);
      xyz2(i)  = wd->riI_uvec(Electron2,Nuclei,i);
      xyz12(i) = wd->rij_uvec(Electron1, Electron2,i);
    }

  U_Function->evaluate(xyz1,  dist1,
		       xyz2,  dist2,
		       xyz12, r12);
  
  sum_U += U_Function->getFunctionValue();

  Array1D<double> * grad1 = U_Function->getElectron1Gradient();
  Array1D<double> * grad2 = U_Function->getElectron2Gradient();
  
  for (int i=0; i<3; i++)
    {
      grad_sum_U(Electron1,i) += (*grad1)(i);
      grad_sum_U(Electron2,i) += (*grad2)(i);
    }
  
  laplacian_sum_U += U_Function->getLaplacianValue();

  if(globalInput.flags.calculate_Derivatives == 1 &&
     globalInput.flags.optimize_NEE_Jastrows == 1)
    {
      for(int ai=0; ai<paramset->getNumberOfTotalParameters(); ai++)
	{
	  paramset->pt_a(ai)    += U_Function->get_p_a(ai);
	  paramset->pt3_xxa(ai) += U_Function->get_p3_xxa(ai);
	  
	  for(int i=0; i<3; i++)
	    {
	      (paramset->pt2_xa(ai))(Electron1,i) += U_Function->get_p2_xa(true,i,ai); 
	      (paramset->pt2_xa(ai))(Electron2,i) += U_Function->get_p2_xa(false,i,ai); 
	    }
	}
    }
}

void QMCThreeBodyJastrow::packageDerivatives()
{
  int numNEE =
    globalInput.JP.getNumberNEupEdnParameters() +
    globalInput.JP.getNumberNEupEupParameters() + 
    globalInput.JP.getNumberNEdnEdnParameters();  
  
  p_a.allocate(numNEE);
  p2_xa.allocate(numNEE);
  p3_xxa.allocate(numNEE);
  for(int i=0; i<p2_xa.dim1(); i++)
    {
      p2_xa(i).allocate(globalInput.WF.getNumberElectrons(),3);
      p2_xa(i) = 0.0;
    }
  p_a = 0.0;
  p3_xxa = 0.0;
  
  int ai = 0;
  int upupstart;

  if(EupEdnNuclear != 0)
    for(int nuc=0; nuc<EupEdnNuclear->dim1(); nuc++)
      {
	(*EupEdnNuclear)(nuc).totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	ai += (*EupEdnNuclear)(nuc).getNumberOfFreeParameters();
      }
  upupstart = ai;

  if(EupEupNuclear != 0)
    for(int nuc=0; nuc<EupEupNuclear->dim1(); nuc++)
      {
	(*EupEupNuclear)(nuc).totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	ai += (*EupEupNuclear)(nuc).getNumberOfFreeParameters();
      }  

  /*
    If upup and dndn are the same, then their derivatives
    are added to the same accumulator.
  */
  if(globalInput.flags.link_Jastrow_parameters == 1)
    ai = upupstart;

  if(EdnEdnNuclear != 0)
    for(int nuc=0; nuc<EdnEdnNuclear->dim1(); nuc++)
      {
	(*EdnEdnNuclear)(nuc).totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	ai += (*EdnEdnNuclear)(nuc).getNumberOfFreeParameters();
      }

  if(ai != numNEE)
    {
      clog << "Error: parameter counting went bad in packageDerivatives" << endl;
      clog << "     ai = " << ai << endl;
      clog << " numNEE = " << numNEE << endl;
    }
}
