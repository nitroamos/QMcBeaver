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
#include "MathFunctions.h"

QMCJastrowElectronNuclear::QMCJastrowElectronNuclear(){}

QMCJastrowElectronNuclear::~QMCJastrowElectronNuclear()
{
  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).deallocate();

  p_a.deallocate();
  p2_xa.deallocate();
  p3_xxa.deallocate();
}

void QMCJastrowElectronNuclear::initialize(QMCInput * input)
{
  int numNE = globalInput.JP.getNumberNEParameters();
  p_a.allocate(numNE);
  p2_xa.allocate(numNE);
  p3_xxa.allocate(numNE);

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).allocate(globalInput.WF.getNumberElectrons(),3);
}

double QMCJastrowElectronNuclear::get_p3_xxa_ln(int ai)
{
  return p3_xxa(ai);
}

Array2D<double> * QMCJastrowElectronNuclear::get_p2_xa_ln(int ai)
{
  return &p2_xa(ai);
}

double QMCJastrowElectronNuclear::get_p_a_ln(int ai)
{
  return p_a(ai);
}

void QMCJastrowElectronNuclear::evaluate(QMCJastrowParameters & JP,
					 QMCWalkerData * wd,
					 Array2D<double> & X)
{
  // initialize the results
  double temp;
  double UiA, UiA_x, UiA_xx;

  p_a              = 0.0;
  p3_xxa           = 0.0;
  for(int ai=0; ai<p2_xa.dim1(); ai++)
    p2_xa(ai)      = 0.0;

  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta  = globalInput.WF.getNumberElectrons(false);

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
  for(int Nuclei=0; Nuclei<globalInput.Molecule.getNumberAtoms(); Nuclei++)
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
	if( globalInput.Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    CurrentNucleiNumber = i;
	    break;
	  } else {
	    nuc_Eup_shift += (*EupNuclear)(i).getTotalNumberOfParameters();
	    nuc_Edn_shift += (*EdnNuclear)(i).getTotalNumberOfParameters();
	  }

      for(int Electron=0; Electron<X.dim1(); Electron++)
        {
	  if(wd->whichE != Electron && wd->whichE != -1)
	    continue;

          r = wd->riA(Electron,Nuclei);

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
	      if(globalInput.flags.link_Jastrow_parameters == 0)
		ai_shift += globalInput.JP.getNumberNEupParameters();
	    }

          U_Function->evaluate(r);

	  UiA    = U_Function->getFunctionValue();
          UiA_x  = U_Function->getFirstDerivativeValue();
	  UiA_xx = 2.0/r * UiA_x + U_Function->getSecondDerivativeValue();

          // Update the values being calculated ...

	  wd->U -= wd->UiA(Electron, Nuclei);
          wd->U += UiA;
	  wd->UiA(Electron, Nuclei) = UiA;

	  wd->U_xx -= wd->UiA_xx(Electron,Nuclei);
          wd->U_xx += UiA_xx;
	  wd->UiA_xx(Electron,Nuclei) = UiA_xx;

          for(int i=0; i<3; i++)
	    {
	      temp = wd->UiA_x(Electron,Nuclei,i);
	      wd->U_x(Electron,i) -= temp;
	      temp = UiA_x * wd->riA_uvec(Electron,Nuclei,i);
	      wd->UiA_x(Electron,Nuclei,i) = temp;
	      wd->U_x(Electron,i) += temp;
	    }

	  if(globalInput.flags.calculate_Derivatives == 1 &&
	     globalInput.flags.optimize_EN_Jastrows == 1)
	    {
	      for(int ai=0; ai<numP; ai++)
		{
		  p_a(ai+ai_shift)    += U_Function->get_p_a(ai);
		  UiA_x           = U_Function->get_p2_xa(ai);
		  p3_xxa(ai+ai_shift) += 2.0/r * UiA_x + U_Function->get_p3_xxa(ai);
		  
		  for(int i=0; i<3; i++)
		    (p2_xa(ai+ai_shift))(Electron,i) += UiA_x * wd->riA_uvec(Electron,Nuclei,i);
		}
	    }
	}
    }
}

double QMCJastrowElectronNuclear::jastrowOnGrid(QMCJastrowParameters & JP,
						int Electron,
						Array2D<double> & R,
						Array2D<double> & grid,
						Array1D<double> & integrand)
{
  // initialize the results
  double denom = 0.0;
  double r;
  Array1D<double> sumU(integrand.dim1());
  sumU = 0.0;

  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta  = globalInput.WF.getNumberElectrons(false);

  // Get values from JP that will be needed during the calc

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  Array1D<QMCCorrelationFunctionParameters> * EupNuclear = 0;
  if (nalpha > 0)
    EupNuclear = JP.getElectronUpNuclearParameters();
  
  Array1D<QMCCorrelationFunctionParameters> * EdnNuclear = 0;
  if (nbeta > 0)
    EdnNuclear = JP.getElectronDownNuclearParameters();

  // Loop over each atom calculating the e-n jastrow function
  for(int Nuclei=0; Nuclei<globalInput.Molecule.getNumberAtoms(); Nuclei++)
    {
      // Find the number of the current nucleus in the nuclei list
      int CurrentNucleiNumber = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( globalInput.Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    CurrentNucleiNumber = i;
	    break;
	  }      

      // Get the correct correlation function to use and evaluate it
      QMCCorrelationFunction *U_Function = 0;
      
      if( Electron < nalpha )
	{
	  U_Function = 
	    (*EupNuclear)(CurrentNucleiNumber).getCorrelationFunction();
	}
      else
	{
	  U_Function =
	    (*EdnNuclear)(CurrentNucleiNumber).getCorrelationFunction();
	}
      
      r = MathFunctions::rij(R,globalInput.Molecule.Atom_Positions,
			     Electron,Nuclei);
      denom += U_Function->getFunctionValue(r);
      
      for(int gr=0; gr<grid.dim1(); gr++){
	r = MathFunctions::rij(grid,globalInput.Molecule.Atom_Positions,
			       gr,Nuclei);
	sumU(gr) += U_Function->getFunctionValue(r);
      }
    }

  for(int gr=0; gr<grid.dim1(); gr++)
    integrand(gr) *= exp(sumU(gr));
  return exp(denom);
}

