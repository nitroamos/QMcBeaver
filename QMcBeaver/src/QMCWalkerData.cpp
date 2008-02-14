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
#include "QMCWalkerData.h"
#include "QMCPotential_Energy.h"

using namespace std;

/**
   The QMCWalkerData data type is meant to hold all the information
   calculated by QMCFunction that is useful to QMCWalker. This should
   in effect decouple QMCFunction and QMCWalker from each other enabling
   QMCFunction to be treated a little bit differently without significant
   modifications to QMCWalker.
*/
QMCWalkerData::QMCWalkerData()
{
}

QMCWalkerData::~QMCWalkerData()
{
}

void QMCWalkerData::initialize(QMCInput * INPUT, int numDimensions,
			       int numNucForceDim1, int numNucForceDim2)
{
  Input = INPUT;
  int numElectrons = Input->WF.getNumberElectrons();
  int numNuclei    = Input->Molecule.getNumberAtoms();
  int numCI        = Input->WF.getNumberDeterminants();

  //Initialization probably requires all to be updated immediately
  whichE = -1;
  if(globalInput.flags.one_e_per_iter)
    {
      int na = Input->WF.getNumberAlphaElectrons();
      int nb = Input->WF.getNumberBetaElectrons();
      int no = Input->WF.OrbitalCoeffs.dim1();

      D_invA.allocate(numCI);
      D_invB.allocate(numCI);

      for(int ci=0; ci<numCI; ci++)
	{
	  D_invA(ci).allocate(na,na);
	  D_invB(ci).allocate(nb,nb);
	}
	 
      Laplacian_DA.allocate(na,no);
      Laplacian_DB.allocate(nb,no);
      Grad_DA.allocate(3);
      Grad_DB.allocate(3);
      for(int i=0; i<3; i++)
	{
	  Grad_DA(i).allocate(na,no);
	  Grad_DB(i).allocate(nb,no);
	}

      PsiA.allocate(numCI);
      PsiB.allocate(numCI);
      Laplacian_PsiRatioA.allocate(numCI);
      Laplacian_PsiRatioB.allocate(numCI);
      Grad_PsiRatioA.allocate(numCI,numElectrons,3);
      Grad_PsiRatioB.allocate(numCI,numElectrons,3);
    }

  gradPsiRatio.allocate(numElectrons,numDimensions);
  modifiedGradPsiRatio.allocate(numElectrons,numDimensions);
  SCF_Grad_PsiRatio.allocate(numElectrons,numDimensions);

  int numAI = Input->getNumberAIParameters();
  rp_a.allocate(numAI);
  p3_xxa.allocate(numAI);

  modificationRatio    = 0.0;
  psi                  = 0.0;
  SCF_Psi              = 0.0;
  isSingular           = false;

  rij.allocate(numElectrons, numElectrons);
  rij_uvec.allocate(numElectrons, numElectrons, 3);
  riI.allocate(numElectrons, numNuclei);
  riI_uvec.allocate(numElectrons, numNuclei, 3);
  
  Uij.allocate(numElectrons, numElectrons);
  Uij_x.allocate(numElectrons, numElectrons);
  Uij_xx.allocate(numElectrons, numElectrons);

  if(Input->flags.nuclear_derivatives != "none")
    nuclearDerivatives.allocate(numNucForceDim1, numNucForceDim2);
  
  if (Input->flags.calculate_bf_density == 1)
    chiDensity.allocate(Input->WF.getNumberBasisFunctions());
}

double QMCWalkerData::getModifiedLocalEnergy()
{
  /*
    I was just playing around with energy_cutoff_type, it didn't
    seem to help much.
   */

  if(Input->flags.energy_cutoff_type == "none")
    {
      return localEnergy;
    }
  else if(Input->flags.energy_cutoff_type == "umrigar93")
    {
      //we could use this factor against the kinetic contribution only
      return localEnergy * modificationRatio;
    }
  else if(Input->flags.energy_cutoff_type == "bdjm-push88")
    {
      /*
	When dt = 0.005,  cutoff ~= 28
	When dt = 0.001,  cutoff ~= 63
	When dt = 0.0001, cutoff  = 200
	
	This probably only matters when the wavefunction is too poor.
       */
      static const double cutoff = 2.0 / sqrt(Input->flags.dt);

      double el = localEnergy - Input->flags.energy_estimated_original;
      if(fabs(el) > cutoff)
	{
	  double sign = el < 0.0 ? -1.0 : 1.0;
	  return Input->flags.energy_estimated_original + sign * cutoff;
	}
      return localEnergy;
    } else {
      clog << "ERROR: Unknown energy cutoff type " << Input->flags.energy_cutoff_type
	   << "!" << endl;
      exit(1);
    }
  return 0.0;
}

void QMCWalkerData::zero()
{
  localEnergy            = 0.0;
  kineticEnergy          = 0.0;
  potentialEnergy        = 0.0;
  neEnergy               = 0.0;
  eeEnergy               = 0.0;
  lnJ                    = 0.0;
  SCF_Laplacian_PsiRatio = 0.0;
}

void QMCWalkerData::writeConfigs(Array2D<double> & Positions, double weight)
{
  Input->outputer.writeCorrelatedSamplingConfiguration(Positions,
						       SCF_Laplacian_PsiRatio,
						       SCF_Grad_PsiRatio,
						       lnJ,
						       potentialEnergy,
						       weight);
}

ostream& operator<<(ostream & strm, const QMCWalkerData & rhs)
{
  strm.precision(15);  strm.width(20);
  strm << "       KE              = " << rhs.kineticEnergy << endl;
  strm.precision(15);  strm.width(20);
  strm << "       PE              = " << rhs.potentialEnergy << endl;
  strm.precision(15);  strm.width(20);
  strm << "       EE              = " << rhs.eeEnergy << endl;
  strm.precision(15);  strm.width(20);
  strm << "       NE              = " << rhs.neEnergy << endl;
  strm.precision(15);  strm.width(20);
  strm << "       Psi             = " << (double)(rhs.psi) << endl;
  strm.precision(15);  strm.width(20);
  strm << "       Energy          = " << rhs.localEnergy << endl;
  return strm;
}

void QMCWalkerData::updateDistances(Array2D<double> & R)
{
  if(whichE == -1)
    {
      for(int i=0; i<R.dim1(); i++)
	{
	  for(int j=0; j<i; j++)
	    {
	      rij(i,j) = 
		QMCPotential_Energy::rij(R,i,j);
	      for(int xyz=0; xyz<3; xyz++)
		rij_uvec(i,j,xyz) = (R(i,xyz) - R(j, xyz)) / rij(i,j);
	    }

	  for(int j=0; j<Input->Molecule.getNumberAtoms(); j++)
	    {
	      double r = QMCPotential_Energy::rij(R,Input->Molecule.Atom_Positions,i,j);
	      riI(i,j) = r;
	      for(int xyz=0; xyz<3; xyz++)
		riI_uvec(i,j,xyz) = (R(i,xyz) - Input->Molecule.Atom_Positions(j,xyz)) / r;
	    }
	}
    } else {
      for(int j=0; j<R.dim1(); j++)
	{
	  if(j == whichE) continue;
	  double r = QMCPotential_Energy::rij(R,j,whichE);
	  int e1 = max(whichE,j);
	  int e2 = min(whichE,j);

	  rij(e1,e2) = r;
	  for(int xyz=0; xyz<3; xyz++)
	    rij_uvec(e1,e2,xyz) = (R(e1,xyz) - R(e2, xyz)) / r;
	}

      for(int j=0; j<Input->Molecule.getNumberAtoms(); j++)
	{
	  riI(whichE,j) =
	    QMCPotential_Energy::rij(R,Input->Molecule.Atom_Positions,whichE,j);
	  for(int xyz=0; xyz<3; xyz++)
	    riI_uvec(whichE,j,xyz) = (R(whichE,xyz) - Input->Molecule.Atom_Positions(j,xyz)) / riI(whichE,j);
	}
    }
}


