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

void QMCWalkerData::partialCopy(const QMCWalkerData & rhs)
{
  localEnergy            = rhs.localEnergy;
  kineticEnergy          = rhs.kineticEnergy;
  potentialEnergy        = rhs.potentialEnergy;
  neEnergy               = rhs.neEnergy;
  eeEnergy               = rhs.eeEnergy;

  D                      = rhs.D;
  D_x                    = rhs.D_x;
  D_xx                   = rhs.D_xx;

  psi                    = rhs.psi;
  singular               = rhs.singular;

  gradPsiRatio           = rhs.gradPsiRatio;
  modifiedGradPsiRatio   = rhs.modifiedGradPsiRatio;
  riA = rhs.riA;
  riA_uvec = rhs.riA_uvec;
}

void QMCWalkerData::initialize(int numDimensions, int numNucForceDim1, int numNucForceDim2)
{
  int numElectrons = globalInput.WF.getNumberElectrons();
  int numNuclei    = globalInput.Molecule.getNumberAtoms();
  int numCI        = globalInput.WF.getNumberDeterminants();

  //Initialization probably requires all to be updated immediately
  whichE = -1;
  if(globalInput.flags.one_e_per_iter)
    {
      int na = globalInput.WF.getNumberAlphaElectrons();
      int nb = globalInput.WF.getNumberBetaElectrons();
      int no = globalInput.WF.OrbitalCoeffs.dim1();

      Dc_invA.allocate(numCI);
      Dc_invB.allocate(numCI);

      for(int ci=0; ci<numCI; ci++)
	{
	  Dc_invA(ci).allocate(na,na);
	  Dc_invB(ci).allocate(nb,nb);
	}
	 
      D_xxA.allocate(na,no);
      D_xxB.allocate(nb,no);
      D_xA.allocate(3);
      D_xB.allocate(3);
      for(int i=0; i<3; i++)
	{
	  D_xA(i).allocate(na,no);
	  D_xB(i).allocate(nb,no);
	}

      DcA.allocate(numCI);
      DcB.allocate(numCI);
      rDc_xxA.allocate(numCI);
      rDc_xxB.allocate(numCI);
      rDc_xA.allocate(numCI,numElectrons,3);
      rDc_xB.allocate(numCI,numElectrons,3);
    }

  gradPsiRatio.allocate(numElectrons,numDimensions);
  modifiedGradPsiRatio.allocate(numElectrons,numDimensions);
  D_x.allocate(numElectrons,numDimensions);

  int numAI = globalInput.getNumberAIParameters();
  rp_a.allocate(numAI);
  p3_xxa.allocate(numAI);

  modificationRatio    = 0.0;
  psi                  = 0.0;
  D                    = 0.0;
  singular             = false;
  
  rij.allocate(numElectrons, numElectrons);
  rij_uvec.allocate(numElectrons, numElectrons, 3);
  riA.allocate(numElectrons, numNuclei);
  riA_uvec.allocate(numElectrons, numNuclei, 3);
  
  Uij.allocate(numElectrons, numElectrons);
  Uij_x.allocate(numElectrons, numElectrons, 3);
  Uij_xx.allocate(numElectrons, numElectrons);
  Uij    = 0.0;
  Uij_x  = 0.0;
  Uij_xx = 0.0;

  UiA.allocate(numElectrons, numNuclei);
  UiA_x.allocate(numElectrons, numNuclei, 3);
  UiA_xx.allocate(numElectrons, numNuclei);
  UiA    = 0.0;
  UiA_x  = 0.0;
  UiA_xx = 0.0;
  
  if(globalInput.flags.use_three_body_jastrow == 1)
    {
      UijA.allocate(numElectrons, numElectrons);
      UijA_xx.allocate(numElectrons, numElectrons);
      UijA    = 0.0;
      UijA_xx = 0.0;

      UijA_x1.allocate(3);
      UijA_x2.allocate(3);      
      for(int xyz=0; xyz<UijA_x1.dim1(); xyz++)
	{
	  UijA_x1(xyz).allocate(numElectrons,numElectrons);
	  UijA_x1(xyz) = 0.0;
	  UijA_x2(xyz).allocate(numElectrons,numElectrons);
	  UijA_x2(xyz) = 0.0;
	}
    }

  U_x.allocate(numElectrons,3);
  U    = 0.0;
  U_x  = 0.0;
  U_xx = 0.0;

  if(globalInput.flags.nuclear_derivatives != "none")
    nuclearDerivatives.allocate(numNucForceDim1, numNucForceDim2);
  
  if (globalInput.flags.calculate_bf_density == 1)
    chiDensity.allocate(globalInput.WF.getNumberBasisFunctions());
}

double QMCWalkerData::getModifiedLocalEnergy()
{
  /*
    I was just playing around with energy_cutoff_type, it didn't
    seem to help much.
   */

  if(globalInput.flags.energy_cutoff_type == "none")
    {
      return localEnergy;
    }
  else if(globalInput.flags.energy_cutoff_type == "umrigar93")
    {
      //we could use this factor against the kinetic contribution only
      return localEnergy * modificationRatio;
    }
  else if(globalInput.flags.energy_cutoff_type == "bdjm-push88")
    {
      /*
	When dt = 0.005,  cutoff ~= 28
	When dt = 0.001,  cutoff ~= 63
	When dt = 0.0001, cutoff  = 200
	
	This probably only matters when the wavefunction is too poor.
       */
      static const double cutoff = 2.0 / sqrt(globalInput.flags.dt);

      double el = localEnergy - globalInput.flags.energy_estimated_original;
      if(fabs(el) > cutoff)
	{
	  double sign = el < 0.0 ? -1.0 : 1.0;
	  return globalInput.flags.energy_estimated_original + sign * cutoff;
	}
      return localEnergy;
    } else {
      clog << "ERROR: Unknown energy cutoff type " << globalInput.flags.energy_cutoff_type
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
  D_xx                   = 0.0;
}

void QMCWalkerData::writeConfigs(Array2D<double> & Positions, double weight)
{
  globalInput.outputer.writeCorrelatedSamplingConfiguration(Positions,
							    D_xx,
							    D_x,
							    U,
							    potentialEnergy,
							    weight);
}

bool QMCWalkerData::isSingular()
{
  if(singular) return true;
  if( IeeeMath::isNaN(kineticEnergy) )
    singular = true;
  if( IeeeMath::isNaN(potentialEnergy) )
    singular = true;
  if( IeeeMath::isNaN(U_xx) )
    singular = true;
  if( IeeeMath::isNaN(D_xx) )
    singular = true;
  return singular;
}

ostream& operator<<(ostream & strm, const QMCWalkerData & rhs)
{
  int w = 25;
  int p = 15;

  strm.width(10);
  strm << "  Psi = "
       << setw(w) << setprecision(p)
       << scientific << (double)(rhs.psi);
  strm.width(10);
  strm << "    E = "
       << setw(w) << setprecision(p)
       << fixed << rhs.localEnergy;
  strm.width(10);
  if(rhs.singular)
    strm << " * ";
  strm << endl;


  strm.width(10);
  strm << "   U = "
       << setw(w) << setprecision(p)
       << scientific << (double)(QMCDouble(1.0,1.0,0.0,rhs.U));
  strm.width(10);
  strm << " U_xx = "
       << setw(w) << setprecision(p)
       << fixed << rhs.U_xx;
  strm << endl;

  strm.width(10);
  strm << "   D = "
       << setw(w) << setprecision(p)
       << scientific << (double)(rhs.D);
  strm.width(10);
  strm << " D_xx = "
       << setw(w) << setprecision(p)
       << fixed << rhs.D_xx;
  strm << endl;

  strm.width(10);
  strm << "   KE = "
       << setw(w) << setprecision(p)
       << fixed << rhs.kineticEnergy;
  strm.width(10);
  strm << "   PE = "
       << setw(w) << setprecision(p)
       << fixed << rhs.potentialEnergy;
  strm << endl;

  strm.width(10);
  strm << "   EE = "
       << setw(w) << setprecision(p)
       << fixed << rhs.eeEnergy;
  strm.width(10);
  strm << "   NE = "
       << setw(w) << setprecision(p)
       << fixed << rhs.neEnergy;
  strm << endl;

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

	  for(int j=0; j<globalInput.Molecule.getNumberAtoms(); j++)
	    {
	      double r = QMCPotential_Energy::rij(R,globalInput.Molecule.Atom_Positions,i,j);
	      riA(i,j) = r;
	      for(int xyz=0; xyz<3; xyz++)
		riA_uvec(i,j,xyz) = (R(i,xyz) - globalInput.Molecule.Atom_Positions(j,xyz)) / r;
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

      for(int j=0; j<globalInput.Molecule.getNumberAtoms(); j++)
	{
	  riA(whichE,j) =
	    QMCPotential_Energy::rij(R,globalInput.Molecule.Atom_Positions,whichE,j);
	  for(int xyz=0; xyz<3; xyz++)
	    riA_uvec(whichE,j,xyz) = (R(whichE,xyz) - globalInput.Molecule.Atom_Positions(j,xyz)) / riA(whichE,j);
	}
    }
}


