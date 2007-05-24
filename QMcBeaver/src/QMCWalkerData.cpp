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
  rp_a.deallocate();
  p3_xxa.deallocate();

  gradPsiRatio.deallocate();
  modifiedGradPsiRatio.deallocate();
  SCF_Grad_PsiRatio.deallocate();
}

void QMCWalkerData::initialize(QMCInput * INPUT, int numDimensions,
			       int numNucForceDim1, int numNucForceDim2,
			       int numAI)
{
  Input = INPUT;
  int numElectrons  = Input->WF.getNumberElectrons();

  gradPsiRatio.allocate(numElectrons,numDimensions);
  modifiedGradPsiRatio.allocate(numElectrons,numDimensions);
  SCF_Grad_PsiRatio.allocate(numElectrons,numDimensions);

  rp_a.allocate(numAI);
  p3_xxa.allocate(numAI);

  modificationRatio    = 0.0;
  psi                  = 0.0;
  isSingular           = false;
    
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
  strm << "       Psi             = " << rhs.psi << endl;
  strm.precision(15);  strm.width(20);
  strm << "       Energy          = " << rhs.localEnergy << endl;
  return strm;
}



