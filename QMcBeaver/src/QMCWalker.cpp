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

#include "QMCWalker.h"

const double QMCWalker::maxFWAsymp = 1.0;
long int QMCWalker::nextID = 0;

QMCWalker::QMCWalker()
{
  TrialWalker    = 0;
  OriginalWalker = 0;

  weight   = 1.0;
  age      = 0;
  dR2      = 0.0;
  ageMoved = -1;
}

void QMCWalker::branchID()
{
  for(int i=numAncestors-1; i > 0; i--)
    genealogy[i] = genealogy[i-1];
  genealogy[0] = ++nextID;
}

void QMCWalker::newID()
{
  genealogy[0] = ++nextID;
  for(int i=1; i < numAncestors; i++)
    genealogy[i] = -1;
}

QMCWalker::QMCWalker( const QMCWalker & rhs )
{
  TrialWalker    = 0;
  OriginalWalker = 0;

  *this = rhs;
}

QMCWalker::~QMCWalker()
{
  TrialWalker = 0;

  delete OriginalWalker;
  OriginalWalker = 0;  
  
  Input = 0;
  
  R.deallocate();
  
  numFWSteps.clear();

  isCollectingFWResults.deallocate();  
  fwNormalization.deallocate();
  fwR12.deallocate();
  fwR2.deallocate();
  fwiR12.deallocate();
  fwiR.deallocate();
  fwEnergy.deallocate();
  fwKineticEnergy.deallocate();
  fwKineticEnergy_grad.deallocate();
  fwPotentialEnergy.deallocate();

  for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
    for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
      fwNuclearForces(d1,d2).deallocate();
  fwNuclearForces.deallocate();  
}

void QMCWalker::operator=( const QMCWalker & rhs )
{
  // In the equality operator DON'T copy the pointers to the child walkers.
  // If you do, then there is the possibility of big memory leak problems.
  // This is not required because these are just temporary variables used in
  // propagating the electrons.
  
  weight   = rhs.weight;
  age      = rhs.age;
  ageMoved = rhs.ageMoved;

  for(int i=0; i<numAncestors; i++)
    genealogy[i] = rhs.genealogy[i];

  Input = rhs.Input;
  move_accepted         = rhs.move_accepted;
  AcceptanceProbability = rhs.AcceptanceProbability;
  localEnergy           = rhs.localEnergy;
  kineticEnergy         = rhs.kineticEnergy;
  potentialEnergy       = rhs.potentialEnergy;
  neEnergy              = rhs.neEnergy;
  eeEnergy              = rhs.eeEnergy;
  
  fwNormalization       = rhs.fwNormalization;
  fwR12                 = rhs.fwR12;
  fwR2                  = rhs.fwR2;
  fwiR12                = rhs.fwiR12;
  fwiR                  = rhs.fwiR;
  fwEnergy              = rhs.fwEnergy;
  fwKineticEnergy       = rhs.fwKineticEnergy;
  fwKineticEnergy_grad  = rhs.fwKineticEnergy_grad;
  fwPotentialEnergy     = rhs.fwPotentialEnergy;
  isCollectingFWResults = rhs.isCollectingFWResults;
  numFWSteps            = rhs.numFWSteps;
  
  if(Input->flags.nuclear_derivatives != "none"){
    fwNuclearForces.allocate(rhs.fwNuclearForces.dim1(),
                             rhs.fwNuclearForces.dim2());
    
    for (int d1=0; d1<fwNuclearForces.dim1(); d1++)
    {
      for (int d2=0; d2<fwNuclearForces.dim2(); d2++)
      {
        fwNuclearForces(d1,d2).allocate(rhs.fwNuclearForces.get(d1,d2).dim1(),2);
        fwNuclearForces(d1,d2) = rhs.fwNuclearForces.get(d1,d2);
      }
    }
  }
       
  distanceMovedAccepted = rhs.distanceMovedAccepted;
  dR2                   = rhs.dR2;
  R                     = rhs.R;
  walkerData            = rhs.walkerData;
}

void QMCWalker::initializePropagation(QMCWalkerData * &data,
                                      Array2D<double> * &rToCalc,
				      int iteration)
{
  this->iteration = iteration;
  createChildWalkers();
  forwardGreensFunction = TrialWalker->moveElectrons();
  data = & TrialWalker->walkerData;
  rToCalc = & TrialWalker->R;
}

void QMCWalker::processPropagation(QMCFunctions & QMF, bool writeConfigs)
{
  QMCGreensRatioComponent reverseGreensFunction =
    calculateReverseGreensFunction();
  double GreensFunctionRatio =
    (double)(reverseGreensFunction/forwardGreensFunction);

  if (IeeeMath::isNaN(GreensFunctionRatio))
    calculateMoveAcceptanceProbability(0.0);
  else
    calculateMoveAcceptanceProbability(GreensFunctionRatio);
  
  acceptOrRejectMove();
  reweight_walker();

  if(getWeight() == 0)
    return;

  if(move_accepted == false)
    {
      dR2                   = OriginalWalker->dR2;
    }

  if(writeConfigs)
    {
      /*
	If we want to be able to exactly recompose the statistics
	off of a cfgs file, then we'll need to print both configurations,
	since both are used in calculateObservables().
      */
      /*
      double p = TrialWalker->AcceptanceProbability;
      double q = 1.0 - p;

      TrialWalker->walkerData.writeConfigs(TrialWalker->R,p);
      OriginalWalker->walkerData.writeConfigs(OriginalWalker->R,q);
      */
    }

  calculateObservables();
  
  if(move_accepted == false)
    {
      R                     = OriginalWalker->R;
      walkerData            = OriginalWalker->walkerData;
      AcceptanceProbability = OriginalWalker->AcceptanceProbability;
    }
  
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Reinitializing a singular walker!!" << endl;
      initializeWalkerPosition(QMF);
    }

  if(writeConfigs)
    {
      /*
	Assuming the AcceptanceProbability is near one, this
	is probably more efficient than the option above.
       */
      walkerData.writeConfigs(R, AcceptanceProbability);
    }
}

QMCGreensRatioComponent QMCWalker::moveElectrons()
{
  // Make sure that dt is positive!
  if( Input->flags.dt <= 0 )
    {
      cerr << "ERROR: Negative dt value! (" << Input->flags.dt << ")" << endl;
      exit(0);
    }
  
  if(Input->flags.sampling_method == "no_importance_sampling")
      return moveElectronsNoImportanceSampling();
  else if(Input->flags.sampling_method == "importance_sampling" )
      return moveElectronsImportanceSampling();
  else if(Input->flags.sampling_method == "umrigar93_importance_sampling")
      return moveElectronsUmrigar93ImportanceSampling();
  else
    {
      cerr << "ERROR: Improper value for sampling_method set ("
	   << Input->flags.sampling_method << ")!" << endl;
      exit(0);
    }
    
  // should never reach this
  QMCGreensRatioComponent TRASH = QMCGreensRatioComponent();
  return TRASH;
}

QMCGreensRatioComponent QMCWalker::moveElectronsNoImportanceSampling()
{
  // Move the electrons of this walker
  double sigma = sqrt(Input->flags.dt);
  dR2 = 0;  
  
  // Don't use the QF  =>  R' = R + gauss_rn
  for(int i=0; i<R.dim1(); i++)
    {
      for(int j=0; j<R.dim2(); j++)
	{
	  // Add the randomness to the displacement
	  double drift = sigma*ran.gasdev();

	  // Calculate the square of the magnitude of the displacement
	  dR2 += drift * drift;

	  // Now update the R
	  R(i,j) += drift;
	}
    }
      
  return QMCGreensRatioComponent(1.0);
}

QMCGreensRatioComponent QMCWalker::moveElectronsImportanceSampling()
{
  Array2D<double> Displacement(R.dim1(),R.dim2());
  double tau    = Input->flags.dt;
  double sigma  = sqrt(tau);
  double greens = 0.0;
  
  // Use the QF => R' = R + dt * QF + gauss_rn
  // We have already counted for the factor of 2.0
  Displacement  = walkerData.modifiedGradPsiRatio;
  Displacement *= tau;
  dR2 = 0.0;  
  for(int i=0; i<R.dim1(); i++)
    {
      for(int j=0; j<R.dim2(); j++)
	{
	  // Add the randomness to the displacement
	  double drift = sigma*ran.gasdev();

	  // Add the randomness to the displacement
	  Displacement(i,j) += drift;

	  // Calculate the square of the magnitude of the displacement
	  dR2 += Displacement(i,j) * Displacement(i,j);

	  // Now update the R
	  R(i,j) += Displacement(i,j);

	  // Calculate the Green's function for the forward move
	  // The 'Quantum Force' term cancels out
	  greens += drift*drift;
	}
    }
    
  //k a^2 exp(c)
  double k = 1.0;
  double a = 2.0*3.14159265359*tau;
  double b = -0.5*R.dim1()*R.dim2();
  double c = -greens/(2.0*tau);
  
  return QMCGreensRatioComponent(k,a,b,c);
}

QMCGreensRatioComponent QMCWalker::moveElectronsUmrigar93ImportanceSampling()
{
  double tau = Input->flags.dt;

  dR2 = 0.0;
  Array2D<double> & Displacement = walkerData.modifiedGradPsiRatio;
  
  Array1D<double> zUnitVector(3);
  Array1D<double> radialUnitVector(3);
  Array1D<double> newPosition(3);
  
  QMCGreensRatioComponent GF(1.0);
  QMCGreensRatioComponent GaussianGF;
  QMCGreensRatioComponent SlaterGF;
  QMCGreensRatioComponent OneE;

  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {
      // Find the nearest nucleus to this electron
      int nearestNucleus = Input->Molecule.findClosestNucleusIndex(R,electron);
      double distanceFromNucleus = 0.0;
      
      // Calculate the unit vector in the Z direction
      for(int i=0; i<3; i++)
        {
          zUnitVector(i) = R(electron,i) -
                           Input->Molecule.Atom_Positions(nearestNucleus,i);
                           
          distanceFromNucleus += zUnitVector(i)*zUnitVector(i);
        }
        
      distanceFromNucleus = sqrt(distanceFromNucleus);
      zUnitVector /= distanceFromNucleus;
      
      
      // Decompose the modified QF into components in the z and radial
      // directions
      
      double zComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
          zComponentQF += zUnitVector(i)*Displacement(electron,i);
        
      double radialComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
        {
          radialUnitVector(i) = Displacement(electron,i) -
                                zComponentQF * zUnitVector(i);
                                
          radialComponentQF += radialUnitVector(i)*radialUnitVector(i);
        }
        
      radialComponentQF = sqrt(radialComponentQF);
      
      if (radialComponentQF > 1e-15)
	radialUnitVector /= radialComponentQF;
      
      // Calculate the gaussian drift components in cylindrical coordinates
      double zCoordinate = max(distanceFromNucleus + zComponentQF * tau,0.0);
                               
      double radialCoordinate = 2.0 * radialComponentQF * tau *
                           zCoordinate / ( distanceFromNucleus + zCoordinate );
                                
      // Calculate the hydrogen like exponential factor
      
      double expParam = Input->Molecule.Z(nearestNucleus);
      expParam = expParam*expParam + 1.0/tau;
      expParam = sqrt(expParam);
      
      // Calculate the probability of moving the electron with respect
      // to a gaussian type distribution with a QF drift or a hydrogenic
      // atom type slater function centered on the closest nucleus
      
      double probabilitySlaterTypeMove = 0.5 *
              MathFunctions::erfc( (distanceFromNucleus + zComponentQF * tau) /
                                                             sqrt(2.0 * tau) );
      if (probabilitySlaterTypeMove > 1.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 1.0;
	}
      if (probabilitySlaterTypeMove < 0.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 0.0;
	}
                                                              
      double probabilityGaussianTypeMove = 1.0 - probabilitySlaterTypeMove;
      // Randomly decide which electron moving method to use
      
      if( probabilityGaussianTypeMove > ran.unidev() )
        {
          // Gaussian Type Move
          
          // Particle is moved in the direction
          // Rnuc + radialCoordinate * radialUnitVector +
          //   zCoordinate * zUnitVector +
          //   gaussian random number with standard deviation sqrt(tau)
          for(int i=0; i<3; i++)
	    newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i)
	      + radialCoordinate * radialUnitVector(i) + zCoordinate * zUnitVector(i) +
	      sqrt(tau)*ran.gasdev();
            }
      else
        {
          // Slater Type Move
          
          for(int i=0; i<3; i++)
	    newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i);
            
          // add random part
          
          double r = ran.randomDistribution1()/(2.0*expParam);
          double phi = 2*3.14159265359*ran.unidev();
          double theta = ran.sindev();
          
          newPosition(0) += r*sin(theta)*cos(phi);
          newPosition(1) += r*sin(theta)*sin(phi);
          newPosition(2) += r*cos(theta);
        }
        
      // Update the greens function
      double distance1Sq = 0.0;
      double temp;
      
      for(int i=0; i<3; i++)
        {
          temp = newPosition(i) -
                        ( Input->Molecule.Atom_Positions(nearestNucleus,i) +
       radialCoordinate * radialUnitVector(i) + zCoordinate * zUnitVector(i) );
                          
          distance1Sq += temp*temp;
        }
        
      double ga = 2*3.14159265359*tau;
      double gb = -1.5;
      double gc = -distance1Sq/(2*tau);
      
      GaussianGF=QMCGreensRatioComponent(probabilityGaussianTypeMove,ga,gb,gc);
      
      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          temp = newPosition(i) -
                              Input->Molecule.Atom_Positions(nearestNucleus,i);
                        
          distance2Sq += temp*temp;
        }
        
      double sk = probabilitySlaterTypeMove/3.14159265359;
      double sa = expParam;
      double sb = 3.0;
      double sc = -2.0*expParam*sqrt(distance2Sq);
      
      SlaterGF = QMCGreensRatioComponent(sk,sa,sb,sc);
      
      // The addition and multiplication operations here have been causing a
      // lot of problems.
      OneE = SlaterGF + GaussianGF;
      
      GF *= OneE;
      
      // Update the distance attempted to move squared and move the electron.
      
      for(int i=0; i<3; i++)
        {
          temp = newPosition(i) - R(electron,i);
          dR2 += temp*temp;
          R(electron,i) = newPosition(i);
        }
    }
  return GF;
}

QMCGreensRatioComponent QMCWalker::calculateReverseGreensFunction()
{
  // Make sure that dt is positive!
  if( Input->flags.dt <= 0 )
    {
      cerr << "ERROR: Negative dt value! (" << Input->flags.dt << ")" << endl;
      exit(0);
    }
    
  if(Input->flags.sampling_method == "no_importance_sampling")
      return calculateReverseGreensFunctionNoImportanceSampling();
  else if(Input->flags.sampling_method == "importance_sampling" )
      return calculateReverseGreensFunctionImportanceSampling();
  else if(Input->flags.sampling_method == "umrigar93_importance_sampling")
      return calculateReverseGreensFunctionUmrigar93ImportanceSampling();
  else
    {
      cerr << "ERROR: Improper value for sampling_method set ("
	   << Input->flags.sampling_method << ")!" << endl;
      exit(0);
    }
    
  // should never reach this
  
  QMCGreensRatioComponent TRASH = QMCGreensRatioComponent();
  return TRASH;
}

QMCGreensRatioComponent \
QMCWalker::calculateReverseGreensFunctionNoImportanceSampling()
{
  return QMCGreensRatioComponent(1.0);
}

QMCGreensRatioComponent \
QMCWalker::calculateReverseGreensFunctionImportanceSampling()
{
  double tau = Input->flags.dt;  
  double greens = 0.0;
  
  for(int i=0; i<R.dim1(); i++)
    {
      for(int j=0; j<R.dim2(); j++)
	{
	  double temp = OriginalWalker->R(i,j) - TrialWalker->R(i,j) -
	    tau*TrialWalker->walkerData.modifiedGradPsiRatio(i,j);
	  
	  greens += temp*temp;
	}
    }

  // k * a^b * exp(c)
  double k = 1.0;
  double a = 2.0*3.14159265359*tau;
  double b = -0.5*R.dim1()*R.dim2();
  double c = -greens/(2.0*tau);

  return QMCGreensRatioComponent(k,a,b,c);
}

QMCGreensRatioComponent \
QMCWalker::calculateReverseGreensFunctionUmrigar93ImportanceSampling()
{
  double tau = Input->flags.dt;
  Array1D<double> zUnitVector(3);
  Array1D<double> radialUnitVector(3);

  QMCGreensRatioComponent GF(1.0);
  QMCGreensRatioComponent GaussianGF;
  QMCGreensRatioComponent SlaterGF;
  QMCGreensRatioComponent OneE;

  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {
      // Find the nearest nucleus to this electron
      int nearestNucleus =
        Input->Molecule.findClosestNucleusIndex(TrialWalker->R,electron);
        
      // Calculate the unit vector in the Z direction
      double distanceFromNucleus = 0.0;
      
      for(int i=0; i<3; i++)
        {
          zUnitVector(i) = TrialWalker->R(electron,i) -
                           Input->Molecule.Atom_Positions(nearestNucleus,i);
                           
          distanceFromNucleus += zUnitVector(i)*zUnitVector(i);
        }
        
      distanceFromNucleus = sqrt(distanceFromNucleus);
      zUnitVector /= distanceFromNucleus;
      
      // Decompose the modified QF into components in the z and radial
      // directions
      
      double zComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
          zComponentQF += zUnitVector(i)*
                          TrialWalker->walkerData.modifiedGradPsiRatio(electron,i);
        
      double radialComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
        {
          radialUnitVector(i) =
            TrialWalker->walkerData.modifiedGradPsiRatio(electron,i) -
            zComponentQF * zUnitVector(i);
            
          radialComponentQF += radialUnitVector(i)*radialUnitVector(i);
        }
        
      radialComponentQF = sqrt(radialComponentQF);
      
      if (radialComponentQF > 1e-15)
	radialUnitVector /= radialComponentQF;
      
      // Calculate the gaussian drift components in cylindrical coordinates
      double zCoordinate = max(distanceFromNucleus + zComponentQF * tau,0.0);
                               
      double radialCoordinate = 2.0 * radialComponentQF * tau *
                           zCoordinate / ( distanceFromNucleus + zCoordinate );
                                
      // Calculate the hydrogen like exponential factor
      
      double expParam = Input->Molecule.Z(nearestNucleus);
      expParam = expParam*expParam + 1.0/tau;
      expParam = sqrt(expParam);
      
      // Calculate the probability of moving the electron with respect
      // to a gaussian type distribution with a QF drift or a hydrogenic
      // atom type slater function centered on the closest nucleus
      
      double probabilitySlaterTypeMove = 0.5 *
              MathFunctions::erfc( (distanceFromNucleus + zComponentQF * tau) /
                                                             sqrt(2.0 * tau) );

      if (probabilitySlaterTypeMove > 1.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 1.0;
	}
      if (probabilitySlaterTypeMove < 0.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 0.0;
	}

      double probabilityGaussianTypeMove = 1.0 - probabilitySlaterTypeMove;
      
      // Update the greens function
      
      double distance1Sq = 0.0;
      double temp;
      
      for(int i=0; i<3; i++)
        {
          temp = OriginalWalker->R(electron,i) -
	    ( Input->Molecule.Atom_Positions(nearestNucleus,i) +
	      radialCoordinate * radialUnitVector(i) + zCoordinate * zUnitVector(i) );
                          
          distance1Sq += temp*temp;
        }
        
      double ga = 2*3.14159265359*tau;
      double gb = -1.5;
      double gc = -distance1Sq/(2*tau);
      
      GaussianGF=QMCGreensRatioComponent(probabilityGaussianTypeMove,ga,gb,gc);
      
      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          temp = OriginalWalker->R(electron,i) - \
                        Input->Molecule.Atom_Positions(nearestNucleus,i);
                        
          distance2Sq += temp*temp;
        }
        
      double sk = probabilitySlaterTypeMove/3.14159265359;
      double sa = expParam;
      double sb = 3.0;
      double sc = -2.0*expParam*sqrt(distance2Sq);
      
      SlaterGF = QMCGreensRatioComponent(sk,sa,sb,sc);

      // The addition and multiplication operations here have been causing a
      // a lot of problems.
      OneE = SlaterGF + GaussianGF;
      
      GF *= OneE;
    }
  return GF;
}

void QMCWalker::reweight_walker()
{
  dW = 0.0;
  double S_trial = 0.0;
  double S_original = 0.0;
  double trialEnergy    = TrialWalker->walkerData.localEnergy;
  double originalEnergy = OriginalWalker->walkerData.localEnergy;
      
  // determine the weighting factor dW so that the new weight = weight*dW
  
  bool weightIsNaN = false;

  if( Input->flags.run_type == "variational" )
      // Keep weights constant for VMC
      dW = 1.0;

  else
    {     
      if( IeeeMath::isNaN(trialEnergy) )
	{
	  cerr << "WARNING: trial energy = ";
	  cerr << TrialWalker->walkerData.localEnergy << "!" << endl;
	  weightIsNaN = true;
	}
      else if( Input->flags.energy_modification_type == "none" )
        {
          S_trial    = Input->flags.energy_trial - trialEnergy;
          S_original = Input->flags.energy_trial - originalEnergy;
        }
      else
        {
          Array2D<double> * tMGPR =
            & TrialWalker->walkerData.modifiedGradPsiRatio;
          Array2D<double> * tGPR = & TrialWalker->walkerData.gradPsiRatio;
          Array2D<double> * oMGPR =
            & OriginalWalker->walkerData.modifiedGradPsiRatio;
          Array2D<double> * oGPR = & OriginalWalker->walkerData.gradPsiRatio;
          
          double lengthGradTrialModified =
            sqrt((*tMGPR).dotAllElectrons(*tMGPR));
          double lengthGradTrialUnmodified =
            sqrt((*tGPR).dotAllElectrons(*tGPR));
          double lengthGradOriginalModified =
            sqrt((*oMGPR).dotAllElectrons(*oMGPR));
          double lengthGradOriginalUnmodified =
            sqrt((*oGPR).dotAllElectrons(*oGPR));

	  if(Input->flags.energy_modification_type=="modified_umrigar93")
            {
              S_trial = (Input->flags.energy_trial - trialEnergy)*
                        (lengthGradTrialModified/lengthGradTrialUnmodified);
                        
              S_original = (Input->flags.energy_trial - originalEnergy) *
                           (lengthGradOriginalModified/lengthGradOriginalUnmodified);
            }
          else if( Input->flags.energy_modification_type == "umrigar93" )
            {
              S_trial =
		(Input->flags.energy_trial - Input->flags.energy_estimated)
		+(Input->flags.energy_estimated - trialEnergy)*
                (lengthGradTrialModified/lengthGradTrialUnmodified);
                
              S_original =
		(Input->flags.energy_trial - Input->flags.energy_estimated)
		+(Input->flags.energy_estimated - originalEnergy) *
                (lengthGradOriginalModified/lengthGradOriginalUnmodified);
            }
          else
            {
              cerr << "ERROR: unknown energy modification method!" << endl;
              exit(0);
            }
        }

      if (IeeeMath::isNaN(S_trial) || IeeeMath::isNaN(S_original))
	{
	  cerr << "Error: S_trial          = " << S_trial << endl;
	  cerr << "       S_original       = " << S_original << endl;
	  cerr << "       energy_trial     = " << Input->flags.energy_trial << endl;
	  cerr << "       energy_estimated = " << Input->flags.energy_estimated << endl;
	  cerr << "       trialEnergy      = " <<    TrialWalker->walkerData.localEnergy << endl;
	  cerr << "       originalEnergy   = " << OriginalWalker->walkerData.localEnergy << endl;
	}

      if( Input->flags.walker_reweighting_method == "simple_symmetric" )
        {
          // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 8
          // This is the "classical" reweighting factor most people use.
          // umrigar says this has larger timestep and statistical errors than
          // umrigar93_probability_weighted
          double temp = 0.5*(S_trial+S_original)*Input->flags.dt_effective;
	  
	  if (IeeeMath::isNaN(temp) || weightIsNaN == true)
	    {
	      cerr << "WARNING: dW = exp(" << temp << ")" << endl;
	      cerr << "Walker's weight is being set to zero so that it does"
		   << " not ruin the whole calculation." << endl;
	      dW = 0.0;
	    }
	  else
          dW = exp(temp);
        }
      else if( Input->flags.walker_reweighting_method == "simple_asymmetric" )
        {
          // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 23
          // supposedly between simple_symmetric and
          // umrigar93_probability_weighted in terms of its timestep error and
          // statistical performance
          double temp;
          if( move_accepted )
            temp = 0.5*(S_trial+S_original)*Input->flags.dt_effective;
          else
            temp = S_original*Input->flags.dt_effective;

	  if (IeeeMath::isNaN(temp) || weightIsNaN == true)
	    {
	      cerr << "WARNING: dW = exp(" << temp << ")" << endl;
	      cerr << "Walker's weight is being set to zero so that it does"
		   << " not ruin the whole calculation." << endl;
	      dW = 0.0;
	    }
	  else
          dW = exp(temp);
        }
      else if( Input->flags.walker_reweighting_method ==
               "umrigar93_probability_weighted" )
        {
          // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993
          // Umrigar claims this has a small time step error and a small
          // statistical error compared to simple_asymmetric and
          // simple_symmetric
          double p = TrialWalker->getAcceptanceProbability();
	  double q = 1.0 - p;
          
	  double temp = (p*0.5*(S_original+S_trial) + q*S_original) *
                        Input->flags.dt_effective;

	  if (IeeeMath::isNaN(temp) || weightIsNaN == true)
	    {
	      cerr << "WARNING: dW = exp(" << temp << ")" << endl;
	      cerr << "Walker's weight is being set to zero so that it does"
		   << " not ruin the whole calculation." << endl;
	      cerr << " p = " << p << "; q = " << q << endl;
	      dW = 0.0;
	    }
	  else
	    dW = exp(temp);
	}
      else
	{
	  cerr << "ERROR: unknown reweighting method!" << endl;
          exit(0);
        }
    }


  // now set the weight of the walker
  setWeight( getWeight() * dW );

  if(getWeight() <= 0.0 || Input->flags.run_type == "variational")
    return;

  double rel_diff;
  if(move_accepted)
    {
      rel_diff = fabs( (TrialWalker->walkerData.localEnergy -
			Input->flags.energy_estimated_original)/
		       Input->flags.energy_estimated_original);
    }
  else
    {
      rel_diff = fabs( (OriginalWalker->walkerData.localEnergy -
			Input->flags.energy_estimated_original)/
		       Input->flags.energy_estimated_original);
    }

  if(iteration < -100)
    {
      /*
	We're equilibrating when iteration < 0.

	If we're equilibrating, we can be more aggressive about
	handling bad walkers.
      */
      
      // we probably don't need warnings for the first couple steps
      int steps = iteration + Input->flags.equilibration_steps;
      // the number of steps until warnings are printed
      int min_steps = 5;

      if(getWeight() > 1.5*Input->flags.branching_threshold)
	{
	  if(steps > min_steps)
	    {
	      cerr << "WARNING: Deleting heavy walker " << ID();
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}

      if(getWeight() < 0.1)
	{
	  if(steps > min_steps)
	    {
	      cerr << "WARNING: Deleting light walker " << ID();
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}

      if(dW > 1.5 || dW < 0.5)
	{
	  if(steps > min_steps)
	    {

	      if(dW > 2.0 || dW < 0.25)
		{
		  double p = TrialWalker->getAcceptanceProbability();
		  double q = 1.0 - p;
		  
		  cerr << "WARNING: Deleting walker with bad dW " << ID();
		  cerr << "       p = " << p << "; q = " << q;
		  cerr << "       energy_trial     = " << Input->flags.energy_trial << endl;
		  cerr << "       energy_estimated = " << Input->flags.energy_estimated << endl;
		  cerr << "       S_trial          = " << S_trial << endl;
		  cerr << "       S_original       = " << S_original << endl;
		  cerr << "TrialWalker:" << endl << TrialWalker->walkerData;
		  cerr << "OriginalWalker:" << endl << OriginalWalker->walkerData;
		  cerr << endl;
		}
	      else
		{
		  //cerr << "WARNING: Deleting fast growth walker " << ID();
		}
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}

      double rel_cutoff = Input->flags.rel_cutoff;
      if(rel_diff > rel_cutoff && steps > min_steps)
	{
	  /*
	    A lot of walkers when initialized will have bad energies. So we
	    want to give them a couple of chances to move first.

	    Since rel_diff uses the SCF energy, it will bias the energy
	    toward zero, so we want a tradeoff between biasing the energy
	    and not allowing horrible energies into the average.

	    rel_diff > 1 corresponds to a positive energy, clearly unacceptable.
	    Is there any argument suggesting that a Local Energy can never be +'ve?
	   */
	  if(rel_diff > rel_cutoff && steps > 2*min_steps){
	    cerr << "WARNING: Deleting walker with bad energy " << ID();
	    cerr.flush();
	  }
	  setWeight(0.0);
	  return;
	}

      int age_cutoff = 30;
      if(age >= age_cutoff)
	{
	  // I don't see how a warning here could be very useful...
	  if(steps > age_cutoff && false)
	    {
	      cerr << "WARNING: Deleting aged walker " << ID();
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}
    }
  else
    {
      /*
	We're not equillibrating anymore, so we want to be more
	careful. We'll only delete walkers we know will mess up
	the calculation.

	How aggressive can we be?
      */
      if(getWeight() > 50.0)
	{
	  cerr << "ERROR: Deleting heavy walker " << ID();
	  cerr.flush();
	  setWeight(0.0);
	  return;
	}

      if(rel_diff > 5.0)
	{
	  /*
	    Maybe this warning is unnecessary since it won't be included
	    in the energy anyway....
	   */
	  cerr << "ERROR: Deleting walker with bad energy " << ID();
	  cerr.flush();
	  setWeight(0.0);
	  return;
	}

      if(dW > 4.0)
	{
	  stringstream strm;
	  double p = TrialWalker->getAcceptanceProbability();
	  double q = 1.0 - p;
	  
	  strm << "ERROR: Deleting fast growth walker " << ID();
	  strm << "       p = " << p << "; q = " << q;
	  strm << "       energy_trial     = " << Input->flags.energy_trial << endl;
	  strm << "       energy_estimated = " << Input->flags.energy_estimated << endl;
	  strm << "       S_trial          = " << S_trial << endl;
	  strm << "       S_original       = " << S_original << endl;
	  strm << "TrialWalker:" << endl << TrialWalker->walkerData;
	  strm << "OriginalWalker:" << endl << OriginalWalker->walkerData;
	  strm << endl;
	  cerr << strm.str();
	  cerr.flush();
	  setWeight(0.0);
	  return;
	}
    }
}

bool QMCWalker::branchRecommended()
{
  /*
    This function will be queried before a branch, meaning the walker is a
    candidate for branching. Since branching is purely an
    efficiency issue, we can do whatever we want and shouldn't have to worry about
    ruining detailed balance.

    We're equilibrating when iteration < 0, so we can be more aggressive in penalizing
    bad walkers... On the other hand, since we're deleting a lot of walkers during
    equilibration, maybe we want to make branching easier to replenish the loses.

    The difficulty is trying to get some warnings, but not too many. If we didn't
    delete all the bad walkers during equilibration, then "badness" can propagate
    through the calculation.
    
    By the time this function is called, OriginalWalker == TrialWalker
  */
  /*
    Some ideas to try:
    1) No branching until you've moved at least once.
    2) No branching if LE - TrialEnergy is outside some range (relative or abs?).
    3) A walker can't contribute to the energy average if its weight is too large
    4) Walkers can only be deleted if they haven't moved (and then raise the threshold?)
    5) 
   */
  bool shouldRecommend = true;
  bool shouldWarn = false;
  int aged = age - 10;

  if(age > 4)
    {
      shouldRecommend = false;
      if(age >= Input->flags.old_walker_acceptance_parameter && age%10 == 0)
	shouldWarn = true;
    }
  if(dW > 1.5)
    {
      shouldRecommend = false;
      if(iteration%10 == 0 && getWeight() > 2.0*Input->flags.branching_threshold)
	shouldWarn = true;
      if(dW > 2.0)
	shouldWarn = true;
    }
  if(getWeight() > 2.0*Input->flags.branching_threshold)
    {
      shouldRecommend = false;
      if(iteration%10 == 0 &&
	 getWeight() > 2.0*Input->flags.branching_threshold &&
	 dW > 1.0)
	shouldWarn = true;
    }

  double virial = -TrialWalker->walkerData.potentialEnergy/TrialWalker->walkerData.kineticEnergy;

  // Relatively large virial ratios do not appear to correlate with bad E_L.
  // It's hard to know how this would be a good indicator given rel_diff
  if(fabs(virial) < 0.1)
    {
      shouldRecommend = false;
      //only warn when we arrive at a bad spot
      if(age == 0)
	shouldWarn = true;
    }

  double rel_diff = fabs( (TrialWalker->walkerData.localEnergy -
			   Input->flags.energy_estimated_original)/
			  Input->flags.energy_estimated_original);
  if(rel_diff > 1.0)
    {
      shouldRecommend = false;
      if(age == 0)
	shouldWarn = true;
    }

  if(shouldWarn && !shouldRecommend && iteration > 0)
    {
      cerr << "WARNING: Not recommending a branch for walker " << ID();
      cerr.flush();
    }

  return shouldRecommend;
}

string QMCWalker::ID()
{
  /*
    The output only depends upon move_accepted if
    ID() is called before the processPropagation function completes
    since the last thing processPropagation does is 
    syncronize the Trial and Original walkerData.
   */
  QMCWalkerData * wd;
  if(move_accepted)
    wd = & TrialWalker->walkerData;
  else
    wd = & OriginalWalker->walkerData;

  double virial = 0;
  if(fabs(wd->kineticEnergy) > 1e-30)
    virial = -wd->potentialEnergy/wd->kineticEnergy;

  double rel_diff = fabs( (wd->localEnergy -
		    Input->flags.energy_estimated_original)/
		   Input->flags.energy_estimated_original);

  stringstream id;

  id << "(" << genealogy[0];
  for(int i=1; i<numAncestors; i++)
    id << "<" << genealogy[i];

  id << "::" << Input->flags.my_rank << ")" << endl;

  id << "     weight = ";
  id.precision(15);
  id.width(20);
  id << getWeight() << endl;

  id << "         dW = ";
  id.precision(15);
  id.width(20);
  id << dW << endl;

  id << "     energy = ";
  id.precision(15);
  id.width(20);
  id << wd->localEnergy << endl;

  /*
  id << "  mod.ratio = ";
  id.precision(15);
  id.width(20);
  id << wd->modificationRatio << endl;
  */

  id << "     virial = ";
  id.precision(15);
  id.width(20);
  id << virial << endl;

  id << "   rel_diff = ";
  id.precision(15);
  id.width(20);
  id << rel_diff << endl;

  id << "        age = ";
  id.width(20);
  id << age << endl;
  
  id << "       iter = ";
  id.width(20);
  id << iteration << endl;
  return id.str();
}

void QMCWalker::calculateMoveAcceptanceProbability(double GreensRatio)
{
  // This tells us the probability of accepting or rejecting a proposed move
  
  double PsiRatio = TrialWalker->walkerData.psi/OriginalWalker->walkerData.psi;
  
  // calculate the probability of accepting the trial move
  double p = PsiRatio * PsiRatio * GreensRatio;
  
  if( !(IeeeMath::isNaN(p)))
    {
      // 1.1 ^ 1000 ~= 1e41
      if(age - Input->flags.old_walker_acceptance_parameter > 1000)
	{
	  //setWeight(0.0);
	  p = 1.0;
	} 
      else if(age - Input->flags.old_walker_acceptance_parameter > 500)
	{
	  p = 1.0;
	} 
      else if(age - Input->flags.old_walker_acceptance_parameter > 0)
	{
	  p *= pow(1.1,age-Input->flags.old_walker_acceptance_parameter);
	}
      
    } else {
      // if the aratio is NaN then reject the move
      cerr << "WARNING: Rejecting trial walker with NaN p!" << endl;
      cerr << "         PsiRatio    = " << PsiRatio << endl;
      cerr << "         GreensRatio = " << GreensRatio << endl;
      p = 0.0;
      // The energies are assigned zero because the later multiplication by 
      // p=0 doesn't result in 0 contribution.
      TrialWalker->walkerData.zero();
    }
  
  // The particular NaN that this is correcting for is not revealed by isinf or
  // isnan...
  double kineticEnergy = TrialWalker->walkerData.kineticEnergy;
  if( IeeeMath::isNaN(kineticEnergy) )
    {
      cerr << "WARNING: Rejecting trial walker with NaN kinetic energy!" << endl;
      p = 0.0;
      TrialWalker->walkerData.zero();
    }

  double potentialEnergy = TrialWalker->walkerData.potentialEnergy;
  if( IeeeMath::isNaN(potentialEnergy) )
    {
      cerr << "WARNING: Rejecting trial walker with NaN potential energy!" << endl;
      p = 0.0;
      TrialWalker->walkerData.zero();
    }
    
  // if the trial position is singular reject the move
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Rejecting singular trial walker!" << endl;
      p = 0.0;
      TrialWalker->walkerData.zero();
    }
    
  // Fixed Node Condition
  if( PsiRatio < 0 && Input->flags.run_type == "diffusion" )
    {
      // Fixed Node Condition
      p = 0.0;
    }
    
  // Limit the probability to 1
  if( p > 1.0)
    {
      p = 1.0;
    }
    
  // Set the probabilities of taking the trial move or keeping the current
  // move
  TrialWalker->setAcceptanceProbability(p);
}

void QMCWalker::acceptOrRejectMove()
{
  if( TrialWalker->getAcceptanceProbability() > ran.unidev() )
    {
      // accept the move
      ageMoved = age;
      age = 0;
      move_accepted = true;
    }
  else
    {
      // reject the move
      //age measures the number of steps for which the move was not accepted
      age++;
      move_accepted = false;
    }

  int aged = age - Input->flags.old_walker_acceptance_parameter;
  if( iteration > 0)
    {
      if(ageMoved > Input->flags.old_walker_acceptance_parameter + 15)
	{
	  cerr << "WARNING: Old walker moved after " << ageMoved << " iterations " << ID();
	}
      /*
      if(aged >= 50 && aged%50 == 0)
	{
	  cerr << "WARNING: Aged walker " << ID();
	}
      */
    }
}

void QMCWalker::createChildWalkers()
{
  /*
    We only ever need 2 walker objects, so AGA modified the code
    so that we don't use 3 walker objects anymore. In order to preserve
    the rest of the code, the TrialWalker pointer just points to "this" now,
    this results in fewer memory copies (since the move is likely to be accepted)
    than if OriginalWalker pointed to "this".

    Correspondingly, I eliminated redundant data.
  */
  if( OriginalWalker == 0 )
    {
      OriginalWalker = new QMCWalker();
    }

  //We create a copy of the data we need
  OriginalWalker->R                     = R;
  OriginalWalker->walkerData            = walkerData;
  OriginalWalker->dR2                   = dR2;
  OriginalWalker->AcceptanceProbability = AcceptanceProbability;

  //And use a psuedonym
  TrialWalker = this;
}

void QMCWalker::initialize(QMCInput *INPUT)
{
  Input = INPUT;

  int numFW = Input->flags.future_walking.size();
  int numElectrons  = Input->WF.getNumberElectrons();
  int numDimensions;
  int numNucForceDim1;
  int numNucForceDim2;
  
  if(Input->flags.trial_function_type == "restricted" ||
     Input->flags.trial_function_type == "unrestricted"){     
    numDimensions = 3;
  } else {
    //We're hijacking the Nbasisfunc parameter input
    numDimensions = Input->flags.Nbasisfunc;
  }

  if(Input->flags.nuclear_derivatives != "none")
    {
      if(Input->flags.nuclear_derivatives != "bin_force_density")
	{
	  numNucForceDim1 = Input->Molecule.getNumberAtoms();
	  numNucForceDim2 = 3;
	} else {
	  numNucForceDim1 = QMCNuclearForces::getNumBins();
	  numNucForceDim2 = 1;
	}
      
      fwNuclearForces.allocate(numNucForceDim1,numNucForceDim2);
      
      for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
	for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
	  fwNuclearForces(d1,d2).allocate(numFW,2);
    }

  R.allocate(numElectrons,numDimensions);

  walkerData.initialize(Input,numDimensions,numNucForceDim1,numNucForceDim2);

  numFWSteps.resize(numFW);
  
  isCollectingFWResults.allocate(numFW,2);
  fwNormalization.allocate(numFW,2);
  fwR12.allocate(numFW,2);
  fwR2.allocate(numFW,2);
  fwiR12.allocate(numFW,2);
  fwiR.allocate(numFW,2);
  fwEnergy.allocate(numFW,2);
  fwKineticEnergy.allocate(numFW,2);
  fwKineticEnergy_grad.allocate(numFW,2);
  fwPotentialEnergy.allocate(numFW,2);
  resetFutureWalking();
        
  //initialize acceptance probability
  setAcceptanceProbability(0.0);
}

void QMCWalker::calculateElectronDensities(double max_pair_distance, double dr,
                          Array1D<double> &pll_spin, Array1D<double> &opp_spin,
		                     Array1D< Array1D<double> > &alpha_density,
                                      Array1D< Array1D<double> > &beta_density)
{
  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();
  
  double dist = 0.0;
  int index = 0;

  // We calculate the distance between each same spin pair and record them in 
  // the histogram.
  if (nalpha > 1)
    {
      for (int i=0; i<nalpha-1; i++)
        for (int j=i+1; j<nalpha; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
	    if (dist < max_pair_distance)
	      {
		index = int(dist/dr);
		pll_spin(index) += weight;
	      }
          }
    }
    
  if (nbeta > 1)
    {
      for (int i=nalpha; i<nalpha+nbeta-1; i++)
        for (int j=i+1; j<nalpha+nbeta; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
	    if (dist < max_pair_distance)
	      {
		index = int(dist/dr);
		pll_spin(index) += weight;
	      }
          }
    }

  // We calculate the distance between each opposite spin pair and record them
  // in the histogram.
  if (nalpha > 0 && nbeta > 0)
    {
      for (int i=0; i<nalpha; i++)
        for (int j=nalpha; j<nalpha+nbeta; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
	    if (dist < max_pair_distance)
	      {
		index = int(dist/dr);
		opp_spin(index) += weight;
	      }
          }
    }

  for (int i=0; i<Input->flags.Natoms; i++)
    {
      string NucleusType = Input->Molecule.Atom_Labels(i);

      int nuc_index = -1;
      for (int j=0; j<Input->Molecule.NucleiTypes.dim1(); j++)
	if (Input->Molecule.NucleiTypes(j) == NucleusType)
	  {
	    nuc_index = j;
	    break;
	  }

      double nuc_x = Input->Molecule.Atom_Positions(i,0);
      double nuc_y = Input->Molecule.Atom_Positions(i,1);
      double nuc_z = Input->Molecule.Atom_Positions(i,2);
      
      for (int k=0; k<nalpha; k++)
	{
	  dist = sqrt( (R(k,0) - nuc_x) * (R(k,0) - nuc_x) + 
		       (R(k,1) - nuc_y) * (R(k,1) - nuc_y) +
		       (R(k,2) - nuc_z) * (R(k,2) - nuc_z) );
	  if (dist < max_pair_distance)
	    {
	      index = int(dist/dr);
	      (alpha_density(nuc_index))(index) += weight;
	    }
	}

      for (int l=nalpha; l<nalpha+nbeta; l++)
	{
	  dist = sqrt( (R(l,0) - nuc_x) * (R(l,0) - nuc_x) +
		       (R(l,1) - nuc_y) * (R(l,1) - nuc_y) +
		       (R(l,2) - nuc_z) * (R(l,2) - nuc_z) );
	  if (dist < max_pair_distance)
	    {
	      index = int(dist/dr);
	      (beta_density(nuc_index))(index) += weight;
	    }
	}
    }
}

void QMCWalker::calculatePllCorrelationDiagram(int coord, double min, 
		    double max, Array1D< Array1D<double> > &CorrelationDiagram)
{
  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();
  
  double dr = (max-min)/CorrelationDiagram.dim1();
  int index1 = -1;
  int index2 = -1;

  for (int i=0; i<nalpha-1; i++)
    for (int j=i+1; j<nalpha; j++)
      if ( (R(i,coord) >= min) && (R(i,coord) <=max) )
	if ( (R(j,coord) >= min) && (R(j,coord) <= max) )
	  {
	    index1 = int((R(i,coord)-min)/dr);
	    index2 = int((R(j,coord)-min)/dr);
	    (CorrelationDiagram(index1))(index2) += weight;
	  }

  for (int i=nalpha; i<nalpha+nbeta-1; i++)
    for (int j=i+1; j<nalpha+nbeta; j++)
      if ( (R(i,coord) >= min) && (R(i,coord) <=max) )
	if ( (R(j,coord) >= min) && (R(j,coord) <= max) )
	  {
	    index1 = int((R(i,coord)-min)/dr);
	    index2 = int((R(j,coord)-min)/dr);
	    (CorrelationDiagram(index1))(index2) += weight;
	  }
}

void QMCWalker::calculateOppCorrelationDiagram(int coord, double min, 
		    double max, Array1D< Array1D<double> > &CorrelationDiagram)
{
  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();
  
  double dr = (max-min)/CorrelationDiagram.dim1();
  int index1 = 0;
  int index2 = 0;

  for (int i=0; i<nalpha; i++)
    for (int j=nalpha; j<nalpha+nbeta; j++)
      if ( (R(i,coord) >= min) && (R(i,coord) <=max) )
	if ( (R(j,coord) >= min) && (R(j,coord) <= max) )
	  {
	    index1 = int((R(i,coord)-min)/dr);
	    index2 = int((R(j,coord)-min)/dr);
	    (CorrelationDiagram(index1))(index2) += weight;
	  }
}

void QMCWalker::toXML(ostream& strm)
{
  strm << "<QMCWalker>" << endl;
  strm << "\t<Position>" <<endl;
  for(int ep=0; ep<R.dim1(); ep++)
    {
      strm << "\t\t";
      for(int j=0;j<R.dim2();j++)
        {
          strm << R(ep,j) << "    ";
        }
      strm << endl;
    }
  strm << "\t</Position>" << endl;
  strm << "\t<Weight>\n\t\t" << getWeight() << "\n\t</Weight>" << endl;
  strm << "\t<Age>\n\t\t" << getAge() << "\n\t</Age>" << endl;
  strm << "\t<Elocal> \n\t\t" << walkerData.psi
  << "\n\t</Elocal>" <<endl;
  strm << "</QMCWalker>\n" << endl;
}

void QMCWalker::readXML(istream& strm, QMCFunctions & QMF)
{
  string temp;
  strm >> temp;
  
  // Read position
  strm >> temp;
  
  for(int ep=0; ep<R.dim1(); ep++)
    {
      for(int j=0;j<R.dim2();j++)
        {
          strm >> R(ep,j);
        }
    }
    
  strm >> temp;
  
  // Read weight
  strm >> temp >> temp;
  weight = atof(temp.c_str());
  strm >> temp;
  
  
  // Read age
  strm >> temp >> temp;
  age = atoi(temp.c_str());
  
  // Read the energy and don't save it
  strm >> temp >> temp >> temp >> temp >> temp;
  
  QMF.evaluate(R,walkerData);
  newID();
}

void QMCWalker::initializeWalkerPosition(QMCFunctions & QMF)
{
  QMCInitializeWalker * IW =
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
        Input->flags.walker_initialization_method);
  
  R = IW->initializeWalkerPosition();
  QMF.evaluate(R,walkerData);
  
  int initilization_try = 1;
  while( isSingular() )
    {
      cerr << "Regenerating Walker..." << endl;
      
      if( initilization_try > 10 )
        {
          cerr << "ERROR: 10 consecutive singular configurations while "
          << "trying to initilize walker!" << endl;
          exit(0);
        }
        
      R = IW->initializeWalkerPosition();
      QMF.evaluate(R,walkerData);
      initilization_try++;
    }
    
  delete IW;
  IW = 0;
}

double QMCWalker::getWeight()
{
  return weight;
}

void QMCWalker::setWeight(double val)
{
  weight = val;
}

int QMCWalker::getAge()
{
  return age;
}

double QMCWalker::getAcceptanceProbability()
{
  return AcceptanceProbability;
}

void QMCWalker::setAcceptanceProbability(double val)
{
  AcceptanceProbability = val;
}

Array2D<double> * QMCWalker::getR()
{
  return &R;
}

bool QMCWalker::setR(Array2D<double>& temp_R)
{
  bool ok = true;
  for(int i=0; i<temp_R.dim1(); i++)
    for(int j=0; j<temp_R.dim2(); j++)
      if(IeeeMath::isNaN(temp_R(i,j)) || temp_R(i,j) == 0 || fabs(temp_R(i,j)) > 500.0 )
	{
	  cout << "(" << i << "," << j << ") = " << temp_R(i,j) << endl;
	  ok = false;
	}
  if(ok) R = temp_R;
  return ok;
}

QMCWalkerData* QMCWalker::getWalkerData()
{
  return &walkerData;
}

void QMCWalker::calculateObservables()
{
  double p = TrialWalker->AcceptanceProbability;
  double q = 1.0 - p;

  // Calculate the Energy ...
  localEnergy = p * TrialWalker->walkerData.localEnergy +
                q * OriginalWalker->walkerData.localEnergy;
  
  // Calculate the kinetic energy...
  kineticEnergy = p * TrialWalker->walkerData.kineticEnergy +
                  q * OriginalWalker->walkerData.kineticEnergy;
  
  // Calculate the potential energy
  potentialEnergy = p * TrialWalker->walkerData.potentialEnergy +
                    q * OriginalWalker->walkerData.potentialEnergy;
  
  // Calculate the ne and ee potential energy
  neEnergy = p * TrialWalker->walkerData.neEnergy + 
             q * OriginalWalker->walkerData.neEnergy;
  eeEnergy = p * TrialWalker->walkerData.eeEnergy + 
             q * OriginalWalker->walkerData.eeEnergy;
  
  // Calculate the DistanceMovedAccepted this is the average distance
  // moved on a step
  distanceMovedAccepted = p * dR2;

  if (Input->flags.calculate_bf_density == 1)
    {
      for (int i=0; i<Input->WF.getNumberBasisFunctions(); i++)
        walkerData.chiDensity(i) =
          p * TrialWalker->walkerData.chiDensity(i) +
          q * OriginalWalker->walkerData.chiDensity(i);
    }
  
  if(Input->flags.nuclear_derivatives != "none")
    {
      for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
        for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
          walkerData.nuclearDerivatives(d1,d2) =
            p * TrialWalker->walkerData.nuclearDerivatives(d1,d2) +
            q * OriginalWalker->walkerData.nuclearDerivatives(d1,d2);
    }
  
  //eventually, i'll remove the r2 and r12 measurement. but for Helium, it's
  //nice to use when verifying FW
  double calcR12_T=0, calcR12_O=0;
  double calcR1_T=0, calcR1_O=0;
  double calcR2_T=0, calcR2_O=0;
  r2 = 0.0;
  for(int i=0; i<R.dim2(); i++)
    {
      //assuming the nucleus is at the origin
      //also, we're measuring the same thing for both electrons
      //and averaging
      calcR1_T += TrialWalker->R(0,i)*TrialWalker->R(0,i);
      calcR2_T += TrialWalker->R(1,i)*TrialWalker->R(1,i);
      calcR1_O += OriginalWalker->R(0,i)*OriginalWalker->R(0,i);
      calcR2_O += OriginalWalker->R(1,i)*OriginalWalker->R(1,i);
      
      double tempT, tempO;
      tempT = TrialWalker->R(0,i)    - TrialWalker->R(1,i);
      tempO = OriginalWalker->R(0,i) - OriginalWalker->R(1,i);
      
      calcR12_T += tempT*tempT;
      calcR12_O += tempO*tempO;
    }
  r12  = p*sqrt(calcR12_T) + q*sqrt(calcR12_O);
  ir12 = p/sqrt(calcR12_T) + q/sqrt(calcR12_O);
  ir = (p*(1.0/sqrt(calcR1_T) + 1.0/sqrt(calcR2_T))
     +  q*(1.0/sqrt(calcR1_O) + 1.0/sqrt(calcR2_O)))/2.0;
  r2 = (p*(calcR1_T + calcR2_T)
     +  q*(calcR1_O + calcR2_O))/2.0;

  /*
    This section of code probably measure anything useful
  */
  double kineticEnergy_grad_O = 0.0;
  double kineticEnergy_grad_T = 0.0;
  kineticEnergy_grad_O = OriginalWalker->walkerData.gradPsiRatio.dotAllElectrons(OriginalWalker->walkerData.gradPsiRatio);
  kineticEnergy_grad_T = TrialWalker->walkerData.gradPsiRatio.dotAllElectrons(TrialWalker->walkerData.gradPsiRatio);
  kineticEnergy_grad = p*kineticEnergy_grad_T*TrialWalker->walkerData.psi
                     + q*kineticEnergy_grad_O*OriginalWalker->walkerData.psi;

  //This is the forward walking portion of the calculation as described in:
  //J. Casulleras and J. Boronat, Phys. Rev. B 52, 3654 (1995) aka "CB95"
  //this is what CB95 suggest
  //double aWeight = getWeight();
  //double pWeight = getWeight()/OriginalWalker->getWeight();
  
  //this produces the same results as CB95; the difference cancels off in the end
  double aWeight = 1.0;
  double pWeight = 1.0;
    
  //This happens in both Formula 15 and 16 from CB95
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
    {
      for(int j=0; j<isCollectingFWResults.dim2(); j++)
	{
	  
	  if(Input->flags.future_walking[i] == 0)
	    continue;
	  
	  if( isCollectingFWResults(i,j) != DONE)
	    {
	      fwNormalization(i,j)     *= pWeight;
	      fwR12(i,j)               *= pWeight;
	      fwR2(i,j)                *= pWeight;
	      fwiR12(i,j)              *= pWeight;
	      fwiR(i,j)                *= pWeight;
	      fwEnergy(i,j)            *= pWeight;
	      fwKineticEnergy(i,j)     *= pWeight;
	      fwKineticEnergy_grad(i,j)*= pWeight;
	      fwPotentialEnergy(i,j)   *= pWeight;

	      /*
		//I probably commented out this section because
		//It doesn't work with harmonic oscillators...
	      if(Input->flags.nuclear_derivatives != "none")
		{
		  for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
		    for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
		      (fwNuclearForces(d1,d2))(i,j) *= pWeight;
		}
	      */
	    }
	  
	  if(isCollectingFWResults(i,j) == ACCUM )
	    {
	      //We are collecting results (Formula 15)
	      fwNormalization(i,j)   += aWeight;
	      fwR12(i,j)             += aWeight * r12;
	      fwR2(i,j)              += aWeight * r2;
	      fwiR12(i,j)            += aWeight * ir12;
	      fwiR(i,j)              += aWeight * ir;
	      fwEnergy(i,j)          += aWeight * localEnergy;
	      fwKineticEnergy(i,j)   += aWeight * kineticEnergy;
	      fwKineticEnergy_grad(i,j) += aWeight * kineticEnergy_grad;
	      fwPotentialEnergy(i,j) += aWeight * potentialEnergy;

	      /*
	      if(Input->flags.nuclear_derivatives != "none")
		{
		  for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
		    for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
		      (fwNuclearForces(d1,d2))(i,j) += aWeight * walkerData.nuclearDerivatives(d1,d2);
		}
	      */
	    }
	}
    }
}

void QMCWalker::calculateObservables( QMCProperties & props )
{
  if(getWeight() <= 0.0)
    return;

  // Add the data from this walker to the accumulating properties
  
  /*
    What we want to measure is how bad the bad ones are.
    So, we want to set this cutoff high enough to cut out
    all the good to decent walkers.
   */
  /*
  if(ageMoved > 1){
    props.walkerAge.newSample(ageMoved,1.0);
  }
  */
  ageMoved = -1;


  props.weightChange.newSample(dW,1.0);
    
  double rel_diff = fabs( (localEnergy -
			   Input->flags.energy_estimated_original)/
			  Input->flags.energy_estimated_original);
  if(rel_diff < Input->flags.rel_cutoff)
    {
      // Calculate the Energy ...
      props.energy.newSample( localEnergy, getWeight() );
      
      // Calculate the Kinetic Energy  
      props.kineticEnergy.newSample( kineticEnergy, getWeight() );
      
      // Calculate the Potential Energy  
      props.potentialEnergy.newSample( potentialEnergy, getWeight() );
      
      // Calculate the ne and ee potential energy
      props.neEnergy.newSample( neEnergy, getWeight() );
      props.eeEnergy.newSample( eeEnergy, getWeight() );
    }
  else
    {
      if(age == 0)
	{
	  bool shouldWarn = false;
	  switch(Input->flags.warn_verbosity)
	    {
	    case 0: break;
	    case 1: if(rel_diff > 5.0) shouldWarn = true; break;
	    case 2: if(rel_diff > 2.0) shouldWarn = true; break;
	    case 3: if(rel_diff > 1.0) shouldWarn = true; break;
	    case 4: shouldWarn = true; break;
	    default: break;
	    }

	  if(shouldWarn && iteration + Input->flags.equilibration_steps > 5){
	    cerr << "ERROR: Not including walker in average " << ID();
	    cerr << "   localEnergy = " << localEnergy << endl;
	    cerr.flush();
	  }
	}
    }

  // Calculate the Acceptance Probability ...
  props.acceptanceProbability.newSample( getAcceptanceProbability(), getWeight() );
  
  // Calculate the DistanceMovedAccepted this is the average distance
  // moved on a step
  props.distanceMovedAccepted.newSample( distanceMovedAccepted, getWeight() );
  
  // Calculate the DistanceMovedTrial this is the average step length for
  // a trial move
  props.distanceMovedTrial.newSample( dR2, getWeight() );
  
  // Calculate the log of the weight
  props.logWeights.newSample(getWeight(),1);
  
  // Calculate the basis function densities
  if (Input->flags.calculate_bf_density == 1)
    for (int i=0; i<Input->WF.getNumberBasisFunctions(); i++)
      props.chiDensity(i).newSample( walkerData.chiDensity(i), getWeight() );
}

void QMCWalker::calculateObservables( QMCFutureWalkingProperties & fwProps )
{
  double fWeight = getWeight();

  // Add the data from this walker to the accumulating properties  
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
    {
      numFWSteps[i]++;
      
      //This collects the non-forward walking results
      if(Input->flags.future_walking[i] == 0)
	{
	  // Calculate the nuclear forces      
	  if(Input->flags.nuclear_derivatives != "none")
	    for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
	      for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
		(fwProps.nuclearForces(i))(d1,d2).newSample( walkerData.nuclearDerivatives(d1,d2), fWeight );
	  
	  (fwProps.props[FW_It])(i).newSample(1.0, fWeight);
	  (fwProps.props[FW_R12])(i).newSample(r12, fWeight);
	  (fwProps.props[FW_R2])(i).newSample(r2, fWeight);
	  (fwProps.props[FW_iR12])(i).newSample(ir12, fWeight);
	  (fwProps.props[FW_iR])(i).newSample(ir, fWeight);
	  (fwProps.props[FW_TE])(i).newSample(localEnergy, fWeight);
	  (fwProps.props[FW_KE])(i).newSample(kineticEnergy, fWeight);
	  (fwProps.props[FW_KEg])(i).newSample(kineticEnergy_grad, fWeight);
	  (fwProps.props[FW_PE])(i).newSample(potentialEnergy, fWeight);

	  (fwProps.props[FW_R12_2])(i).newSample(r12*r12, fWeight);
	  (fwProps.props[FW_R2_2])(i).newSample(r2*r2, fWeight);
	  (fwProps.props[FW_iR12_2])(i).newSample(ir12*ir12, fWeight);
	  (fwProps.props[FW_iR_2])(i).newSample(ir*ir, fWeight);
	  (fwProps.props[FW_TE_2])(i).newSample(localEnergy*localEnergy, fWeight);
	  (fwProps.props[FW_KE_2])(i).newSample(kineticEnergy*kineticEnergy, fWeight);
	  (fwProps.props[FW_KEg_2])(i).newSample(kineticEnergy_grad*kineticEnergy_grad, fWeight);
	  (fwProps.props[FW_PE_2])(i).newSample(potentialEnergy*potentialEnergy, fWeight);
	  continue;
	}
      
      //this should be >=, but with > we can allow 1.0 for maxFWAsymp
      //so we actually do 1 more step than was input
      if(numFWSteps[i] > Input->flags.future_walking[i])
	{
	  int whichIsDone = isCollectingFWResults(i,0) == DONE ? 0 : 1;
	  int otherStage = (whichIsDone+1)%2;
	  isCollectingFWResults(i,whichIsDone) = ACCUM;
	  isCollectingFWResults(i,otherStage)  = ASYMP;
	  numFWSteps[i] = 0;  
	}
      else if(numFWSteps[i] == int(maxFWAsymp*Input->flags.future_walking[i]+0.5))
	{
	  //CB95 do not suggest how long to propagate the FW steps into the next
	  //block, so maxFWAsymp is the percentage of the next block to use.
	  //I've found that a value of maxFWAsymp=1 works best.

	  int whichIsDone = isCollectingFWResults(i,0) == ACCUM ? 1 : 0;
	  
	  if(isCollectingFWResults(i,whichIsDone) == DONE)
	    break;
	  isCollectingFWResults(i,whichIsDone) = DONE;
	  
	  double norm = 1.0/fwNormalization(i,whichIsDone);
	  (fwProps.props[FW_It])(i).newSample(fwNormalization(i,whichIsDone)/numFWSteps[i],fWeight);
	  (fwProps.props[FW_R12])(i).newSample(fwR12(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_R2])(i).newSample(fwR2(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_iR12])(i).newSample(fwiR12(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_iR])(i).newSample(fwiR(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_TE])(i).newSample(fwEnergy(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_KE])(i).newSample(fwKineticEnergy(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_KEg])(i).newSample(fwKineticEnergy_grad(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_PE])(i).newSample(fwPotentialEnergy(i,whichIsDone)*norm,fWeight);

	  (fwProps.props[FW_R12_2])(i).newSample(fwR12(i,whichIsDone)*norm*r12,fWeight);
	  (fwProps.props[FW_R2_2])(i).newSample(fwR2(i,whichIsDone)*norm*r2,fWeight);
	  (fwProps.props[FW_iR12_2])(i).newSample(fwiR12(i,whichIsDone)*norm*ir12,fWeight);
	  (fwProps.props[FW_iR_2])(i).newSample(fwiR(i,whichIsDone)*norm*ir,fWeight);
	  (fwProps.props[FW_TE_2])(i).newSample(fwEnergy(i,whichIsDone)*norm*localEnergy,fWeight);
	  (fwProps.props[FW_KE_2])(i).newSample(fwKineticEnergy(i,whichIsDone)*norm*kineticEnergy,fWeight);
	  (fwProps.props[FW_KEg_2])(i).newSample(fwKineticEnergy_grad(i,whichIsDone)*norm*kineticEnergy,fWeight);
	  (fwProps.props[FW_PE_2])(i).newSample(fwPotentialEnergy(i,whichIsDone)*norm*potentialEnergy,fWeight);
	  
	  // Calculate the nuclear forces      
	  if(Input->flags.nuclear_derivatives != "none")//what about binning?
	    for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
	      for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
		(fwProps.nuclearForces(i))(d1,d2).newSample( (fwNuclearForces(d1,d2))(i,whichIsDone)*norm, fWeight );
	  
	  resetFutureWalking(i,whichIsDone);
	} 
    }
}

void QMCWalker::resetFutureWalking(int whichBlock, int whichIsDone)
{
  if(whichIsDone > 1 || whichIsDone < 0)
  {
    cerr << "Error: There are only 2 FW stages!\n";
    exit(1);
  }
  
  fwNormalization(whichBlock,whichIsDone)     = 0;
  fwR12(whichBlock,whichIsDone)               = 0;
  fwR2(whichBlock,whichIsDone)                = 0;
  fwiR12(whichBlock,whichIsDone)              = 0;
  fwiR(whichBlock,whichIsDone)                = 0;
  fwEnergy(whichBlock,whichIsDone)            = 0;
  fwKineticEnergy(whichBlock,whichIsDone)     = 0;
  fwKineticEnergy_grad(whichBlock,whichIsDone)     = 0;
  fwPotentialEnergy(whichBlock,whichIsDone)   = 0;
      
  if(Input->flags.nuclear_derivatives != "none")
    for (int d1=0; d1<fwNuclearForces.dim1(); d1++)
      for (int d2=0; d2<fwNuclearForces.dim2(); d2++)
        (fwNuclearForces(d1,d2))(whichBlock,whichIsDone) = 0;
}

void QMCWalker::resetFutureWalking()
{
  numFWSteps = 0;
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
  {
    isCollectingFWResults(i,0) = ACCUM;
    isCollectingFWResults(i,1) = DONE;
    for(int j=0; j<isCollectingFWResults.dim2(); j++)
      resetFutureWalking(i,j);
  }
}


bool QMCWalker::isSingular()
{
  return walkerData.isSingular;
}

double QMCWalker::getLocalEnergyEstimator()
{
  return localEnergy;
}

