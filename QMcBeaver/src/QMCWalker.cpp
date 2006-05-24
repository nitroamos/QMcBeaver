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

QMCWalker::QMCWalker()
{
  TrialWalker    = 0;
  OriginalWalker = 0;
  
  weight = 1.0;
  age    = 0;
  dR2 = 0.0;
}

QMCWalker::QMCWalker( const QMCWalker & rhs )
{
  TrialWalker    = 0;
  OriginalWalker = 0;
  
  *this = rhs;
}

QMCWalker::~QMCWalker()
{
  delete TrialWalker;
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
  fwEnergy.deallocate();
  fwKineticEnergy.deallocate();
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
  
  weight = rhs.weight;
  age    = rhs.age;
  
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
  fwEnergy              = rhs.fwEnergy;
  fwKineticEnergy       = rhs.fwKineticEnergy;
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
                                      Array2D<double> * &rToCalc)
{
  createChildWalkers();
  forwardGreensFunction = TrialWalker->moveElectrons();
  data = & TrialWalker->walkerData;
  rToCalc = & TrialWalker->R;
}

void QMCWalker::processPropagation(QMCFunctions & QMF)
{
  QMCGreensRatioComponent reverseGreensFunction =
    calculateReverseGreensFunction();
  double GreensFunctionRatio =
    reverseGreensFunction/forwardGreensFunction;
    
  if (IeeeMath::isNaN(GreensFunctionRatio))
    calculateMoveAcceptanceProbability(0.0);
  else
  calculateMoveAcceptanceProbability(GreensFunctionRatio);
  
  acceptOrRejectMove();
  reweight_walker();
  calculateObservables();
  
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Reinitializing a singular walker!!" << endl;
      initializeWalkerPosition(QMF);
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
  
  Array2D<double> Displacement(Input->WF.getNumberElectrons(),3);
  
  // Add the drift to the displacement
  
  // Don't use the QF  =>  R' = R + gauss_rn
  Displacement = 0;
  
  // Add the randomness to the displacement
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<Displacement.dim2(); j++)
      Displacement(i,j) += sigma*gasdev(&Input->flags.iseed);
    
  // Calculate the square of the magnitude of the displacement
  dR2 = 0;
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<Displacement.dim2(); j++)
      dR2 += Displacement(i,j) * Displacement(i,j);
    
  // Now update the R
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<Displacement.dim2(); j++)
      R(i,j) += Displacement(i,j);
    
  // Return the greens function for the forward move
  
  QMCGreensRatioComponent GF = QMCGreensRatioComponent(1.0);
  return GF;
}

QMCGreensRatioComponent QMCWalker::moveElectronsImportanceSampling()
{
  // Move the electrons of this walker
  double sigma = sqrt(Input->flags.dt);
  
  Array2D<double> Displacement(Input->WF.getNumberElectrons(),3);
  
  // Add the drift to the displacement
  
  // Use the QF => R' = R + dt * QF + gauss_rn
  Displacement = walkerData.modifiedGradPsiRatio;
  Displacement *= Input->flags.dt;
  
  // Add the randomness to the displacement
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<Displacement.dim2(); j++)
      Displacement(i,j) += sigma*gasdev(&Input->flags.iseed);
    
  // Calculate the square of the magnitude of the displacement
  dR2 = 0.0;
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<Displacement.dim2(); j++)
      dR2 += Displacement(i,j) * Displacement(i,j);
    
  // Calculate the Green's function for the forward move
  
  double greens = 0.0;
  double tau = Input->flags.dt;
  
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<3; j++)
      {
	double temp = Displacement(i,j)-tau*Displacement(i,j);
	greens += temp*temp;
    }
    
  greens = greens/(2*tau);
  
  double k = 1.0;
  double a = 2.0*3.14159265359*tau;
  double b = -1.5*Displacement.dim1();
  double c = -greens;
  
  QMCGreensRatioComponent GF(k,a,b,c);
  
  // Now update the R
  for(int i=0; i<Displacement.dim1(); i++)
    for(int j=0; j<Displacement.dim2(); j++)
      R(i,j) += Displacement(i,j);
    
  return GF;
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
      
      if( probabilityGaussianTypeMove > ran1(&Input->flags.iseed) )
        {
          // Gaussian Type Move
          
          // Particle is moved in the direction
          // Rnuc + radialCoordinate * radialUnitVector +
          //   zCoordinate * zUnitVector +
          //   gaussian random number with standard deviation sqrt(tau)
          for(int i=0; i<3; i++)
              newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i)
      + radialCoordinate * radialUnitVector(i) + zCoordinate * zUnitVector(i) +
                               sqrt(tau)*gasdev(&Input->flags.iseed);
            }
      else
        {
          // Slater Type Move
          
          for(int i=0; i<3; i++)
	    newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i);
            
          // add random part
          
          double r = randomDistribution1(&Input->flags.iseed)/(2.0*expParam);
          double phi = 2*3.14159265359*ran1(&Input->flags.iseed);
          double theta = sindev(&Input->flags.iseed);
          
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
  QMCGreensRatioComponent GF = QMCGreensRatioComponent(1.0);
  return GF;
}

QMCGreensRatioComponent \
QMCWalker::calculateReverseGreensFunctionImportanceSampling()
{
  double tau = Input->flags.dt;
  
  double greens = 0.0;
  
  for(int i=0; i< R.dim1(); i++)
      for(int j=0; j<3; j++)
        {
          double temp = OriginalWalker->R(i,j) - TrialWalker->R(i,j) -
                        tau*TrialWalker->walkerData.modifiedGradPsiRatio(i,j);
                        
          greens += temp*temp;
        }
    
  greens = greens/(2*tau);
  
  greens = pow(2*3.14159265359*tau,-1.5*R.dim1()) * exp(-greens);
  
  double k = 1.0;
  double a = 2.0*3.14159265359*tau;
  double b = -1.5*R.dim1();
  double c = -greens;
  QMCGreensRatioComponent GF = QMCGreensRatioComponent(k,a,b,c);
  return GF;
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
  double dW = 0.0;
  // determine the weighting factor dW so that the new weight = weight*dW
  
  bool weightIsNaN = false;

  if( Input->flags.run_type == "variational" )
      // Keep weights constant for VMC
      dW = 1.0;

  else
    {
      double S_trial = 0.0;
      double S_original = 0.0;
      
      double trialEnergy = TrialWalker->walkerData.localEnergy;
      double originalEnergy = OriginalWalker->walkerData.localEnergy;
      
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
            
	  if (IeeeMath::isNaN(lengthGradTrialModified) || IeeeMath::isNaN(lengthGradTrialUnmodified))
	    {
	      cerr << "WARNING: modified Grad Psi Ratio is NaN in "; 
	      cerr << "QMCWalker::reweight_walker()" << endl;
	      weightIsNaN = true;
	    }
	  if (IeeeMath::isNaN(lengthGradOriginalModified) || IeeeMath::isNaN(lengthGradOriginalUnmodified))
	    {
	      cerr << "WARNING: original Grad Psi Ratio is NaN in "; 
	      cerr << "QMCWalker::reweight_walker()" << endl;
	      weightIsNaN = true;
	    }
	  else if(Input->flags.energy_modification_type=="modified_umrigar93")
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
          double q = OriginalWalker->getAcceptanceProbability();
          
	  double temp = (p*0.5*(S_original+S_trial) + q*S_original) *
                        Input->flags.dt_effective;
                        
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
      else
	{
	  cerr << "ERROR: unknown reweighting method!" << endl;
          exit(0);
        }
    }
    
#ifdef DELETE_LARGE_WEIGHT_WALKERS
    
  // If a very erratic point is discovered during reweighting,
  // give it a statistical weight of zero so that it will be not
  // mess up the calculation.  The cutoff here needs to be played with
  // to get the best value.
  
  if( dW > 10 )
    {
      cerr << "WARNING: Walkers weight is being increased by a factor of "
      << dW << "!" << endl;
      cerr << "Walker's weight is being set to zero so that it does not"
      << " ruin the whole calculation or use all of the "
      << "machine's memory." << endl;
      dW = 0.0;
    }
#endif
    
  // now set the weight of the walker
  setWeight( getWeight() * dW );
}

void QMCWalker::calculateMoveAcceptanceProbability(double GreensRatio)
{
  // This tells us the probability of accepting or rejecting a proposed move
  
  double PsiRatio = TrialWalker->walkerData.psi/OriginalWalker->walkerData.psi;
  
  // increase the probability of accepting a move if the walker has not
  // moved in a long time
  
  double MoveOldWalkerFactor = 1.0;
  
  if( age > Input->flags.old_walker_acceptance_parameter )
    MoveOldWalkerFactor =
      pow(1.1,age-Input->flags.old_walker_acceptance_parameter);
      
  // calculate the probability of accepting the trial move
  double p = PsiRatio * PsiRatio * GreensRatio * MoveOldWalkerFactor;
  
  if( !(IeeeMath::isNaN(p)) && age > 2*Input->flags.old_walker_acceptance_parameter )
    p = 1;
    
  // if the aratio is NaN then reject the move
  if( IeeeMath::isNaN(p) )
    {
      cerr << "WARNING: Rejecting trial walker with NaN aratio!" << endl;
      p = 0.0;
      // The energies are assigned zero because the later multiplication by 
      // p=0 doesn't result in 0 contribution.
      TrialWalker->walkerData.kineticEnergy   = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy     = 0;
      TrialWalker->walkerData.neEnergy        = 0;
      TrialWalker->walkerData.eeEnergy        = 0;
    }
    
  // The particular NaN that this is correcting for is not revealed by isinf or
  // isnan...
  double kineticEnergy = TrialWalker->walkerData.kineticEnergy;
  if( IeeeMath::isNaN(kineticEnergy) )
    {
      cerr << "WARNING: Rejecting trial walker with NaN kinetic energy!" << endl;
      p = 0.0;
      TrialWalker->walkerData.kineticEnergy   = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy     = 0;
      TrialWalker->walkerData.neEnergy        = 0;
      TrialWalker->walkerData.eeEnergy        = 0;
    }

  double potentialEnergy = TrialWalker->walkerData.potentialEnergy;
  if( IeeeMath::isNaN(potentialEnergy) )
    {
      cerr << "WARNING: Rejecting trial walker with NaN potential energy!" << endl;
      p = 0.0;
      TrialWalker->walkerData.kineticEnergy   = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy     = 0;
      TrialWalker->walkerData.neEnergy        = 0;
      TrialWalker->walkerData.eeEnergy        = 0;
    }
    
  // if the trial position is singular reject the move
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Rejecting singular trial walker!" << endl;
      p = 0.0;
      TrialWalker->walkerData.kineticEnergy   = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy     = 0;
      TrialWalker->walkerData.neEnergy        = 0;
      TrialWalker->walkerData.eeEnergy        = 0;
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
  OriginalWalker->setAcceptanceProbability(1.0-p);
}

void QMCWalker::acceptOrRejectMove()
{
  if( TrialWalker->getAcceptanceProbability() > ran1(&Input->flags.iseed) )
    {
      // accept the move
      *this = *TrialWalker;
      age = 0;
      move_accepted = true;
    }
  else
    {
      // reject the move
      *this = *OriginalWalker;
      age++;
      move_accepted = false;
    }
    
  //if( age > Input->flags.old_walker_acceptance_parameter )
  // {
  //commented this out to time more accurately
  //CHIP we will want this back in
  //cerr << "WARNING: Walker older than "
  //<< Input->flags.old_walker_acceptance_parameter << " steps!"
  //    << " age = " << age << endl;
  //}
}

void QMCWalker::createChildWalkers()
{
  if( TrialWalker == 0 )
    {
      TrialWalker = new QMCWalker();
    }
  *TrialWalker = *this;
  
  if( OriginalWalker == 0 )
    {
      OriginalWalker = new QMCWalker();
    }
  *OriginalWalker = *this;
}

void QMCWalker::initialize(QMCInput *INPUT)
{
  Input = INPUT;
  walkerData.gradPsiRatio.allocate(Input->WF.getNumberElectrons(),3);
  walkerData.modifiedGradPsiRatio.allocate(Input->WF.getNumberElectrons(),3);
  walkerData.localEnergy     = 0.0;
  walkerData.kineticEnergy   = 0.0;
  walkerData.potentialEnergy = 0.0;
  walkerData.neEnergy        = 0.0;
  walkerData.eeEnergy        = 0.0;
  walkerData.psi             = 0.0;
  walkerData.configOutput    = new stringstream();
  
  int numFW = Input->flags.future_walking.size();
  numFWSteps.resize(numFW);
  
  isCollectingFWResults.allocate(numFW,2);
  fwNormalization.allocate(numFW,2);
  fwR12.allocate(numFW,2);
  fwR2.allocate(numFW,2);
  fwEnergy.allocate(numFW,2);
  fwKineticEnergy.allocate(numFW,2);
  fwPotentialEnergy.allocate(numFW,2);
  
  if(Input->flags.nuclear_derivatives != "none")
  {    
    if(Input->flags.nuclear_derivatives != "bin_force_density")
    {
      walkerData.nuclearDerivatives.allocate(Input->Molecule.getNumberAtoms(), 3);
    } else {
      walkerData.nuclearDerivatives.allocate(QMCNuclearForces::getNumBins(), 1);      
    }
  
    fwNuclearForces.allocate(walkerData.nuclearDerivatives.dim1(),
                             walkerData.nuclearDerivatives.dim2());

    for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
      for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
        fwNuclearForces(d1,d2).allocate(numFW,2);
  }
      
  resetFutureWalking();
    
  if (Input->flags.calculate_bf_density == 1)
    walkerData.chiDensity.allocate(Input->WF.getNumberBasisFunctions());
    
  //initialize acceptance probability
  AcceptanceProbability = 0.0;
    
  // Make current position in 3N space
  // N by 3 matrix, R[which electron][x,y, or z]
  R.allocate(Input->WF.getNumberElectrons(),3);
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

void QMCWalker::toXML(ostream& strm)
{
  strm << "<QMCWalker>" << endl;
  strm << "\t<Position>" <<endl;
  for(int ep=0; ep<Input->WF.getNumberElectrons(); ep++)
    {
      strm << "\t\t";
      for(int j=0;j<3;j++)
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

void QMCWalker::readXML(istream& strm)
{
  string temp;
  strm >> temp;
  
  // Read position
  strm >> temp;
  
  for(int ep=0; ep<Input->WF.getNumberElectrons(); ep++)
    {
      for(int j=0;j<3;j++)
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
  
  QMCFunctions QMF(Input);
  QMF.evaluate(R,walkerData);
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

void QMCWalker::setR(Array2D<double>& temp_R)
{
  R = temp_R;
}

QMCWalkerData* QMCWalker::getWalkerData()
{
  return &walkerData;
}

void QMCWalker::calculateObservables()
{
  double p = TrialWalker->AcceptanceProbability;
  double q = OriginalWalker->AcceptanceProbability;
  
  if( fabs(p + q - 1.0) > 1e-10 )
    {
      cerr << "ERROR: In QMCWalker::calculateObservables probabilities "
	   << "do not sum correctly!" << endl;
      exit(0);
    }
  
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

  // Calculate the acceptance probability
  AcceptanceProbability = p;
  
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
  r2 = 0.0;
  for(int i=0; i<3; i++)
    {
      r2 += p*TrialWalker->R(0,i)*TrialWalker->R(0,i)
	+  p*TrialWalker->R(1,i)*TrialWalker->R(1,i)
	+  q*OriginalWalker->R(0,i)*OriginalWalker->R(0,i)
	+  q*OriginalWalker->R(1,i)*OriginalWalker->R(1,i);
      
      double tempT, tempO;
      tempT = TrialWalker->R(0,i)    - TrialWalker->R(1,i);
      tempO = OriginalWalker->R(0,i) - OriginalWalker->R(1,i);
      
      calcR12_T += tempT*tempT;
      calcR12_O += tempO*tempO;
    }
  r12 = p*sqrt(calcR12_T) + q*sqrt(calcR12_O);
  r2 /= 2.0;
  
  //This is the forward walking portion of the calculation as described in:
  //J. Casulleras and J. Boronat, Phys. Rev. B 52, 3654 (1995) aka "CB95"
  //this is what CB95 suggest
  //double aWeight = getWeight();
  
  //but this is what seems to be working
  double aWeight = 1.0;
  
  //this is what CB95 suggest
  //double pWeight = getWeight()/OriginalWalker->getWeight();
  
  //although it doesn't make much difference, this seems to work better
  double pWeight = 1.0;
  
  //this might create nan if run with diffusion
  //Don't use for that reason!
  //double pWeight = getWeight();
  
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
	      fwEnergy(i,j)            *= pWeight;
	      fwKineticEnergy(i,j)     *= pWeight;
	      fwPotentialEnergy(i,j)   *= pWeight;
	      
	      for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
		for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
		  (fwNuclearForces(d1,d2))(i,j) *= pWeight;
	    }
	  
	  if(isCollectingFWResults(i,j) == ACCUM )
	    {
	      //We are collecting results (Formula 15)
	      fwNormalization(i,j)   += aWeight;
	      fwR12(i,j)             += aWeight * r12;
	      fwR2(i,j)              += aWeight * r2;
	      fwEnergy(i,j)          += aWeight * localEnergy;
	      fwKineticEnergy(i,j)   += aWeight * kineticEnergy;
	      fwPotentialEnergy(i,j) += aWeight * potentialEnergy;
	      
	      for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
		for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
		  (fwNuclearForces(d1,d2))(i,j) += aWeight * walkerData.nuclearDerivatives(d1,d2);
	    }
	}
    }
}

void QMCWalker::calculateObservables( QMCProperties & props )
{
  // Add the data from this walker to the accumulating properties
  
  // Calculate the Energy ...
  props.energy.newSample( localEnergy, getWeight() );
  
  // Calculate the Kinetic Energy  
  props.kineticEnergy.newSample( kineticEnergy, getWeight() );
  
  // Calculate the Potential Energy  
  props.potentialEnergy.newSample( potentialEnergy, getWeight() );
  
  // Calculate the ne and ee potential energy
  props.neEnergy.newSample( neEnergy, getWeight() );
  props.eeEnergy.newSample( eeEnergy, getWeight() );
  
  // Calculate the Acceptance Probability ...
  props.acceptanceProbability.newSample( AcceptanceProbability, getWeight() );
  
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

void QMCWalker::calculateObservables( QMCFutureWalkingProperties & props )
{
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
		(props.nuclearForces(i))(d1,d2).newSample( walkerData.nuclearDerivatives(d1,d2), getWeight() );
	  
	  props.r12(i).newSample(r12, getWeight());
	  props.r2(i).newSample(r2, getWeight());
	  props.fwEnergy(i).newSample(localEnergy, getWeight());
	  props.fwKineticEnergy(i).newSample(kineticEnergy, getWeight());
	  props.fwPotentialEnergy(i).newSample(potentialEnergy, getWeight());
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
	  props.r12(i).newSample(fwR12(i,whichIsDone)*norm,getWeight());
	  props.r2(i).newSample(fwR2(i,whichIsDone)*norm,getWeight());
	  props.fwEnergy(i).newSample(fwEnergy(i,whichIsDone)*norm,getWeight());
	  props.fwKineticEnergy(i).newSample(fwKineticEnergy(i,whichIsDone)*norm,getWeight());
	  props.fwPotentialEnergy(i).newSample(fwPotentialEnergy(i,whichIsDone)*norm,getWeight());
	  
	  // Calculate the nuclear forces      
	  if(Input->flags.nuclear_derivatives != "none")
	    for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
	      for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
		(props.nuclearForces(i))(d1,d2).newSample( (fwNuclearForces(d1,d2))(i,whichIsDone)*norm, 1.0 );
	  
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
  fwEnergy(whichBlock,whichIsDone)            = 0;
  fwKineticEnergy(whichBlock,whichIsDone)     = 0;
  fwPotentialEnergy(whichBlock,whichIsDone)   = 0;
      
  if(Input->flags.nuclear_derivatives != "none")
    for (int d1=0; d1<fwNuclearForces.dim1(); d1++)
      for (int d2=0; d2<fwNuclearForces.dim2(); d2++)
        (fwNuclearForces(d1,d2))(whichBlock,whichIsDone) = 0;
}

void QMCWalker::resetFutureWalking()
{
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
  {
    isCollectingFWResults(i,0) = ACCUM;
    isCollectingFWResults(i,1) = DONE;
    numFWSteps[i] = 0;
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

