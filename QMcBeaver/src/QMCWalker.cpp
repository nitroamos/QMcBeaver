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

void QMCWalker::processPropagation()
{
  QMCGreensRatioComponent reverseGreensFunction =
    calculateReverseGreensFunction();
  double GreensFunctionRatio =
    reverseGreensFunction/forwardGreensFunction;
    
  calculateMoveAcceptanceProbability(GreensFunctionRatio);
  
  acceptOrRejectMove();
  reweight_walker();
  calculateObservables();
  
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Reinitializing a singular walker!!" << endl;
      initializeWalkerPosition();
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
    {
      return moveElectronsNoImportanceSampling();
    }
  else if(Input->flags.sampling_method == "importance_sampling" )
    {
      return moveElectronsImportanceSampling();
    }
  else if(Input->flags.sampling_method == "umrigar93_importance_sampling")
    {
      return moveElectronsUmrigar93ImportanceSampling();
    }
  else
    {
      cerr << "ERROR: Improper value for sampling_method set ("
      << Input->flags.sampling_method
      << ")!" << endl;
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
    {
      for(int j=0; j<Displacement.dim2(); j++)
        {
          Displacement(i,j) += sigma*gasdev(&Input->flags.iseed);
        }
    }
    
  // Calculate the square of the magnitude of the displacement
  dR2 = 0;
  for(int i=0; i<Displacement.dim1(); i++)
    {
      for(int j=0; j<Displacement.dim2(); j++)
        {
          dR2 += Displacement(i,j) * Displacement(i,j);
        }
    }
    
  // Now update the R
  for(int i=0; i<Displacement.dim1(); i++)
    {
      for(int j=0; j<Displacement.dim2(); j++)
        {
          R(i,j) += Displacement(i,j);
        }
    }
    
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
    {
      for(int j=0; j<Displacement.dim2(); j++)
        {
          Displacement(i,j) += sigma*gasdev(&Input->flags.iseed);
        }
    }
    
  // Calculate the square of the magnitude of the displacement
  dR2 = 0.0;
  for(int i=0; i<Displacement.dim1(); i++)
    {
      for(int j=0; j<Displacement.dim2(); j++)
        {
          dR2 += Displacement(i,j) * Displacement(i,j);
        }
    }
    
  // Calculate the Green's function for the forward move
  
  double greens = 0.0;
  double tau = Input->flags.dt;
  
  for(int i=0; i<Displacement.dim1(); i++)
    {
      for(int j=0; j<3; j++)
        {
          double temp = Displacement(i,j)-
                        tau*Displacement(i,j);
                        
          greens += temp*temp;
        }
    }
    
  greens = greens/(2*tau);
  
  double k = 1.0;
  double a = 2.0*3.14159265359*tau;
  double b = -1.5*Displacement.dim1();
  double c = -greens;
  
  QMCGreensRatioComponent GF(k,a,b,c);
  
  // Now update the R
  for(int i=0; i<Displacement.dim1(); i++)
    {
      for(int j=0; j<Displacement.dim2(); j++)
        {
          R(i,j) += Displacement(i,j);
        }
    }
    
  return GF;
}

QMCGreensRatioComponent QMCWalker::moveElectronsUmrigar93ImportanceSampling()
{
  double tau = Input->flags.dt;
  QMCGreensRatioComponent GF(1.0);
  dR2 = 0.0;
  Array2D<double> & Displacement = walkerData.modifiedGradPsiRatio;
  
  Array1D<double> zUnitVector(3);
  Array1D<double> radialUnitVector(3);
  Array1D<double> newPosition(3);
  
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
        {
          zComponentQF += zUnitVector(i)*Displacement(electron,i);
        }
        
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
      double zCoordinate = max(distanceFromNucleus + zComponentQF * tau,
                               0.0);
                               
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
            {
              newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i)
                               + radialCoordinate * radialUnitVector(i) +
                               zCoordinate * zUnitVector(i) +
                               sqrt(tau)*gasdev(&Input->flags.iseed);
            }
        }
      else
        {
          // Slater Type Move
          
          for(int i=0; i<3; i++)
            {
              newPosition(i) =
                Input->Molecule.Atom_Positions(nearestNucleus,i);
            }
            
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
      
      for(int i=0; i<3; i++)
        {
          double temp = newPosition(i) -
                        ( Input->Molecule.Atom_Positions(nearestNucleus,i) +
                          radialCoordinate * radialUnitVector(i) +
                          zCoordinate * zUnitVector(i) );
                          
          distance1Sq += temp*temp;
        }
        
      double ga = 2*3.14159265359*tau;
      double gb = -1.5;
      double gc = -distance1Sq/(2*tau);
      
      QMCGreensRatioComponent GaussianGF(probabilityGaussianTypeMove,ga,gb,gc);
      
      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          double temp = newPosition(i) -
                        Input->Molecule.Atom_Positions(nearestNucleus,i);
                        
          distance2Sq += temp*temp;
        }
        
      double sk = probabilitySlaterTypeMove/3.14159265359;
      double sa = expParam;
      double sb = 3.0;
      double sc = -2.0*expParam*sqrt(distance2Sq);
      
      QMCGreensRatioComponent SlaterGF(sk,sa,sb,sc);
      
      QMCGreensRatioComponent OneE = GaussianGF + SlaterGF;
      
      GF *= OneE;
      
      // Update the distance attempted to move squared
      
      for(int i=0; i<3; i++)
        {
          double temp = newPosition(i) - R(electron,i);
          dR2 += temp*temp;
        }
        
      // Move the electron
      
      for(int i=0; i<3; i++)
        {
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
    {
      return calculateReverseGreensFunctionNoImportanceSampling();
    }
  else if(Input->flags.sampling_method == "importance_sampling" )
    {
      return calculateReverseGreensFunctionImportanceSampling();
    }
  else if(Input->flags.sampling_method == "umrigar93_importance_sampling")
    {
      return calculateReverseGreensFunctionUmrigar93ImportanceSampling();
    }
  else
    {
      cerr << "ERROR: Improper value for sampling_method set ("
      << Input->flags.sampling_method
      << ")!" << endl;
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
    {
      for(int j=0; j<3; j++)
        {
          double temp = OriginalWalker->R(i,j) - TrialWalker->R(i,j) -
                        tau*TrialWalker->walkerData.modifiedGradPsiRatio(i,j);
                        
          greens += temp*temp;
        }
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
  QMCGreensRatioComponent GF(1.0);
  
  Array1D<double> zUnitVector(3);
  Array1D<double> radialUnitVector(3);
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
        {
          zComponentQF += zUnitVector(i)*
                          TrialWalker->walkerData.modifiedGradPsiRatio(electron,i);
        }
        
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
      double zCoordinate = max(distanceFromNucleus + zComponentQF * tau,
                               0.0);
                               
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
      double probabilityGaussianTypeMove = 1.0 - probabilitySlaterTypeMove;
      
      
      // Update the greens function
      
      double distance1Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          double temp = OriginalWalker->R(electron,i) -
                        ( Input->Molecule.Atom_Positions(nearestNucleus,i) +
                          radialCoordinate * radialUnitVector(i) +
                          zCoordinate * zUnitVector(i) );
                          
          distance1Sq += temp*temp;
        }
        
      double ga = 2*3.14159265359*tau;
      double gb = -1.5;
      double gc = -distance1Sq/(2*tau);
      
      QMCGreensRatioComponent GaussianGF(probabilityGaussianTypeMove,ga,gb,gc);
      
      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          double temp = OriginalWalker->R(electron,i) - \
                        Input->Molecule.Atom_Positions(nearestNucleus,i);
                        
          distance2Sq += temp*temp;
        }
        
      double sk = probabilitySlaterTypeMove/3.14159265359;
      double sa = expParam;
      double sb = 3.0;
      double sc = -2.0*expParam*sqrt(distance2Sq);
      
      QMCGreensRatioComponent SlaterGF(sk,sa,sb,sc);
      
      QMCGreensRatioComponent OneE = GaussianGF + SlaterGF;
      
      GF *= OneE;
    }
  return GF;
}

void QMCWalker::reweight_walker()
{
  double dW = 0.0;
  // determine the weighting factor dW so that the new weight = weight*dW
  
  if( Input->flags.run_type == "variational" )
    {
      // Keep weights constant for VMC
      dW = 1.0;
    }
    
  else
    {
      double S_trial;
      double S_original;
      
      double trialEnergy = TrialWalker->walkerData.localEnergy;
      double originalEnergy = OriginalWalker->walkerData.localEnergy;
      
      if( Input->flags.energy_modification_type == "none" )
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
            
          if( Input->flags.energy_modification_type == "modified_umrigar93" )
            {
              S_trial = (Input->flags.energy_trial - trialEnergy)*
                        (lengthGradTrialModified/lengthGradTrialUnmodified);
                        
              S_original = (Input->flags.energy_trial - originalEnergy) *
                           (lengthGradOriginalModified/lengthGradOriginalUnmodified);
            }
          else if( Input->flags.energy_modification_type == "umrigar93" )
            {
              S_trial =
                (Input->flags.energy_trial - Input->flags.energy_estimated) +
                (Input->flags.energy_estimated - trialEnergy)*
                (lengthGradTrialModified/lengthGradTrialUnmodified);
                
              S_original =
                (Input->flags.energy_trial - Input->flags.energy_estimated) +
                (Input->flags.energy_estimated - originalEnergy) *
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
          
          double temp = (p/2.0*(S_original+S_trial) + q*S_original) *
                        Input->flags.dt_effective;
                        
          dW = exp(temp);
        }
      else
        {
          cerr << "ERROR: unknown branching method!" << endl;
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
  
  // Force the walker to move if it is too old
  if( age > 2*Input->flags.old_walker_acceptance_parameter )
    p = 1;
    
  // if the aratio is NaN then reject the move
  if( IeeeMath::isNaN(p) != 0 )
    {
      cerr << "WARNING: Rejecting trial walker with NaN aratio!" << endl;
      p = 0.0;
      //the energies are assigned zero because the later multiplication by p = 0
      //doesn't result in 0 contribution
      TrialWalker->walkerData.kineticEnergy = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy = 0;
    }
    
  //the particular NaN that this is correcting for is not revealed by isinf or isnan...
  double energy = TrialWalker->walkerData.kineticEnergy;
  if( energy++ == energy)
    {
      cerr << "WARNING: Rejecting trial walker with NaN kinetic energy!" << endl;
      p = 0.0;
      TrialWalker->walkerData.kineticEnergy = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy = 0;
    }
    
  // if the trial position is singular reject the move
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Rejecting singular trial walker!" << endl;
      p = 0.0;
      TrialWalker->walkerData.kineticEnergy = 0;
      TrialWalker->walkerData.potentialEnergy = 0;
      TrialWalker->walkerData.localEnergy = 0;
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
  walkerData.psi             = 0.0;
  walkerData.configOutput    = new stringstream();
  if (Input->flags.calculate_bf_density == 1)
    walkerData.chiDensity.allocate(Input->WF.getNumberBasisFunctions());
    
  //initialize acceptance probability
  AcceptanceProbability = 0.0;
  
  // Make current position in 3N space
  // N by 3 matrix, R[which electron][x,y, or z]
  R.allocate(Input->WF.getNumberElectrons(),3);
}

void QMCWalker::writeCorrelatedSamplingConfiguration(ostream& strm)
{
  strm << "&" << endl;
  int R_moved=R.dim1();
  
  strm << "R\t" << R_moved << endl;
  for(int i=0;i<R_moved;i++)
    {
      //now we are printing out all the electrons
      
      strm << "\t" << i << "\t";
      strm << R(i,0) << "\t";
      strm << R(i,1) << "\t";
      strm << R(i,2) << endl;
    }
    
  strm << (*walkerData.configOutput).str();
}

void QMCWalker::calculatePairDistances(double max_pair_distance, double dr,
     Array1D<double> &pll_spin, Array1D<double> &opp_spin, double &totalWeight)
{
  int nalpha = Input->WF.getNumberAlphaElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();
  
  double dist = 0.0;
  int index = 0;

  totalWeight += weight;
  
  if (nalpha > 1)
    {
      for (int i=0; i<nalpha-1; i++)
        for (int j=i+1; j<nalpha; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
            index = int(dist/dr);
            pll_spin(index) += weight;
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
            index = int(dist/dr);
            pll_spin(index) += weight;
          }
    }
    
  if (nalpha > 0 && nbeta > 0)
    {
      for (int i=0; i<nalpha; i++)
        for (int j=nalpha; j<nalpha+nbeta; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
            index = int(dist/dr);
            opp_spin(index) += weight;
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
  
  QMCFunctions QMF;
  QMF.initialize(Input);
  QMF.evaluate(R,walkerData);
}

void QMCWalker::initializeWalkerPosition()
{
  QMCInitializeWalker * IW =
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
        Input->flags.walker_initialization_method);
  QMCFunctions QMF;
  QMF.initialize(Input);
  
  R = IW->initializeWalkerPosition();
  QMF.evaluate(R,walkerData);
  
  int initilization_try = 1;
  while( isSingular() )
    {
      cerr << "Regenerating Walker..." << endl;
      
      if( initilization_try > 100 )
        {
          cerr << "ERROR: 100 consecutive singular configurations while "
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

bool QMCWalker::isSingular()
{
  return walkerData.isSingular;
}

double QMCWalker::getLocalEnergyEstimator()
{
  return localEnergy;
}

