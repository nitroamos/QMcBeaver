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
  // In the equality operator DON'T copy the 
  // pointers to the child walkers.  If you do,
  // then there is the possibility of big memory leak
  // problems.  This is not required because these are 
  // just temporary variables used in propagating the 
  // electrons.

  weight = rhs.weight;
  age    = rhs.age;

  QMF = rhs.QMF;
  Input = rhs.Input;

  move_accepted         = rhs.move_accepted;
 
  AcceptanceProbability = rhs.AcceptanceProbability;
  localEnergy           = rhs.localEnergy;
  potentialEnergy       = rhs.potentialEnergy;
  kineticEnergy         = rhs.kineticEnergy;
  distanceMovedAccepted = rhs.distanceMovedAccepted;
  dR2                   = rhs.dR2;

  R = rhs.R;
}

void QMCWalker::propagateWalker()
{
  createChildWalkers();
  double forwardGreensFunction = TrialWalker->moveElectrons();
  TrialWalker->evaluate();
  double reverseGreensFunction = calculateReverseGreensFunction();
  calculateMoveAcceptanceProbability(forwardGreensFunction,
				     reverseGreensFunction);
  acceptOrRejectMove();
  reweight_walker();
  calculateObservables();

  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Reinitializing a singular walker!!" << endl;
      initializeWalkerPosition();
    }
}

double QMCWalker::moveElectrons()
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
  return 0.0;
}

double QMCWalker::moveElectronsNoImportanceSampling()
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
  return 1.0;
}

double QMCWalker::moveElectronsImportanceSampling()
{
  // Move the electrons of this walker
  double sigma = sqrt(Input->flags.dt);
  
  Array2D<double> Displacement(Input->WF.getNumberElectrons(),3); 

  // Add the drift to the displacement

  // Use the QF => R' = R + dt * QF + gauss_rn
  Displacement = *QMF.getModifiedGradPsiRatio();
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
  dR2 = 0;
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
	    tau*(*QMF.getModifiedGradPsiRatio())(i,j);

	  greens += temp*temp;
	}
    }

  greens = greens/(2*tau);

  greens = pow(2.0*3.14159265359*tau,-1.5*Displacement.dim1()) * exp(-greens);

  // Now update the R
  for(int i=0; i<Displacement.dim1(); i++)
    {
      for(int j=0; j<Displacement.dim2(); j++)
	{
	  R(i,j) += Displacement(i,j);
	}
    }

  return greens;
}

double QMCWalker::moveElectronsUmrigar93ImportanceSampling()
{
  double tau = Input->flags.dt;
  double greens = 1.0;
  dR2 = 0.0;

  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {
      // Find the nearest nucleus to this electron
      int nearestNucleus = Input->Molecule.findClosestNucleusIndex(R,electron);

      // Calculate the unit vector in the Z direction
      Array1D<double> zUnitVector(3);

      double distanceFromNucleus = 0.0;

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
	  zComponentQF += zUnitVector(i)*
	    (*QMF.getModifiedGradPsiRatio())(electron,i);
	}

      Array1D<double> radialUnitVector(3);
      double radialComponentQF = 0.0;

      for(int i=0; i<3; i++)
	{
	  radialUnitVector(i) = (*QMF.getModifiedGradPsiRatio())(electron,i) - 
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

      Array1D<double> newPosition(3);

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

      double g1 = pow(2*3.14159265359*tau,-1.5)*exp( -distance1Sq/(2*tau) );

      double distance2Sq = 0.0;

      for(int i=0; i<3; i++)
	{
	  double temp = newPosition(i) - 
	    Input->Molecule.Atom_Positions(nearestNucleus,i);

	  distance2Sq += temp*temp;
	}

      double g2 = pow(expParam,3.0) / 3.14159265359 *
	exp( -2.0 * expParam * sqrt( distance2Sq ) );

      greens *= (probabilityGaussianTypeMove * g1 + 
		 probabilitySlaterTypeMove * g2 );


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

  return greens;
}

double QMCWalker::calculateReverseGreensFunction()
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
  return 0.0;
}

double QMCWalker::calculateReverseGreensFunctionNoImportanceSampling()
{
  return 1.0;
}

double QMCWalker::calculateReverseGreensFunctionImportanceSampling()
{
  double tau = Input->flags.dt;

  double greens = 0.0;

  for(int i=0; i< R.dim1(); i++)
    {
      for(int j=0; j<3; j++)
	{
	  double temp = OriginalWalker->R(i,j) - TrialWalker->R(i,j) - 
	    tau*(*TrialWalker->QMF.getModifiedGradPsiRatio())(i,j);

	  greens += temp*temp;
	}
    }

  greens = greens/(2*tau);

  greens = pow(2*3.14159265359*tau,-1.5*R.dim1()) * exp(-greens);

  return greens;
}

double QMCWalker::calculateReverseGreensFunctionUmrigar93ImportanceSampling()
{
  double tau = Input->flags.dt;
  double greens = 1.0;

  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {
      // Find the nearest nucleus to this electron
      int nearestNucleus = 
	Input->Molecule.findClosestNucleusIndex(TrialWalker->R,electron);

      // Calculate the unit vector in the Z direction
      Array1D<double> zUnitVector(3);

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
	    (*TrialWalker->QMF.getModifiedGradPsiRatio())(electron,i);
	}

      Array1D<double> radialUnitVector(3);
      double radialComponentQF = 0.0;

      for(int i=0; i<3; i++)
	{
	  radialUnitVector(i) = 
	    (*TrialWalker->QMF.getModifiedGradPsiRatio())(electron,i) - 
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
	  double temp = TrialWalker->R(electron,i) - 
	    ( Input->Molecule.Atom_Positions(nearestNucleus,i) + 
	      radialCoordinate * radialUnitVector(i) + 
	      zCoordinate * zUnitVector(i) );

	  distance1Sq += temp*temp;
	}

      double g1 = pow(2*3.14159265359*tau,-1.5)*exp( -distance1Sq/(2*tau) );

      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
	{
	  double temp = TrialWalker->R(electron,i) -
	    Input->Molecule.Atom_Positions(nearestNucleus,i);

	  distance2Sq += temp*temp;
	}

      double g2 = pow(expParam,3.0) / 3.14159265359 *
	exp( -2.0 * expParam * sqrt( distance2Sq ) );

      greens *= (probabilityGaussianTypeMove * g1 + 
		 probabilitySlaterTypeMove * g2 );

    }

  return greens;
}

void QMCWalker::reweight_walker()
{
  double S_trial;
  double S_original;

  if( Input->flags.energy_modification_type == "none" )
    {
      S_trial    = Input->flags.energy_trial - 
	TrialWalker->QMF.getLocalEnergy();
      S_original = Input->flags.energy_trial - 
	OriginalWalker->QMF.getLocalEnergy();
    }
  else if( Input->flags.energy_modification_type == "modified_umrigar93" )
    {
      Array2D<double> * gradTrialModified = 
	TrialWalker->QMF.getModifiedGradPsiRatio();

      double temp = 0;
      for(int i=0; i<gradTrialModified->dim1(); i++ )
	{
	  for(int j=0; j<gradTrialModified->dim2(); j++)
	    {
	      temp += (*gradTrialModified)(i,j) *
		(*gradTrialModified)(i,j); 
	    }
	}

      double lengthGradTrialModified = sqrt(temp);

      Array2D<double> * gradTrialUnmodified = 
	TrialWalker->QMF.getGradPsiRatio();

      temp = 0;
      for(int i=0; i<gradTrialUnmodified->dim1(); i++ )
	{
	  for(int j=0; j<gradTrialUnmodified->dim2(); j++)
	    {
	      temp += (*gradTrialUnmodified)(i,j) *
		(*gradTrialUnmodified)(i,j); 
	    }
	}

      double lengthGradTrialUnmodified = sqrt(temp);


      S_trial = (Input->flags.energy_trial - 
		 TrialWalker->QMF.getLocalEnergy())*
	(lengthGradTrialModified/lengthGradTrialUnmodified);




      Array2D<double> * gradOriginalModified = 
	OriginalWalker->QMF.getModifiedGradPsiRatio();

      temp = 0;
      for(int i=0; i<gradOriginalModified->dim1(); i++ )
	{
	  for(int j=0; j<gradOriginalModified->dim2(); j++)
	    {
	      temp += (*gradOriginalModified)(i,j) *
		(*gradOriginalModified)(i,j); 
	    }
	}

      double lengthGradOriginalModified = sqrt(temp);

      Array2D<double> * gradOriginalUnmodified = 
	OriginalWalker->QMF.getGradPsiRatio();

      temp = 0;
      for(int i=0; i<gradOriginalUnmodified->dim1(); i++ )
	{
	  for(int j=0; j<gradOriginalUnmodified->dim2(); j++)
	    {
	      temp += (*gradOriginalUnmodified)(i,j) *
		(*gradOriginalUnmodified)(i,j); 
	    }
	}

      double lengthGradOriginalUnmodified = sqrt(temp);
      
      S_original = (Input->flags.energy_trial - 
		    OriginalWalker->QMF.getLocalEnergy()) *
	(lengthGradOriginalModified/lengthGradOriginalUnmodified);
    }
  else if( Input->flags.energy_modification_type == "umrigar93" )
    {
      Array2D<double> * gradTrialModified = 
	TrialWalker->QMF.getModifiedGradPsiRatio();

      double temp = 0;
      for(int i=0; i<gradTrialModified->dim1(); i++ )
	{
	  for(int j=0; j<gradTrialModified->dim2(); j++)
	    {
	      temp += (*gradTrialModified)(i,j) *
		(*gradTrialModified)(i,j); 
	    }
	}

      double lengthGradTrialModified = sqrt(temp);

      Array2D<double> * gradTrialUnmodified = 
	TrialWalker->QMF.getGradPsiRatio();

      temp = 0;
      for(int i=0; i<gradTrialUnmodified->dim1(); i++ )
	{
	  for(int j=0; j<gradTrialUnmodified->dim2(); j++)
	    {
	      temp += (*gradTrialUnmodified)(i,j) *
		(*gradTrialUnmodified)(i,j); 
	    }
	}

      double lengthGradTrialUnmodified = sqrt(temp);


      S_trial = (Input->flags.energy_trial - Input->flags.energy_estimated) +
	(Input->flags.energy_estimated - 
		 TrialWalker->QMF.getLocalEnergy())*
	(lengthGradTrialModified/lengthGradTrialUnmodified);



      Array2D<double> * gradOriginalModified = 
	OriginalWalker->QMF.getModifiedGradPsiRatio();

      temp = 0;
      for(int i=0; i<gradOriginalModified->dim1(); i++ )
	{
	  for(int j=0; j<gradOriginalModified->dim2(); j++)
	    {
	      temp += (*gradOriginalModified)(i,j) *
		(*gradOriginalModified)(i,j); 
	    }
	}

      double lengthGradOriginalModified = sqrt(temp);

      Array2D<double> * gradOriginalUnmodified = 
	OriginalWalker->QMF.getGradPsiRatio();

      temp = 0;
      for(int i=0; i<gradOriginalUnmodified->dim1(); i++ )
	{
	  for(int j=0; j<gradOriginalUnmodified->dim2(); j++)
	    {
	      temp += (*gradOriginalUnmodified)(i,j) *
		(*gradOriginalUnmodified)(i,j); 
	    }
	}

      double lengthGradOriginalUnmodified = sqrt(temp);
      
      S_original = (Input->flags.energy_trial - Input->flags.energy_estimated)
	+ (Input->flags.energy_estimated - 
	   OriginalWalker->QMF.getLocalEnergy()) *
	(lengthGradOriginalModified/lengthGradOriginalUnmodified);
    }  
  else
    {
      cerr << "ERROR: unknown energy modification method!" << endl;
      exit(0);
    }

  // determine the weighting factor dw so that the new weight = Weight*dw
  double dW = 0;

  if( Input->flags.run_type == "variational" )
    {
      // Keep all of the weights constant for VMC
      dW = 1.0;
    }
  else if( Input->flags.walker_reweighting_method == "simple_symmetric" )
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
      // supposedly between simple_symmetric and umrigar93_probability_weighted
      // in terms of its timestep and statistical performance

      double temp;

      if( move_accepted )
	{
	  temp = 0.5*(S_trial+S_original)*Input->flags.dt_effective;
	}
      else
	{
	  temp = S_original*Input->flags.dt_effective;
	}

      dW = exp(temp);
    }
  else if( Input->flags.walker_reweighting_method == 
	   "umrigar93_probability_weighted" )
    {
      // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993
      // Umrigar claims this has a small time step error and a small 
      // statistical error compared to simple_asymmetric and simple_symmetric

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

void QMCWalker::evaluate()
{
  QMF.evaluate(R);
}

void QMCWalker::calculateMoveAcceptanceProbability(double forwardGreens,
						   double reverseGreens)
{
  // This tells us the probability of accepting or rejecting a proposed move

  double PsiRatio = TrialWalker->QMF.getPsi()/OriginalWalker->QMF.getPsi();

  // increase the probability of accepting a move if the walker has not
  // moved in a long time

  double MoveOldWalkerFactor = 1.0;

  if( age > Input->flags.old_walker_acceptance_parameter )
    {
      MoveOldWalkerFactor = 
	pow(1.1,age-Input->flags.old_walker_acceptance_parameter);
    }

  // calculate the probability of accepting the trial move
  double p = PsiRatio * PsiRatio * reverseGreens / forwardGreens * 
    MoveOldWalkerFactor;

  // Force the walker to move if it is too old
  if( age > 2*Input->flags.old_walker_acceptance_parameter )
    {
      p = 1;
    }


  // if the aratio is NaN then reject the move
  if( isnan(p) != 0 )
    {
      cerr << "WARNING: Rejecting trial walker with NaN aratio!" << endl;
      p = 0.0;
    }

  // if the trial position is singular reject the move
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Rejecting singular trial walker!" << endl;
      p = 0.0;
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
      //	   << " age = " << age << endl;
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
  QMF.initialize(Input);

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

   QMF.writeCorrelatedSamplingConfiguration(strm);
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
  strm << "\t<Elocal> \n\t\t" << QMF.getLocalEnergy() 
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

  evaluate();
}

void QMCWalker::initializeWalkerPosition()
{
  QMCInitializeWalker * IW = 
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
				   Input->flags.walker_initialization_method);

  R = IW->initializeWalkerPosition();
  evaluate();

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
      evaluate();
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
  localEnergy = p * TrialWalker->QMF.getLocalEnergy() + 
    q * OriginalWalker->QMF.getLocalEnergy();

  // Calculate the kinetic energy...
  kineticEnergy = p * TrialWalker->QMF.getKineticEnergy() +
    q * OriginalWalker->QMF.getKineticEnergy();

  // Calculate the potential energy
  potentialEnergy = p * TrialWalker->QMF.getPotentialEnergy() +
    q * OriginalWalker->QMF.getPotentialEnergy();

  // Calculate the acceptance probability
  AcceptanceProbability = p;

  // Calculate the DistanceMovedAccepted this is the average distance
  // moved on a step
  distanceMovedAccepted = p * dR2;

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
}

bool QMCWalker::isSingular()
{
  return QMF.isSingular();
}


double QMCWalker::getLocalEnergyEstimator()
{
  return localEnergy;
}




