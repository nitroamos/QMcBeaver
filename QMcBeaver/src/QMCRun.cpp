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

#include "QMCRun.h"

QMCRun::QMCRun()
{
  populationSizeBiasCorrectionFactor = 1.0;
}

void QMCRun::propagateWalkers()
{
  // Propigate all of the walkers

  for(list<QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
    {
      wp->propagateWalker(); 
    }
}

void QMCRun::branchWalkers()
{
  if( Input->flags.run_type == "variational" )
    {
      // No branching in VMC
    }
  else if( Input->flags.branching_method == "non_branching" )
    {
      // Non branching DMC -- The weights are varied instead of branching
    }
  else if( Input->flags.branching_method == "unit_weight_branching" )
    {
      // This is the method described in Lester and Hammonds book
      // It is Chip Kent's (my) experience, that this method 
      // often has an exponentially growing or shrinking population
      // and has a large time step error

      unitWeightBranching();
    }
  else if( Input->flags.branching_method == "nonunit_weight_branching" )
    {
      // This method is described in Umrigar, Nightingale, and Runge 
      // JCP 99(4) 2865; 1993 and elsewhere.  A walker is reweighted
      // and is divided if it's weight exceeds a threshold and is fused
      // with another walker if it's weight falls below a threshold.

      nonunitWeightBranching();
    }
  else
    {
      cerr << "ERROR: Unknown branching method!" << endl;
      exit(0);
    }
}

void QMCRun::zeroOut()
{
  Properties.zeroOut();
}

void QMCRun::initialize(QMCInput *INPUT)
{
  Input = INPUT;

  Properties.zeroOut();
}

void QMCRun::randomlyInitializeWalkers()
{
  for(int i=0;i<Input->flags.number_of_walkers;i++)
    {
      QMCWalker w;
      w.initialize(Input);

      // randomly initialize the walkers position
      w.initializeWalkerPosition();
   
      wlist.push_back(w);
    }
}

void QMCRun::calculateObservables()
{
  // Get the observables calculated in each walker

  // Pre block all of the statistics from one time step together

  QMCProperties timeStepProperties;

  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      wp->calculateObservables( timeStepProperties );
    }

  // Add the pre blocked data from this time step to the accumulated
  // statistics

  double totalWeights = getWeights() * populationSizeBiasCorrectionFactor;

  // Calculate the Energy ...
  Properties.energy.newSample( 
	timeStepProperties.energy.getAverage(), totalWeights );

  // Calculate the Kinetic Energy
  Properties.kineticEnergy.newSample( 
        timeStepProperties.kineticEnergy.getAverage(), totalWeights );

  // Calculate the Potential Energy
  Properties.potentialEnergy.newSample( 
	timeStepProperties.potentialEnergy.getAverage(), totalWeights );

  // Calculate the Acceptance Probability ...
  Properties.acceptanceProbability.newSample( 
	timeStepProperties.acceptanceProbability.getAverage(), totalWeights );

  // Calculate the DistanceMovedAccepted this is the average distance
  // moved on a step
  Properties.distanceMovedAccepted.newSample( 
	timeStepProperties.distanceMovedAccepted.getAverage(), totalWeights );

  // Calculate the DistanceMovedTrial this is the average step length for
  // a trial move
  Properties.distanceMovedTrial.newSample( 
	timeStepProperties.distanceMovedTrial.getAverage(), totalWeights );

  // Calculate the log of the weights
  Properties.logWeights.newSample(
        timeStepProperties.logWeights.getAverage(), getNumberOfWalkers() );
}

void QMCRun::writeEnergies(ostream& strm)
{
  // Get the energy calculated in each walker and write out

  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      strm << wp->getLocalEnergyEstimator() << endl;
    }
}

void QMCRun::writeCorrelatedSamplingConfigurations(ostream& strm)  
{
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      wp->writeCorrelatedSamplingConfiguration(strm);
    }
}

void QMCRun::unitWeightBranching()
{
  // Make a list of walkers to add and delete
  list<QMCWalker> WalkersToAdd;
  list<list<QMCWalker>::iterator> WalkersToDelete;

  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      int times_to_branch = int(wp->getWeight() + ran1(&Input->flags.iseed))-1;

      // Set the walkers weight back to unity
      wp->setWeight(1.0);

      if( times_to_branch == 0 )
	{
	  // Don't need to copy the walker
	}
      else if( times_to_branch < 0 )
	{
	  // Mark the walker to delete

	  WalkersToDelete.push_back( wp );
	}
      else
	{
	  // Mark the walker to add times_to_branch times

	  for(int i=0; i<times_to_branch; i++)
	    {
	      WalkersToAdd.push_back( *wp );
	    }
	}
    }

  // Delete the unneeded walkers
#ifdef DEBUG_BRANCHING
  if( WalkersToDelete.size() > 0 )
    {
      cerr << "Deleting " << WalkersToDelete.size() << " walker(s) from "
	   << wlist.size() << " walkers" << endl;
    }
#endif

  for( list< list<QMCWalker>::iterator >::iterator  
	 wtd=WalkersToDelete.begin(); wtd!=WalkersToDelete.end(); ++wtd)
    {
      wlist.erase( *wtd );
    }

  // Add the new walkers
#ifdef DEBUG_BRANCHING
  if( WalkersToAdd.size() > 0 )
    {
      cerr << "Adding " << WalkersToAdd.size() << " walker(s) to "
	   << wlist.size() << " walkers" << endl;
    }
#endif

  wlist.splice( wlist.end(), WalkersToAdd );
}

void QMCRun::nonunitWeightBranching()
{
  // Make a list of walkers to add and delete
  list<QMCWalker> WalkersToAdd;
  list<list<QMCWalker>::iterator> WalkersToDelete;

  list<QMCWalker>::iterator TempWalkerToDelete; 
  bool isTempWalkerToDelete = false;

  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      if( wp->getWeight() > Input->flags.branching_threshold )
	{
	  // Split this walker into two; each with half the weight

	  wp->setWeight( wp->getWeight()/2.0 );

	  WalkersToAdd.push_back( *wp );
	}
      else if( wp->getWeight() < Input->flags.fusion_threshold )
	{
	  // This walker needs to be fused with another

	  if( !isTempWalkerToDelete )
	    {
	      // No walker to immediately fuse this with so wait for one

	      TempWalkerToDelete = wp;
	      isTempWalkerToDelete = true;
	    }
	  else
	    {
	      double weight1 = TempWalkerToDelete->getWeight();
	      double weight2 = wp->getWeight();
	      double weight3 = weight1 + weight2;

	      if( ran1(&Input->flags.iseed) < weight1/weight3 )
		{
		  // Keep TempWalkerToDelete and delete other

		  TempWalkerToDelete->setWeight( weight3 );
		  WalkersToDelete.push_back( wp );
		}
	      else
		{
		  // opposite of above

		  wp->setWeight( weight3 );
		  WalkersToDelete.push_back( TempWalkerToDelete );
		}

	      isTempWalkerToDelete = false;
	    }
	}
      else
	{
	  // walker is ok and don't do anything to it
	}
    }

  // Delete the unneeded walkers
#ifdef DEBUG_BRANCHING
  if( WalkersToDelete.size() > 0 )
    {
      cerr << "Deleting " << WalkersToDelete.size() 
	   << " walker(s) from " << wlist.size() << " walkers" << endl;
    }
#endif

  for(list< list<QMCWalker>::iterator >::iterator wtd=WalkersToDelete.begin();
      wtd!=WalkersToDelete.end(); ++wtd)
    {
      wlist.erase( *wtd );
    }

  // Add the new walkers
#ifdef DEBUG_BRANCHING
  if( WalkersToAdd.size() > 0 )
    {
      cerr << "Adding " << WalkersToAdd.size() << " walker(s) to "
	   << wlist.size() << " walkers" << endl;
    }
#endif

  wlist.splice( wlist.end(), WalkersToAdd );
}

double QMCRun::getWeights()
{
  // determine the total of all the walker weights
  double total_weights = 0;

  for( list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp )
    {
      total_weights += wp->getWeight();
    }

  return total_weights;
}

void QMCRun::toXML(ostream& strm)
{
  strm << "<populationSizeBiasCorrectionFactor>\n\t" 
       << populationSizeBiasCorrectionFactor 
       << "\n</populationSizeBiasCorrectionFactor>" << endl;

  // writes out the properties
  Properties.toXML(strm);

  //prints all the walkers
  strm << "<walkers>" << endl;

  for(list <QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
    {
      wp->toXML(strm);
    }

  strm << "</walkers>\n" << endl;
}

void QMCRun::readXML(istream& strm)
{
  string temp;

  // read populationSizeBiasCorrectionFactor
  strm >> temp >> temp;
  populationSizeBiasCorrectionFactor = atof(temp.c_str());
  strm >> temp;

  // read the properties
  Properties.readXML(strm);

  // read the walkers

  strm >> temp;

  for(int i=0;i<Input->flags.number_of_walkers;i++)
    {
      QMCWalker w;
      w.initialize(Input);

      // read a walker
      w.readXML(strm);

      wlist.push_back(w);
    }

  strm >> temp;

}

QMCProperties * QMCRun::getProperties()
{
  return &Properties;
}

int QMCRun::getNumberOfWalkers()
{
  return wlist.size();
}

void QMCRun::step()
{
  propagateWalkers();
  calculatePopulationSizeBiasCorrectionFactor();
  calculateObservables();
  branchWalkers();
}

void QMCRun::calculatePopulationSizeBiasCorrectionFactor()
{
  if( Input->flags.correct_population_size_bias && 
      Input->flags.run_type != "variational" )
    {
      double temp = Input->flags.energy_trial - 
	Input->flags.energy_estimated_original;
      
      temp *= -Input->flags.dt;
      
      temp = exp(temp);

      populationSizeBiasCorrectionFactor *= temp;
    }
}



