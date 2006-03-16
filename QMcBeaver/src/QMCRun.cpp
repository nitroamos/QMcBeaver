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

void QMCRun::propagateWalkers(bool writeConfigs)
{
  int count = 0;
  int index = 0;
  int wpp = Input->flags.walkers_per_pass;
  // Propagate all of the walkers
  Array1D<QMCWalkerData *> dataPointers = 0;
  Array1D<Array2D<double> * > rPointers = 0;
  dataPointers.allocate(wpp);
  rPointers.allocate(wpp);
  
  /*The point here is to collect WALKERS_PER_PASS amount of walkers
    to process at once. Once we have that many (or the last remaining), we 
    finally run QMF.evaluate which fills up the dataPointers data structure 
    with values.  The actual data for dataPointers is stored in each walker, so
    initializePropagation asks for a pointer to it. The collection of pointers 
    is then passed around to everybody.
    The advantage to filling the array dynamically is that we don't have to 
    worry about branching -- walkers being deleted and created.
  */
  for(list<QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
    {
      wp->initializePropagation(dataPointers(index),rPointers(index));
      count++;
      index = count%wpp;
      if(index == 0 || count == (int)(wlist.size()))
        {
          QMF.evaluate(dataPointers, rPointers,
                       index==0?wpp:index, writeConfigs);
        }
    }
    
  /*At this point, all of the dataPointers have been filled, so we
    can tell each walker to finish processing each move. Note this will
    give different answers than before (when all processing was done before
    the next walker was evaluated) because random numbers are drawn in a
    different order.
   */
  for(list<QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
    {
      wp->processPropagation(QMF);
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
      // and is divided if its weight exceeds a threshold and is fused
      // with another walker if its weight falls below a threshold.
      
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
  if (Input->flags.use_equilibration_array == 1) EquilibrationArray.zeroOut();
  else Properties.zeroOut();
}

void QMCRun::initialize(QMCInput *INPUT)
{
  Input = INPUT;
  QMF.initialize(Input,&HartreeFock);
  
  if (Input->flags.calculate_bf_density == 1)
    {
      bool calcDensity = true;
      timeStepProperties.setCalcDensity(calcDensity,
                                        Input->WF.getNumberBasisFunctions());
      if (Input->flags.use_equilibration_array == 1)
        EquilibrationArray.setCalcDensity(calcDensity,
                                          Input->WF.getNumberBasisFunctions());
      else
        Properties.setCalcDensity(calcDensity,
                                  Input->WF.getNumberBasisFunctions());
    }
    
  if (Input->flags.use_equilibration_array == 1) EquilibrationArray.zeroOut();
  else Properties.zeroOut();
  
  if (Input->flags.write_electron_densities == 1)
    {
      // The pair density for each pair of particles is recorded in a 
      // histogram.  There are histograms for opposite and parallel spin
      // electrons, and for alpha and beta electrons and each unique nucleus in
      // the molecule.  The maximum distance is 20 divided by the largest
      // atomic charge in the molecule, and there are 5,000 bins.
      
      int max_Z = 0;
      
      for (int i=0; i<Input->Molecule.getNumberAtoms(); i++)
        if ( Input->Molecule.Z(i) > max_Z )
          max_Z = Input->Molecule.Z(i);
          
      max_pair_distance = 20.0/max_Z;
      dr = max_pair_distance/5000;
      
      // Histograms for parallel and opposite spin electron densities.
	  int nalpha = Input->WF.getNumberAlphaElectrons();
	  int nbeta = Input->WF.getNumberBetaElectrons();

	  if (nalpha > 1 || nbeta > 1)
	    {
	      pllSpinHistogram.allocate(5000);
	      pllSpinHistogram = 0.0;
	    }
	  
	  if (nalpha > 0 && nbeta > 0)
	    {
	      oppSpinHistogram.allocate(5000);
	      oppSpinHistogram = 0.0;
	    }

      // Histograms for the one electron densities.
      int nucleiTypes = Input->Molecule.NucleiTypes.dim1();

	  if (nalpha > 0)
	    {
	      alphaHistograms.allocate(nucleiTypes);
	      for (int k=0; k<nucleiTypes; k++)
		{
		  alphaHistograms(k).allocate(5000);
		  alphaHistograms(k) = 0.0;
		}
	    }

	  if (nbeta > 0)
	    {
	      betaHistograms.allocate(nucleiTypes);
	      for (int k=0; k<nucleiTypes; k++)
		{
		  betaHistograms(k).allocate(5000);
		  betaHistograms(k) = 0.0;
		}
	    }
    }

  if (Input->flags.use_hf_potential == 1)
    HartreeFock.Initialize(Input);
}

void QMCRun::randomlyInitializeWalkers()
{
  Array2D<double> temp_R;
  int initialization_try = 0;

  QMCInitializeWalker * IW =
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
        Input->flags.walker_initialization_method);
  
  for (int i=0; i<Input->flags.number_of_walkers; i++)
    {
      QMCWalker w;
      w.initialize(Input);

      temp_R = IW->initializeWalkerPosition();

      w.setR(temp_R);
      QMF.evaluate(*w.getR(),*w.getWalkerData());
  
      initialization_try = 1;
      while( w.isSingular() )
	{
	  cerr << "Regenerating Walker..." << endl;
      
	  if( initialization_try > 10 )
	    {
	      cerr << "ERROR: 10 consecutive singular configurations while "
		   << "trying to initialize walker!" << endl;
	      exit(0);
	    }
        
	  temp_R = IW->initializeWalkerPosition();
	  w.setR(temp_R);
	  QMF.evaluate(*w.getR(),*w.getWalkerData());
	  initialization_try++;
	}
      wlist.push_back(w);
    }
  delete IW;
  IW = 0;
}

void QMCRun::calculateObservables()
{
  // Get the observables calculated in each walker
  // Pre block all of the statistics from one time step together
  
  timeStepProperties.zeroOut();
  
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    wp->calculateObservables( timeStepProperties );
    
  // Add the pre blocked data from this time step to the accumulated
  // statistics
  
  double totalWeights = getWeights() * populationSizeBiasCorrectionFactor;
  
  if (Input->flags.use_equilibration_array == 1)
    EquilibrationArray.newSample(&timeStepProperties, totalWeights,
                                 getNumberOfWalkers());
  else
    Properties.newSample(&timeStepProperties, totalWeights,
                         getNumberOfWalkers());
}

void QMCRun::writeEnergies(ostream& strm)
{
  // Get the energy calculated in each walker and write out
  
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
    strm << wp->getLocalEnergyEstimator() << "\t" << wp->getWeight() << endl;
}

void QMCRun::calculateElectronDensities()
{
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
        wp->calculateElectronDensities(max_pair_distance, dr, pllSpinHistogram,
		            oppSpinHistogram, alphaHistograms, betaHistograms);
}

Array1D<double>* QMCRun::getPllSpinHistogram()
{
  return &pllSpinHistogram;
}

Array1D<double>* QMCRun::getOppSpinHistogram()
{
  return &oppSpinHistogram;
}

Array1D< Array1D<double> >* QMCRun::getAlphaHistograms()
{
  return &alphaHistograms;
}

Array1D< Array1D<double> >* QMCRun::getBetaHistograms()
{
  return &betaHistograms;
}

double QMCRun::getdr()
{
  return dr;
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
  double total_weights = 0.0;
  
  for( list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp )
    total_weights += wp->getWeight();
    
  return total_weights;
}

double QMCRun::getPopulationSizeBiasCorrectionFactor()
{
  return populationSizeBiasCorrectionFactor;
}

void QMCRun::toXML(ostream& strm)
{
  strm << "<populationSizeBiasCorrectionFactor>\n\t"
  << populationSizeBiasCorrectionFactor
  << "\n</populationSizeBiasCorrectionFactor>" << endl;
  
  // writes out the properties
  if (Input->flags.use_equilibration_array == 1)
    {
      EquilibrationArray.toXML(strm);
    }
  else Properties.toXML(strm);
  
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
  if (Input->flags.use_equilibration_array == 1)
    {
      EquilibrationArray.readXML(strm);
    }
  else Properties.readXML(strm);
  
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

void QMCRun::startTimers()
{
  EquilibrationArray.startTimers();
}

void QMCRun::stopTimers()
{
  EquilibrationArray.stopTimers();
}

Stopwatch * QMCRun::getPropagationStopwatch()
{
  return EquilibrationArray.getPropagationStopwatch();
}

Stopwatch * QMCRun::getEquilibrationStopwatch()
{
  return EquilibrationArray.getEquilibrationStopwatch();
}

QMCProperties * QMCRun::getTimeStepProperties()
{
  return &timeStepProperties;
}

QMCProperties * QMCRun::getProperties()
{
  if (Input->flags.use_equilibration_array == 1)
    return EquilibrationArray.chooseDecorrObject();
  else
    return &Properties;
}

int QMCRun::getNumberOfWalkers()
{
  return (int)wlist.size();
}

void QMCRun::step(bool writeConfigs)
{
  propagateWalkers(writeConfigs);
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

// Updates hartree-fock potential, adding all electrons over all walkers to the
// object.
void QMCRun::updateHFPotential()
{
  int numelecs = Input->WF.getNumberElectrons();

  for (list<QMCWalker>::iterator wp = wlist.begin(); wp != wlist.end(); ++wp)
    {
      Array2D<double> positions = *(wp->getR());
      for (int i = 0; i < numelecs; i++)
	HartreeFock.AddElectron(i,wp->getWeight(),positions(i,0),
				positions(i,1), positions(i, 2));
      HartreeFock.IncrementSample();
    }
}
