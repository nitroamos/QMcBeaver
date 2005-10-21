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

#ifndef QMCRUN_H
#define QMCRUN_H

#include <fstream>
#include <list>

#include "QMCWalker.h"
#include "QMCEquilibrationArray.h"
#include "QMCProperties.h"
#include "QMCStopwatches.h"
#include "QMCHartreeFock.h"

using namespace std;

/**
  Collection of walkers (QMCWalker) with the functionality to do the basic
  operations from which a QMC algorithm is built.
*/

class QMCRun
{
 public:
 
  /**
    Creates an uninitialized instance of this class.
  */
  QMCRun();

  /**
    Initializes this object.
    @param input input data for the calculation
  */
  void initialize(QMCInput *input);  

  /**
    Sets all of the data in the object to zero.
  */
  void zeroOut();

  /**
    Propagate the QMC calculation one time step forward.
  */
  void step(bool writeConfigs);

  /**
    Gets the statistics for the properties that have been calculated.
    @return statistics for the properties that have been calculated.
  */
  QMCProperties * getProperties();

  /**
    Gets the statistics for the properties that have been calculated at this 
    time step.
    @return statistics for the properties that have been calculated at this
    time step.
  */
  QMCProperties * getTimeStepProperties();

  /**
    Starts the timers in the EquilibrationArray.
  */
  void startTimers();

  /**
    Stops the timers in the EquilibrationArray.
  */
  void stopTimers();

  /**
    Gets the the propagation stopwatch from the appropriate element of the
    EquilibrationArray.
    @return propagation Stopwatch from the appropriate element of the 
    EquilibrationArray.
  */
  Stopwatch * getPropagationStopwatch();

  /** 
    Gets the equilibration stopwatch from the appropriate element of the
    EquilibrationArray.
    @return equilibration Stopwatch from the appropriate element of the 
    EquilibrationArray.
  */
  Stopwatch * getEquilibrationStopwatch();

  /**
    Gets the total statistical weights for all the current living walkers.
    @return total weights for current walkers.
  */
  double getWeights();

  /** 
    Gets the factor that corrects the population size bias in the statistics.
    @return population size bias correction factor
  */
  double getPopulationSizeBiasCorrectionFactor();

  /**
    Gets the current number of walkers.
    @return number of walkers.
  */
  int getNumberOfWalkers();

  /**
    Generates all of the walkers by initializing the electronic 
    configurations for the walkers using an algorithm from 
    QMCInitializeWalkerFactory.  
  */
  void randomlyInitializeWalkers();

  /**
    Writes the energies of all the walkers to a stream.
    @param strm stream to write energies to.
  */
  void writeEnergies(ostream& strm);

  /**
    Writes the state of this group of walkers to a stream in a 
    format that is suitable for correlated sampling calculations.  
    This writes out more information than 
    <code>toXML</code> so that parts of the wavefunction do not have to
    be reevaluated every time properties are calculated using correlated
    sampling.
    @param strm stream to write correlated sampling information to.
  */ 
  void writeCorrelatedSamplingConfigurations(ostream& strm);  

  /**
    Calculates the distances between all pairs of electrons and records them in
    parallel and opposite spin histograms.
  */
  void calculateElectronDensities();

  /**
    Writes the parallel and opposite spin pair distance histograms to files.
  */
  void writeElectronDensityHistograms();
  
  /**
    Writes the state of this object to an XML stream.
    @param strm XML stream
  */
  void toXML(ostream& strm);  

  /**
    Reads the state of this object from an XML stream.
    @param strm XML stream
  */
  void readXML(istream& strm);

  void updateHFPotential();

private:

  /**
    List of all the walkers.
  */
  list<QMCWalker> wlist;

  QMCFunctions QMF;

  /**
    The array of Decorrelation objects for this group of walkers.
  */
  QMCEquilibrationArray EquilibrationArray;

  /**
    The statistics for this group of walkers if QMCEquilibrationArray is not
    used.
  */
  QMCProperties Properties;

  /**
    The statistics for these walkers at this time step.
  */
  QMCProperties timeStepProperties;

  /**
    Input data to control the calculation.
  */
  QMCInput *Input;

  /**
    A factor that is used in removing the bias introduced into a calculation
    by using a finite number of walkers.  The bias is only
    present in calculations that use branching.  
  */
  double populationSizeBiasCorrectionFactor;

  /**
    During a DMC calculation branch the walkers while keeping the weights
    equal to one.  This was developed by Lester.  It is of Chip Kent's (my)
    experience that this method often has an exponentially growing or 
    shrinking population and has a large time step error.
  */
  void unitWeightBranching();

  /**
    During a DMC calculation branch a walker if its weight exceeds a 
    threshold and fuse walkers when their weights fall below a threshold.
  */
  void nonunitWeightBranching();

  /**
    Proposes trial walker moves and accepts or rejects them. This function
    steps QMCWalker in two stages. First, it calculates the forward green's
    function by calling initializePropagation, and then it calculates the
    reverse green's function by calling processPropgatation. The purpose of
    this is to enable this function to call QMCFunction for several walkers
    at once (i.e. in chunks) by pausing the rest of the tasks that QMCWalker
    has to do.
  */
  void propagateWalkers(bool writeConfigs);

  /**
    Creates and destroys walkers based on their weights.
  */
  void branchWalkers();

  /**
    Adds the observable data calculated during this step to the data that
    has already been recorded.
  */
  void calculateObservables();

  /**
    Calculates a factor that is used in removing the bias introduced into 
    a calculation by using a finite number of walkers.  The bias is only
    present in calculations that use branching.  
  */
  void calculatePopulationSizeBiasCorrectionFactor();

  /**
    The maximum distance between a pair of particles that will be recorded in
    the histograms.
  */
  double max_pair_distance;

  /**
    The size of a bin in the pair histograms.
  */
  double dr;

  /** 
    The histogram of distances between parallel spin electrons.
  */
  Array1D<double> pll_spin_histogram;

  /**
    The histogram of distances between opposite spin electrons.
  */
  Array1D<double> opp_spin_histogram;

  /**
    The alpha one electron density histograms.  There is one histogram for each
    unique nucleus.
  */
  Array2D<double> alpha_density_histogram;

  /**
    The beta one electron density histograms.  There is one histogram for each
    unique nucleus.
  */
  Array2D<double> beta_density_histogram;

  /**
    The total weight of all samples recorded in the histograms.
  */
  double total_sample_weight;

  /**
    These objects allow HF calculations to be done with QMC.
  */
  QMCHartreeFock HartreeFock;
};

#endif
