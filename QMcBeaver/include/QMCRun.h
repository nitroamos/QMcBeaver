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
#include <queue>

#include "QMCWalker.h"
#include "QMCEquilibrationArray.h"
#include "QMCProperties.h"
#include "QMCStopwatches.h"
#include "QMCHartreeFock.h"
#include "QMCInitializeWalkerFactory.h"
#include "QMCFunctionsFactory.h"

using namespace std;

/**
  Collection of walkers (QMCWalker) with the functionality to do the basic
  operations from which a QMC algorithm is built.
*/

class QMCRun
{
 public:

  ~QMCRun();

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
     Initalize the QMC Function. This is available
     separately because it will need to be reset everytime
     we change the number of parameters we're optimizing.
  */
  void initializeFunction();

  /**
    Sets all of the data in the object to zero.
  */
  void zeroOut();

  /**
    Propagate the QMC calculation one time step forward.

    @param iteration our current iteration number. If equilibrating, less than zero
    @return if the calculation is still ok, it will return true
  */
  bool step(bool writeConfigs, int iteration);

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
    Gets the statistics for the future walking properties that have been calculated.
    @return statistics for the properties that have been calculated.
  */
  QMCPropertyArrays * getFWProperties();
  
  /**
    Gets the statistics for the future walking properties that have
    been calculated at this time step.
    @return statistics for the properties that have been calculated at this
    time step.
  */
  QMCPropertyArrays * getFWTimeStepProperties();

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
    Calculates the distances between all pairs of particles and records them in
    histograms.
  */
  void calculateElectronDensities();

  /** 
    Gets a pointer to the the histogram of distances between parallel spin
    electrons.
    @return parallel spin histogram.
  */
  Array1D<double>* getPllSpinHistogram();

  /**
    Gets a pointer to the histogram of distances between opposite spin 
    electrons.
    @return opposite spin histogram.
  */
  Array1D<double>* getOppSpinHistogram();

  /**
    Gets a pointer to the alpha one electron density histograms.
    @return alpha histograms.
  */
  Array1D< Array1D<double> >* getAlphaHistograms();

  /**
    Gets a pointer to the beta one electron density histograms.
    @return beta histograms.
  */
  Array1D< Array1D<double> >* getBetaHistograms();

  /**
    Gets a pointer to the 2D x coordinate parallel pair correlation diagram.
    @return 2D x coordinate parallel pair correlation diagram.
  */
  Array1D< Array1D<double> >* getPllxCorrelationDiagram();

  /**
    Gets a pointer to the 2D y coordinate parallel pair correlation diagram.
    @return 2D y coordinate parallel pair correlation diagram.
  */
  Array1D< Array1D<double> >* getPllyCorrelationDiagram();

  /**
    Gets a pointer to the 2D z coordinate parallel pair correlation diagram.
    @return 2D z coordinate parallel pair correlation diagram.
  */
  Array1D< Array1D<double> >* getPllzCorrelationDiagram();

  /**
    Gets a pointer to the 2D x coordinate opposite pair correlation diagram.
    @return 2D x coordinate opposite pair correlation diagram.
  */
  Array1D< Array1D<double> >* getOppxCorrelationDiagram();

  /**
    Gets a pointer to the 2D y coordinate opposite pair correlation diagram.
    @return 2D y coordinate opposite pair correlation diagram.
  */
  Array1D< Array1D<double> >* getOppyCorrelationDiagram();

  /**
    Gets a pointer to the 2D z coordinate opposite pair correlation diagram.
    @return 2D z coordinate opposite pair correlation diagram.
  */
  Array1D< Array1D<double> >* getOppzCorrelationDiagram();

  /**
    Gets the size of a bin in the electron density histograms.
  */
  double getdr();

  /**
    Writes the state of this object to an XML stream.
    @param strm XML stream
  */
  void toXML(ostream& strm);  

  /**
    Reads the state of this object from an XML stream.
    @param strm XML stream
    @return whether the read was successful
  */
  bool readXML(istream& strm);

  void updateHFPotential();

  int getNumTimers();
  void aggregateTimers(Array1D<Stopwatch> & timers,
		       int & idx);

private:

  Array1D<QMCWalkerData *> dataPointers;
  Array1D<Array2D<double> *> rPointers;

  /**
    List of all the walkers.
  */
  list<QMCWalker> wlist;

  QMCFunctions * QMF;

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
    The statistics for this group of walkers if QMCEquilibrationArray is not
    used.
  */
  QMCPropertyArrays fwProperties;

  /**
    The statistics for these walkers at this time step.
  */
  QMCProperties timeStepProperties;
  
  /**
    The statistics for these walkers at this time step, using
    future walking.
  */
  QMCPropertyArrays fwTimeStepProperties;

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
     This queue is only used by calculatePopulationSizeBiasCorrectionFactor used
     to store the history of populationSizeBiasCorrectionFactor factors. See UNR93
     for more details.
   */
  queue<double> correctionDivisor;

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
     A new kind of branching where the number of walkers is a constant.
     Basically, bad walkers are replaced with good walkers. I'm not sure
     if the implementation is 100% correct.

     See Phys Rev E, vol 16, no 4, pg 4566
     Diffusion Monte Carlo methods with a fixed number of walkers
     Assaraf, Caffarel, Khelif
   */
  void ack_reconfiguration();

  /**
    Proposes trial walker moves and accepts or rejects them. This function
    steps QMCWalker in two stages. First, it calculates the forward green's
    function by calling initializePropagation, and then it calculates the
    reverse green's function by calling processPropgatation. The purpose of
    this is to enable this function to call QMCFunction for several walkers
    at once (i.e. in chunks) by pausing the rest of the tasks that QMCWalker
    has to do.

    @param iteration our current iteration number. If equilibrating, less than zero
  */
  void propagateWalkers(bool writeConfigs, int iteration);

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

    This factor is described in Umrigar, Nightingale, Runge 1993 paper (UNR93).
  */
  void calculatePopulationSizeBiasCorrectionFactor();

  /**
     Growth = Births - Deaths
   */
  int growthRate;

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
  Array1D<double> pllSpinHistogram;

  /**
    The histogram of distances between opposite spin electrons.
  */
  Array1D<double> oppSpinHistogram;

  /**
    The alpha one electron density histograms.  There is one histogram for each
    unique nucleus.
  */
  Array1D< Array1D<double> > alphaHistograms;

  /**
    The beta one electron density histograms.  There is one histogram for each
    unique nucleus.
  */
  Array1D< Array1D<double> > betaHistograms;

  /**
    The 2D x coordinate parallel pair correlation diagram.
  */
  Array1D< Array1D<double> > pllxCorrelationDiagram;

  /**
    The 2D y coordinate parallel pair correlation diagram.
  */
  Array1D< Array1D<double> > pllyCorrelationDiagram;

  /**
    The 2D z coordinate parallel pair correlation diagram.
  */
  Array1D< Array1D<double> > pllzCorrelationDiagram;

  /**
    The 2D x coordinate opposite pair correlation diagram.
  */
  Array1D< Array1D<double> > oppxCorrelationDiagram;

  /**
    The 2D y coordinate opposite pair correlation diagram.
  */
  Array1D< Array1D<double> > oppyCorrelationDiagram;

  /**
    The 2D z coordinate opposite pair correlation diagram.
  */
  Array1D< Array1D<double> > oppzCorrelationDiagram;

  /**
    These objects allow HF calculations to be done with QMC.
  */
  QMCHartreeFock HartreeFock;

  Array1D<Stopwatch> swTimers;
};

#endif
