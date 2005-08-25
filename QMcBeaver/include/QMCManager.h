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

#ifndef QMCManager_H
#define QMCManager_H

#define WORK_STEP 0
#define REDUCE 1
#define TERMINATE 2
#define SYNCHRONIZE 3

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "QMCRun.h"
#include "QMCDerivativeProperties.h"
#include "QMCCorrelatedSamplingVMCOptimization.h"
#include "QMCCopyright.h"
#include "QMCProperties.h"
#include "QMCProperty.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

using namespace std;

/**
   Controls the major sections of a QMC calculation.  This allows a 
   QMC calculation to be run and parameters to be optimized.  
*/

class QMCManager
{
 public:
  /**
     Creates an uninitialized instance of this class.
  */
  QMCManager();

  /**
     Destroys this object, cleans up the memory, and closes all open streams.
  */
  ~QMCManager();

  /**
     Initializes this object and loads the input data for the calculation.
     @param argc number of command line arguments.
     @param argv command line arguments.
  */
  void initialize(int argc, char **argv);

  /**
     Prepares the calculation to terminate.
  */
  void finalize();

  /**
     Performs a QMC calculation.  The specifics of the calculation are
     prescribed in the input.
  */
  void run();

  /**
     Optimizes the parameters in a variational QMC (VMC) calculation using 
     the correlated sampling method.
  */
  void optimize();

  /**
     Zeroes out all of the statistical data calculated by this object.
  */
  void zeroOut();

  /**
     Writes the restart file for the calculation.
  */
  void writeRestart();

  /**
     Writes the basis function density for the calculation.
  */
  void writeBFDensity();

  /**
     Writes the timing data to a stream. This is only valid after
     <code>finalize</code> is called and only on the root node.
     @param strm stream to write timing information to.
  */
  void writeTimingData(ostream & strm);

  /**
     Gets the input data for the calculation.
     @return input data for the calculation.
  */
  QMCInput * getInputData();

  /**
     Gets the stream for outputting results from a calculation.
     @return output stream for results.
  */
  ostream * getResultsOutputStream();

  /**
     Writes the current QMC results calculated by this object to
     an output stream in a human readable format.
  */
  friend ostream & operator<<(ostream & strm, QMCManager & rhs);

 private:
  /**
     Walkers plus the basic functionality needed to build a QMC calculation.
  */
  QMCRun QMCnode;

  /**
     Input data for the calculation.
  */
  QMCInput Input;

  /**
     Local code timers for this node.
  */
  QMCStopwatches localTimers;  

  /**
     Global code timers for the calculation.  This is only valid after
     <code>finalize</code> is called and only on the root node.
  */
  QMCStopwatches globalTimers;

  /**
     <code>true</code> if the calculation is equilibrating and 
     <code>false</code> otherwise.
  */
  bool equilibrating;

  /**
     <code>true</code> if the current QMC calculation is done and 
     <code>flase</code> otherwise.  This only applies to the QMC calculation
     and not the optimization.
  */
  bool done;

  /**
     Iteration in the QMC calculation this processor is on.  The equilibration
     time is included in the iteration number.
  */
  unsigned long iteration;

  /**
     Storage location for global statistical properties on the root node
     and scratch space on other nodes.
  */
  QMCProperties Properties_total;

  /**
     Local statistical properties calculated during the equilibration phase
     of the calculation.
  */
  QMCProperties equilibrationProperties;

  /**
     An output stream for writing the results of the calculation to.
  */
  ofstream *qmcRslts;

  /**
     An output stream for writing the progress of the calculation to.
  */
  ofstream *qmcOut;

  /**
     Initializes MPI on this processor for a parallel calculation.
  */ 
  void initializeMPI();

  /**
     Initialize the output streams for the calculation.
  */
  void initializeOutputs();

  /**
     Finalize the output streams for the calculation.
  */
  void finalizeOutputs();

  /**
     Writes the checkpoint file for this node to disk.
  */
  void writeCheckpoint();

  /**
     Writes the current statistics for the calculated properties to a file.
     @param label label to indicate which set of transient properties this is.
  */
  void writeTransientProperties(int label);

  /**
     Initializes the state of this object.  If checkpoint files are provided
     and designated to be used, the state is loaded from them; otherwise,
     the state is randomly generated according to the prescription in the
     input data.
  */
  void initializeCalculationState();

  /**
     Send all of the processors the input file name.
     @param argv input from the command line where the first element is the 
     file name.
     @return input file name.
  */
  string sendAllProcessorsInputFileName(char **argv);

  /**
     Send all of the processors a command in the form of an integer.
     @param command command to be sent to all processors.
  */
  void sendAllProcessorsACommand(int command);

  /**
     Gathers the statistics of the calculated properties from all the 
     processors and saves the result in Properties_total.
  */ 
  void gatherProperties();

  /**
     Gathers the statistics for the basis function densities from all the 
     processors and saves the result in Properties_total.
  */
  void gatherDensities();

  /**
     The global results are collected and sent to all processors, where they 
     are used to update the estimated energy, trial energy, and effective time
     step.  
  */
  void synchronizeDMCEnsemble();

  /**
     Checks to see if a command is waiting for this processor.  If there is a
     command waiting, the integer value of the command is returned, if not, 
     -1 is returned.

     @return integer value of the command waiting for this processor or -1 if 
     there is no waiting command.
  */
  int pollForACommand();

  /**
     Writes a summary of the energy statistics to a stream.
     @param strm output stream to write the energy statistics summary to.
  */
  void writeEnergyResultsSummary(ostream & strm);

  /**
     Checks to see if any of the termination criteria for the calculation have
     been satisfied.  
  */
  void checkTerminationCriteria();

  /**
     Checks to see if the total number of QMC steps performed by the
     calculation exceeds the maximum value set in the input.
  */
  void checkMaxStepsTerminationCriteria();

  /**
     Checks to see if the energy standard deviation has converged sufficiently
     to terminate.  The desired convergence is specified in the input.
  */
  void checkConvergenceBasedTerminationCriteria();

  void equilibration_step();

  /**
     Updates the estimated value of the energy.
  */
  void updateEstimatedEnergy(QMCProperties * Properties);

  /**
     Updates the diffusion QMC (DMC) trial energy.
  */
  void updateTrialEnergy(double totalWeight, int nwalkers_original);

  /**
     Updates the effective time step used for diffusion QMC (DMC).
  */
  void updateEffectiveTimeStep(QMCProperties * Properties);

  /**
     Forces all processors to synchronize at this point.
  */
  void synchronizationBarrier();

  /**
     Writes the state of this object to an XML stream.
     @param strm XML stream
  */
  void writeXML(ostream & strm);

  /**
     Loads the state of this object from an XML stream.  The input stream must
     be formatted exactly like the output from <code>toXML</code> because
     it is not intelligent.
     @param strm XML stream
  */
  void readXML(istream & strm);
};

#endif
