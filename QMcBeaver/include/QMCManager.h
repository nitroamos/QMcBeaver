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

//The old macros conflicted with some Windows macros
#define QMC_WORK_STEP   0
#define QMC_REDUCE      1
#define QMC_TERMINATE   2
#define QMC_SYNCHRONIZE 3
#define QMC_REDUCE_ALL  4

//These are all the signals that QMCManager has been programmed
//to handle.
enum signalType {SIG_REDUCE, SIG_INCREASE, SIG_QUIT, SIG_NOTHING};

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <signal.h>

#include "QMCRun.h"
#include "QMCDerivativeProperties.h"
#include "QMCCorrelatedSamplingVMCOptimization.h"
#include "QMCCopyright.h"
#include "QMCProperties.h"
#include "QMCProperty.h"
#include "QMCPropertyArrays.h"
#include "QMCNuclearForces.h"

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

     @return whether the run was cut short (e.g. kill signal)
  */
  bool run(bool equilibrate);

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
     Useful for resetting the timers for every optimization step.
  */
  void resetTimers()
    {
      localTimers.reset();
    }

  /**
     Writes the restart file for the calculation.
     @param filename the name of the file to write the input into
  */
  void writeRestart(string filename);

  /**
     Writes the restart file for the calculation.
     Uses the Input.flags.restart_file_name as
     filename.
  */
  void writeRestart();

  /**
     Writes the basis function density for the calculation.
  */
  void writeBFDensity();

  /**
     Writes the calculated forces for the calculation.
  */
  void writeForces();

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

  /**
     Call this function when you want QMCManager to change its
     signal flags.
  */
  static void receiveSignal(signalType signal);

 private:
 
  /**
     These are signal flags. Since system signaling can only
     call static functions, we translate any signal information
     into these flags, which QMCManager can look at during some
     decision.
  */
  static bool SIGNAL_SAYS_QUIT;
  static bool REDUCE_ALL_NOW;
  static bool INCREASE_TIME;
  static bool PRINT_SIG_INFO;

  /**
     Walkers plus the basic functionality needed to build a QMC calculation.
  */
  QMCRun QMCnode;

  /**
     Input data for the calculation.

  QMCInput Input;
  */

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
     Same purpose as Properties_total, except that this one holds
     the properties collected as future walking properties. It turned
     out to be too annoying to include dynamic arrays in QMCProperties.
     This permits the amount of data in QMCProperties to be known at
     compile time while permitting fwProperties_total to depend on
     information from the input file. This difference is important
     for how MPI reduces the properties.
  */
  QMCPropertyArrays fwProperties_total;

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
     Storage location for the global parallel spin histogram on the root node 
     and scratch space on other nodes.
  */
  Array1D<double> pllSpinHistogram_total;

  /**
     Storage location for the global opposite spin histogram on the root node
     and scratch space on other nodes.
  */
  Array1D<double> oppSpinHistogram_total;

  /**
     Storage location for the global alpha one electron histograms on the root
     node and scratch space on other nodes.
  */
  Array1D< Array1D<double> > alphaHistograms_total;

  /**
     Storage location for the global beta one electron histograms on the root
     node and scratch space on other nodes.
  */
  Array1D< Array1D<double> > betaHistograms_total;

  /**
     Storage location for the global parallel pair x coordinate correlation
     diagram on the root node and scratch space on the other nodes.
  */
  Array1D< Array1D<double> > pllxCorrelationDiagram_total;

  /**
     Storage location for the global parallel pair y coordinate correlation
     diagram on the root node and scratch space on the other nodes.
  */
  Array1D< Array1D<double> > pllyCorrelationDiagram_total;

  /**
     Storage location for the global parallel pair z coordinate correlation
     diagram on the root node and scratch space on the other nodes.
  */
  Array1D< Array1D<double> > pllzCorrelationDiagram_total;

  /**
     Storage location for the global opposite pair x coordinate correlation
     diagram on the root node and scratch space on the other nodes.
  */
  Array1D< Array1D<double> > oppxCorrelationDiagram_total;

  /**
     Storage location for the global opposite pair y coordinate correlation
     diagram on the root node and scratch space on the other nodes.
  */
  Array1D< Array1D<double> > oppyCorrelationDiagram_total;

  /**
     Storage location for the global opposite pair z coordinate correlation
     diagram on the root node and scratch space on the other nodes.
  */
  Array1D< Array1D<double> > oppzCorrelationDiagram_total;

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
     @param iseed the value read in from the input file
  */
  void initializeCalculationState(long int iseed);

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
     Gathers the statistics of the calculated properties from all the 
     processors.

     These are the non-critical properties, no decisions are made based
     on this data, so we shouldn't collect the data so often.
  */ 
  void gatherExtraProperties();

  /**
     Gathers the statistics for the basis function densities from all the 
     processors and saves the result in Properties_total.
  */
  void gatherDensities();

  /**
     Gathers the statistics for the forces from all the 
     processors and saves the result in fwProperties_total.
  */
  void gatherForces();

  /**
     Gathers the electron density histograms from all the processors.
  */
  void gatherHistograms();

  /**
     Writes out the electron density histograms to files.
  */
  void writeElectronDensityHistograms();

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
  void writeEnergyResultsHeader(ostream & strm);

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
     Checks if we have run the calculation for a long enough
     period of time, as specified by max_time in the input.
  */
  void checkMaxTimeTerminationCriteria();

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
     @return whether the read was successful
  */
  bool readXML(istream & strm);
};

#endif
