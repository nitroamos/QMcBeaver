//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2002.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef QMCStopwatches_H
#define QMCStopwatches_H

#include "Stopwatch.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

using namespace std;

/**
  A collection of Stopwatch objects used to record information
  relevant to the timing of a QMC calculation.
*/

class QMCStopwatches
{
 private:
  Stopwatch Initialization;
  Stopwatch Equilibration;
  Stopwatch Propagation;
  Stopwatch Communication_send;
  Stopwatch Communication_reduce;
  Stopwatch Communication_synch;
  Stopwatch Communication_poll;
  Stopwatch Optimization;
  Stopwatch Total;

 public:

  /**
    Creates a new instance of this class with all timers stopped.
  */

  QMCStopwatches();

  /**
    Stops all stopwatches in this object which are running.
  */

  void stop();

  /**
    Resets all stopwatches in this object and leaves the stopwatches stopped.
  */

  void reset();

  /**
    Gets the stopwatch which times the initialization of the calculation.
  */

  Stopwatch * getInitializationStopwatch();

  /**
    Gets the stopwatch that times the equilibration phase of the calc.
  */

  Stopwatch * getEquilibrationStopwatch();

  /**
    Gets the stopwatch which times the useful propagation of walkers.
    The time required to initialize the walkers is not included.

    @return the stopwatch.
  */

  Stopwatch * getPropagationStopwatch();

  /**
    Gets the stopwatch which times the sending of commands between
    processors.

    @return the stopwatch.
  */

  Stopwatch * getSendCommandStopwatch();

  /**
    Gets the stopwatch which times the gathering of QMCProperties from
    all processors.

    @return the stopwatch.
  */

  Stopwatch * getGatherPropertiesStopwatch();

  /**
    Gets the stopwatch which times the synchronization of all the processors.

    @return the stopwatch.
  */

  Stopwatch * getCommunicationSynchronizationStopwatch();

  /**
    Gets the stopwatch which times how long is devoted to seeing if a
    processor has a command waiting for it.

    @return the stopwatch.
  */

  Stopwatch * getCommandPollingStopwatch();

  /**
    Gets the stopwatch which times the VMC optimization.

    @return the stopwatch.
  */

  Stopwatch * getOptimizationStopwatch();

  /**
    Gets the stopwatch which records the total time of the calculation.

    @return the stopwatch.
  */

  Stopwatch * getTotalTimeStopwatch();

  /**
     Sets two object equal.
  */

  void operator = ( const QMCStopwatches &rhs);

  /**
    Returns a QMCStopwatches which is the sum of two QMCStopwatches objects.
  */

  QMCStopwatches operator+(QMCStopwatches & rhs);

  /**
    Writes the timing results of this class to a human readable
    stream.
  */ 

  friend ostream& operator <<(ostream& strm, QMCStopwatches & rhs);

#ifdef PARALLEL

private:

  /**
    A flag which tells if MPI_TYPE has been generated.
  */

  static bool mpiTypeCreated;
  
  /**
    Build MPI_TYPE.
  */

  static void buildMpiType();

  /**
    Build MPI_REDUCE.
  */

  static void buildMpiReduce();

  /**
    An MPI function which allows MPI_Reduce to be used in adding Stopwatches.
  */

  static void Reduce_Function(QMCStopwatches *in, QMCStopwatches *inout, 
                              int *len, MPI_Datatype *dptr);

public:

  /**
    The MPI data type for a QMCStopwatches.
  */

  static MPI_Datatype MPI_TYPE;

  /** 
    The MPI operation for performing MPI_Reduce on QMCStopwatches objects.
  */
  
  static MPI_Op MPI_REDUCE;

#endif


};

#endif





