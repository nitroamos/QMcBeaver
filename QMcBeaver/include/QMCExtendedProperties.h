// Header file for QMCExtendedProperties.

#ifndef QMCExtendedProperties_H
#define QMCExtendedProperties_H

#include "QMCProperties.h"
#include "Stopwatch.h"

using namespace std;

/**
  Includes the data for a calculation as well as timing data as to when it was
  gathered.
*/

class QMCExtendedProperties
{
 private:

  /**
    All of the data for the calculation.
  */
  QMCProperties properties;

  /** 
    A timer for the equilibration phase of the calculation.
  */
  Stopwatch equilibrationStopwatch;

  /**
    A timer for the propagation phase of the calculation.
  */
  Stopwatch propagationStopwatch;

  /**
    The iteration on which this object started to collect statistics.
  */
  long startingStep;

 public:

  /**
    Creates and initializes the object.
  */
  QMCExtendedProperties();

  /**
    Gets the accumulated statistics in this object.
    @return accumulated statistics in this object.
  */
  QMCProperties * getProperties();

  /**
    Gets the stopwatch for the equilibration phase of the calculation.
    @return equilibration stopwatch.
  */
  Stopwatch * getEquilibrationStopwatch();

  /**
    Gets the stopwatch for the propogation phase of the calculation.
    @return propagation stopwatch.
  */
  Stopwatch * getPropagationStopwatch();

  /**
    Zeros out statistics and resets the timers.
  */
  void zeroOut();

  /**
    Gets the iteration on which this object started to collect statistics.
    @return starting step for this object.
  */
  long getStartingStep();

  /**
    Sets the step on which this object started to collect statistics.
    @param starting step for this object.
  */
  void setStartingStep(long i);
};

#endif






