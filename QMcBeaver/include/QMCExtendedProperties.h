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

  QMCExtendedProperties();

  QMCProperties * getProperties();

  Stopwatch * getEquilibrationStopwatch();

  Stopwatch * getPropagationStopwatch();

  void zeroOut();

  long getStartingStep();

  void setStartingStep(long i);
};

#endif






