#ifndef QMCEquilibrationArray_H
#define QMCEquilibrationArray_H

#define EQ 30

#include "QMCExtendedProperties.h"
#include <math.h>

using namespace std;

/**
  This class contains an array of QMCExtendedProprties objects, where the ith
  element starts collecting statistics on the (2^i)th step.  The object with 
  the lowest variance in the energy is the one with the optimal number of 
  equilibration steps.
*/

class QMCEquilibrationArray
{
 private:

  /** 
    The array of QMCExtendedProperties objects.
  */
  QMCExtendedProperties Eq_Array[EQ];

  /**
    The number of active objects in the array.
  */
  int decorr_objects;

  /**
    Returns the index of the element of the array with the lowest variance for
    the total energy.
  */
  int getDecorrObjectIndex();

  /**
    Returns a^b.  This function had to be included because some compilers
    don't allow pow to be called with two integer arguments.
  */

  long power(int a,int b);

 public:

  /**
    Creates a zeroed out instance of the class.
  */
  QMCEquilibrationArray();

  /**
    Sets all of the data in the object to zero.
  */
  void zeroOut(); 

  /**
    Adds a new data sample to the live objects in the array and updates
    decorr_objects if necessary. 
  */
  void newSample
  (QMCProperties * timeStepProperties, double totalWeight, int nWalkers);

  /**
    Returns the element of the array that has the lowest variance for the total
    energy and updates the current_index.
  */
  QMCProperties * chooseDecorrObject();

  /**
    Returns the propagation stopwatch for the element of the array with the
    lowest variance for the total energy.
  */
  Stopwatch * getPropagationStopwatch();

  /**
    Returns the equilibration stopwatch for the element of the array with the
    lowest variance for the total energy.
  */
  Stopwatch * getEquilibrationStopwatch();

  /**
    Starts the propagation stopwatches in the active elements of the array and 
    the equilibration stopwatches in the inactive elements.
  */
  void startTimers();

  /**
    Stops the timers.
  */
  void stopTimers();

  void toXML(ostream& strm);

  void readXML(istream& strm);

};

#endif
