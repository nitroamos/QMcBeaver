#ifndef QMCEquilibrationArray_H
#define QMCEquilibrationArray_H

#define EQ 30

#include "QMCExtendedProperties.h"
#include "IeeeMath.h"

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
    @return index of the element with the lowest variance for the energy.
  */
  int getDecorrObjectIndex();

  /**
    Returns a^b.  This function had to be included because some compilers
    don't allow pow to be called with two integer arguments.
    @param a base.
    @param b power.
    @return a^b.
  */
  long power(int a,int b);

  /** 
    True if basis function density is being calculated.
  */
  bool calc_density;

  /**
    Only has a value if calc_density = true
  */
  int nBasisFunc;

  /** 
    True if nuclear forces are being calculated.
  */
  bool calc_forces;
  
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
    Tells the object if basis function densities are being calculated.
  */
  void setCalcDensity(bool calcDensity, int nbasisfunctions);

  /**
    Tells the object if basis function densities are being calculated.
  */
  void setCalcForces(bool calcForces, int dim1, int dim2);
  
  /**
    Adds a new data sample to the live objects in the array and updates
    decorr_objects if necessary.
    @param timeStepProperties new sample.
    @param totalWeight total weight of the new sample.
    @param nWalkers number of walkers that generated the sample.
  */
  void newSample
  (QMCProperties * timeStepProperties, double totalWeight, int nWalkers);

  /**
    Gets the element of the array that has the lowest variance for the total
    energy and updates the current_index.
    @return the element of the array with the lowest variance for the energy.
  */
  QMCProperties * chooseDecorrObject();

  /**
    Gets the propagation stopwatch for the element of the array with the
    lowest variance for the total energy.
    @return the propagation stopwatch for the element of the array with the 
    lowest variance for the energy.
  */
  Stopwatch * getPropagationStopwatch();

  /**
    Gets the equilibration stopwatch for the element of the array with the
    lowest variance for the total energy.
    @return the equilibration stopwatch for the element of the array with the 
    lowest variance for the energy.
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

  /**
    Writes the state of the object to an XML stream.
    @param strm stream to write data to.
  */
  void toXML(ostream& strm);

  /**
    Reads the state of the object from an XML stream.
    @param strm stream to read data from.
  */
  void readXML(istream& strm);

};

#endif
