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

#ifndef QMCPROPERTIES_H
#define QMCPROPERTIES_H

#include <iostream>
#include <string>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "QMCProperty.h"
#include "Array1D.h"

using namespace std;

/**
  All of the quantities and properties evaluated during a calculation.
*/

class QMCProperties
{
public:

  /**
    Total energy of the system.
  */
  QMCProperty energy;

  /**
    Kinetic energy of the system.
  */
  QMCProperty kineticEnergy;

  /**
    Potential energy of the system.
  */
  QMCProperty potentialEnergy;

  /**
    Log of the weights on the walkers.
  */
  QMCProperty logWeights;

  /**
    Probability a trial move is accepted.
  */
  QMCProperty acceptanceProbability;

  /**
    Average distance an accepted move travels.
  */
  QMCProperty distanceMovedAccepted;

  /**
    Average distance for a trial move.
  */
  QMCProperty distanceMovedTrial;

  /**
    Tells if basis function density is being calculated.
  */
  bool calc_density;

  /**
    The number of basis functions.  Only has a value if calc_density is true.
  */
  int nBasisFunc;

  /**
     Densities for the basis functions.
  */
  Array1D<QMCProperty> chiDensity;

  /**
    Creates an instance of the class.
  */
  QMCProperties();

  /**
    Sets all of the data in the object to zero.
  */
  void zeroOut();

  /**
    Sets two objects equal.
  */
  void operator = ( const QMCProperties &rhs);

  /**
    Returns the sum of two QMCProperties.
    @return sum of two QMCProperties
  */
  QMCProperties operator + ( QMCProperties &rhs );

  /**
    Tells the object if basis function densities are being calculated.

    @param calcDensity- true if densities are being calculated.
    @param nbasisfunctions- number of basis functions.
  */
  void setCalcDensity(bool calcDensity, int nbasisfunctions);

  /**
    Writes the state of this object to an XML stream.
    @param strm XML stream
  */
  void toXML(ostream& strm);

  /**
    Loads the state of this object from an XML stream.
    @param strm XML stream
  */
  void readXML(istream& strm);

  /**
    Formats and prints the properties to a stream.
  */
  friend ostream& operator <<(ostream& strm, QMCProperties &rhs);

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
    Build MPI_REDUCE.  This must be changed when a property is added to or 
    removed from this class.
  */
  static void buildMpiReduce();

  /**
    An MPI function which allows MPI_Reduce to be used in adding QMCProperties.
  */
  static void Reduce_Function(QMCProperties *in, QMCProperties *inout, 
                              int *len, MPI_Datatype *dptr);

public:

  /**
    The MPI data type for a QMCProperties.
  */
  static MPI_Datatype MPI_TYPE;

  /** 
    The MPI operation for performing MPI_Reduce on QMCProperties.
  */
  static MPI_Op MPI_REDUCE;

#endif

};

#endif








