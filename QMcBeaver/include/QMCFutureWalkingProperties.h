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

#ifndef QMCFWPROPERTIES_H
#define QMCFWPROPERTIES_H

#include <iostream>
#include <string>
#include <math.h>

#include "QMCProperty.h"
#include "QMCInput.h"
#include "Array1D.h"

using namespace std;

/**
  All of the quantities and properties evaluated during a calculation.
*/

class QMCFutureWalkingProperties
{
public:

  /**
    this will always be at least 1 with the smallest element
    0 representing future walking turned off
  */
  static int numFutureWalking;

  Array1D< QMCProperty > fwEnergy;
  Array1D< QMCProperty > fwKineticEnergy;
  Array1D< QMCProperty > fwPotentialEnergy;  
  Array1D< QMCProperty > r12;  
  Array1D< QMCProperty > r2;  

  /**
    Tells if nuclear forces are being calculated.
  */
  bool calc_forces;

  /**
     Forces for each atom.
  */
  Array1D< Array2D<QMCProperty> > nuclearForces;
   
  /**
    Creates an instance of the class.
  */
  QMCFutureWalkingProperties();
  ~QMCFutureWalkingProperties();
  
  /**
    Sets all of the data in the object to zero.
  */
  void zeroOut();

  void matchParametersTo( const QMCFutureWalkingProperties &rhs );

  /**
    Sets two objects equal.
  */
  void operator = ( const QMCFutureWalkingProperties &rhs );

  /**
    Returns the sum of two QMCFutureWalkingProperties.
    @return sum of two QMCFutureWalkingProperties
  */
  QMCFutureWalkingProperties operator + ( QMCFutureWalkingProperties &rhs );

  /**
    Adds the statistics calculated for a time step to the object as new 
    samples.  This is not the same as adding the two QMCFutureWalkingProperties objects 
    together.
    @param newProperties the statistics calculated at the time step.
    @param weight the weight of the new samples.
    @param nwalkers the number of walkers used to make this sample.
  */
  void newSample(QMCFutureWalkingProperties* newProperties, double weight, int nwalkers);

  /**
    Tells the object if basis function densities are being calculated.

    @param calcDensity- true if densities are being calculated.
    @param nbasisfunctions- number of basis functions.
  */
  void setCalcForces(bool calcForces, int dim1, int dim2);

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
  friend ostream& operator <<(ostream& strm, QMCFutureWalkingProperties &rhs);
};

#endif









