//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef QMCDerivativeProperties_H
#define QMCDerivativeProperties_H

#include "QMCProperties.h"

/**
  All of the calculated quantities and properties that are derived 
  from quantities and properties evaluated during a calculation.
  */

class QMCDerivativeProperties
{
public:
  /**
    Creates and initializes an instance of this class.

    @param properties calculated properties for the system.
    @param dt time step for the calculation.
    */
  QMCDerivativeProperties(QMCProperties * properties, double dt);

  /**
    Gets the effective time step for the calculation.

    @return effective time step for the calculation.
    */
  double getEffectiveTimeStep();

  /**
    Gets the variance of the calculated effective time step for the 
    calculation.

    @return variance of the effective time step for the calculation.
    */
  double getEffectiveTimeStepVariance();

  /**
    Gets the standard deviation of the calculated effective time step for the 
    calculation.

    @return standard deviation of the effective time step for the calculation.
    */
  double getEffectiveTimeStepStandardDeviation();

  /**
    Gets the virial ratio for the calculation.  The virial ratio is 
    \f$-\left<V\right>/\left<T\right>\f$ where \f$\left<V\right>\f$ is the
    expectation value of the potential energy and \f$\left<T\right>\f$ is
    the expectation value of the kinetic energy.

    @return virial ratio.
    */
  double getVirialRatio();

  /**
    Gets the variance of the calculated virial ratio for the calculation.

    @return variance of the virial ratio.
    */
  double getVirialRatioVariance();

  /**
    Gets the standard deviation of the calculated virial ratio for the 
    calculation.

    @return standard deviation of the virial ratio.
    */
  double getVirialRatioStandardDeviation();

  /**
    Formats and prints the properties to a stream in human readable fromat.
    */
  friend ostream& operator <<(ostream& strm, QMCDerivativeProperties &rhs);

private:
  QMCProperties * properties;
  double dt;
};

#endif
