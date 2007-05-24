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
#include "QMCFutureWalkingProperties.h"

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
    @param fwProperties the properties calculated with future walking
    @param dt time step for the calculation.
    */
  QMCDerivativeProperties(QMCProperties * properties,
			  QMCFutureWalkingProperties * fwProperties,
			  double dt);

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

	 @param whichFW which future walking index we want to query
    @return virial ratio.
    */
  double getVirialRatio(int whichFW);

  /**
    Gets the variance of the calculated virial ratio for the calculation.

	 @param whichFW which future walking index we want to query
    @return variance of the virial ratio.
    */
  double getVirialRatioVariance(int whichFW);

  /**
    Gets the standard deviation of the calculated virial ratio for the 
    calculation.

	 @param whichFW which future walking index we want to query
    @return standard deviation of the virial ratio.
    */
  double getVirialRatioStandardDeviation(int whichFW);

  /**
     Convert the der terms collected in QMCProperties
     into the actual gradient.
  */
  Array1D<double> getParameterGradient();

  /**
     Convert the hess terms collected in QMCProperties
     into the actual hessian.
  */
  Array2D<double> getParameterHessian();

  /**
    Formats and prints the properties to a stream in human readable fromat.
    */
  friend ostream& operator <<(ostream& strm, QMCDerivativeProperties &rhs);

private:
  QMCProperties * properties;
  QMCFutureWalkingProperties * fwProperties;
  double dt;
};

#endif
