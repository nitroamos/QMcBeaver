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

#include "QMCDerivativeProperties.h"


QMCDerivativeProperties::QMCDerivativeProperties(QMCProperties * properties, 
						 double dt)
{
  this->properties = properties;
  this->dt         = dt;
}

double QMCDerivativeProperties::getEffectiveTimeStep()
{
  double value = dt * properties->distanceMovedAccepted.getAverage()/
    properties->distanceMovedTrial.getAverage();

  return value;
}

double QMCDerivativeProperties::getEffectiveTimeStepVariance()
{
  double distanceMovedAcceptAve =
    properties->distanceMovedAccepted.getAverage();
  double distanceMovedAcceptVar = 
    properties->distanceMovedAccepted.getVariance();

  double distanceMovedTrialAve = properties->distanceMovedTrial.getAverage();
  double distanceMovedTrialVar = properties->distanceMovedTrial.getVariance();

  double temp1 = 1.0/distanceMovedTrialAve/distanceMovedTrialAve;
  double temp2 = distanceMovedAcceptAve * temp1;
  temp2 = temp2*temp2;

  double result = dt * ( temp1 * distanceMovedAcceptVar + 
			 temp2 * distanceMovedTrialVar );

  return result;
}

double QMCDerivativeProperties::getEffectiveTimeStepStandardDeviation()
{
  return sqrt( getEffectiveTimeStepVariance() );
}

double QMCDerivativeProperties::getVirialRatio()
{
  double value = -properties->potentialEnergy.getAverage() /
    properties->kineticEnergy.getAverage();

  return value;
}

double QMCDerivativeProperties::getVirialRatioVariance()
{
  double vave = properties->potentialEnergy.getAverage();
  double vvar = properties->potentialEnergy.getVariance();

  double tave = properties->kineticEnergy.getAverage();
  double tvar = properties->kineticEnergy.getVariance();

  double temp1 = 1.0/tave/tave;
  double temp2 = vave*temp1;
  temp2 = temp2*temp2;

  double result = temp1 * vvar + temp2 * tvar;

  return result;
}

double QMCDerivativeProperties::getVirialRatioStandardDeviation()
{
  return sqrt( getVirialRatioVariance() );
}

ostream& operator <<(ostream& strm, QMCDerivativeProperties &rhs)
{
  strm << "dt_effective: " << rhs.getEffectiveTimeStep() << " +/- " 
       << rhs.getEffectiveTimeStepStandardDeviation() << endl;

  strm << "Virial Ratio (-<V>/<T>): " << rhs.getVirialRatio() << " +/- "
       << rhs.getVirialRatioStandardDeviation() << endl;

  return strm;
}


