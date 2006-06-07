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
						 QMCFutureWalkingProperties * fwProperties,
						 double dt)
{
  this->properties   = properties;
  this->fwProperties = fwProperties;
  this->dt           = dt;
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

double QMCDerivativeProperties::getVirialRatio(int whichFW)
{
  double value = -(fwProperties->props[FW_PE])[whichFW].getAverage() /
    (fwProperties->props[FW_KE])[whichFW].getAverage();

  return value;
}

double QMCDerivativeProperties::getVirialRatioVariance(int whichFW)
{
  double vave = (fwProperties->props[FW_PE])[whichFW].getAverage();
  double vvar = (fwProperties->props[FW_PE])[whichFW].getVariance();

  double tave = (fwProperties->props[FW_KE])[whichFW].getAverage();
  double tvar = (fwProperties->props[FW_KE])[whichFW].getVariance();

  double temp1 = 1.0/tave/tave;
  double temp2 = vave*temp1;
  temp2 = temp2*temp2;

  double result = temp1 * vvar + temp2 * tvar;

  return result;
}

double QMCDerivativeProperties::getVirialRatioStandardDeviation(int whichFW)
{
  return sqrt( getVirialRatioVariance(whichFW) );
}

ostream& operator <<(ostream& strm, QMCDerivativeProperties &rhs)
{
  strm << "dt_effective: " << rhs.getEffectiveTimeStep() << " +/- " 
       << rhs.getEffectiveTimeStepStandardDeviation() << endl;

  strm << endl << "--------------- Virial Ratio (-<V>/<T>) ---------------" << endl;
  for(int fw=0; fw<rhs.fwProperties->numFutureWalking; fw++)
    {
      strm.width(7);
      strm << globalInput.flags.future_walking[fw] << ": " << rhs.getVirialRatio(fw) <<
	" +/- " << rhs.getVirialRatioStandardDeviation(fw) << endl;

      double vave = (rhs.fwProperties->props[FW_PE])[fw].getAverage();
      double vvar = (rhs.fwProperties->props[FW_PE])[fw].getVariance();

      double tave = (rhs.fwProperties->props[FW_KE])[fw].getAverage();
      double tvar = (rhs.fwProperties->props[FW_KE])[fw].getVariance();

      strm << "         Total Energy Estimator (KE/Virial): " << -tave << " +/- " << sqrt(tvar) << endl; 
      strm << "         Total Energy Estimator (PE/Virial): " << vave/2.0 << " +/- " << sqrt(vvar)/2.0 << endl; 
    }
     
  return strm;
}


