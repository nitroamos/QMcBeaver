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

Array1D<double> QMCDerivativeProperties::getParameterGradient()
{
  Array1D<double> gradient;
  if(properties->der.dim1() == 0 ||
     properties->der.dim2() == 0)
    {
      //clog << "Error: the parameter gradients are not available in this QMCProperties object.\n";
      return gradient;
    }

  int numCI = properties->der.dim1();

  gradient.allocate(numCI);
  for(int ci=0; ci<numCI; ci++)
    {
      double t1 = properties->der(ci,0).getAverage();  // < \frac{ \Psi_i }{ \Psi } >
      double t2 = properties->der(ci,1).getAverage();  // < \frac{ \Psi_i }{ \Psi } E_L>
      double t3 = properties->der(ci,2).getAverage();  // < \frac{ \Psi_i }{ \Psi } E_L^2>
      double t4 = properties->der(ci,3).getAverage();  // < E_{L,i} > (= 0 in limit)
      double t5 = properties->der(ci,4).getAverage();  // < E_{L,i} E_L >
      double t6 = properties->energy.getAverage();     // < E_L >
      double t7 = properties->energy2.getAverage();    // < E_L^2 >
      
      if(!true)
	{
	  t1 = 0;
	  t2 = 0;
	  t3 = 0;
	}
      double der;

      //Looking at the PRL 94, 150201 (2005) paper, formula 2
      // < E_{L,i} ( E_L - < E_L> ) >
      der  = t5 - t4*t6;

      // < \frac{ \Psi_i }{ \Psi } E_L^2>
      der += t3;

      // < \frac{ \Psi_i }{ \Psi } > < E_L^2 >
      der -= t1*t7;

      // 2 < E_L > < \frac{ \Psi_i }{ \Psi } ( E_L - < E_L> ) > 
      der -= 2.0*t6*(t2 - t1*t6);

      der *= 2.0;

      gradient(ci) = der;
    }
  return gradient;
}

Array2D<double> QMCDerivativeProperties::getParameterHessian()
{
  Array2D<double> hessian;
  if(properties->hess.dim1() == 0 ||
     properties->hess.dim2() == 0)
    {
      //clog << "Error: the parameter hessian is not available in this QMCProperties object.\n";
      return hessian;
    }
  
  int numCI = properties->der.dim1();
  hessian.allocate(numCI,numCI);
  
  for(int ci=0; ci<numCI; ci++)
    {
      for(int cj=0; cj<numCI; cj++)
	{
	  // < E_{L,i} E_{L,j} >
	  double h1 = properties->hess(ci,cj).getAverage();
	  // < E_{L,i} > < E_{L,j} >
	  double h2 = properties->der(ci,3).getAverage()
	            * properties->der(cj,3).getAverage();
	  
	  hessian(ci,cj) = 2.0 * (h1 - h2);
	}
    }
  return hessian;
}

ostream& operator <<(ostream& strm, QMCDerivativeProperties &rhs)
{
  strm << "dt_effective: " << rhs.getEffectiveTimeStep() << " +/- " 
       << rhs.getEffectiveTimeStepStandardDeviation() << endl;

  strm << endl << "--------------- Virial Ratio (-<V>/<T>) ---------------" << endl;
  for(int fw=0; fw<rhs.fwProperties->numFutureWalking; fw++)
    {
      if((rhs.fwProperties->props[FW_PE])[fw].getNumberSamples() <= 0) break;
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


