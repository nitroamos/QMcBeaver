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

double QMCDerivativeProperties::getParameterValue()
{
  if(globalInput.flags.optimize_Psi_criteria == "generalized_eigenvector" ||
     globalInput.flags.optimize_Psi_criteria == "energy_average")
    {
      return properties->energy.getAverage();
    }
  return properties->energy.getSeriallyCorrelatedVariance();
}

Array1D<double> QMCDerivativeProperties::getParameterGradient()
{
  Array1D<double> gradient;
  if(fwProperties->der.dim1() == 0 ||
     fwProperties->der.dim2() == 0)
    {
      //clog << "Error: the parameter gradients are not available in this QMCFwProperties object.\n";
      return gradient;
    }

  int numCI = fwProperties->der.dim1();

  gradient.allocate(numCI);
  for(int ci=0; ci<numCI; ci++)
    {
      double pi    = fwProperties->der(ci,0).getAverage();  // < \frac{ \Psi_i }{ \Psi } >
      double pi_e  = fwProperties->der(ci,1).getAverage();  // < \frac{ \Psi_i }{ \Psi } E_L>
      double pi_e2 = fwProperties->der(ci,2).getAverage();  // < \frac{ \Psi_i }{ \Psi } E_L^2>
      double ei    = fwProperties->der(ci,3).getAverage();  // < E_{L,i} > (= 0 in limit)
      double ei_e  = fwProperties->der(ci,4).getAverage();  // < E_{L,i} E_L >
      double e     = properties->energy.getAverage();       // < E_L >
      double e2    = properties->energy2.getAverage();      // < E_L^2 >

      //References are from:
      //PRL 94, 150201 (2005)
      double der;

      if(globalInput.flags.optimize_Psi_criteria == "generalized_eigenvector") 
	{
	  //Energy Minimization

	  //Formula 6
	  der = pi_e - pi*e;

	} else {
	  //Variance Minimization

	  if(!true)
	    {
	      //Formula 3
	      pi    = 0;
	      pi_e  = 0;
	      pi_e2 = 0;
	    }
	  
	  // < E_{L,i} ( E_L - < E_L> ) >
	  der  = ei_e - ei*e;
	  
	  // < \frac{ \Psi_i }{ \Psi } E_L^2>
	  der += pi_e2;
	  
	  // < \frac{ \Psi_i }{ \Psi } > < E_L^2 >
	  der -= pi*e2;
	  
	  // 2 < E_L > < \frac{ \Psi_i }{ \Psi } ( E_L - < E_L> ) > 
	  der -= 2.0*e*(pi_e - pi*e);
	}
	  
      der *= 2.0;
      gradient(ci) = der;
    }
  return gradient;
}

Array2D<double> QMCDerivativeProperties::getParameterHessian()
{
  Array2D<double> hessian;
  if(fwProperties->hess.dim1() == 0 ||
     fwProperties->hess(0).dim1() == 0 ||
     fwProperties->hess(0).dim2() == 0)
    {
      //clog << "Error: the parameter hessian is not available in this QMCFwProperties object.\n";
      return hessian;
    }

  if(globalInput.flags.optimize_Psi_criteria != "analytical_energy_variance")
    return hessian;

  int numAI = fwProperties->der.dim1();
  hessian.allocate(numAI,numAI);

  for(int ai=0; ai<numAI; ai++)
    {
      for(int aj=0; aj<=ai; aj++)
	{
	  // < E_{L,i} E_{L,j} >
	  double h1 = (fwProperties->hess(0))(ai,aj).getAverage();
	  // < E_{L,i} > < E_{L,j} >
	  double h2 = fwProperties->der(ai,3).getAverage()
	            * fwProperties->der(aj,3).getAverage();

	  double h3 = 2.0 * (h1 - h2);
	  hessian(ai,aj) = h3;
	  hessian(aj,ai) = h3;
	}
    }

  return hessian;
}

Array2D<double> QMCDerivativeProperties::getParameterHamiltonian()
{
  Array2D<double> hamiltonian;
  if(fwProperties->hess.dim1() == 0 ||
     fwProperties->hess(0).dim1() == 0 ||
     fwProperties->hess(0).dim2() == 0)
    {
      clog << "Error: the parameter hamiltonian is not available in this QMCFwProperties object.\n";
      return hamiltonian;
    }

  if(globalInput.flags.optimize_Psi_criteria != "generalized_eigenvector")
    return hamiltonian;

  int numAI = fwProperties->der.dim1();
  hamiltonian.allocate(numAI+1,numAI+1);
  hamiltonian = 0.0;

  // < E_L >
  double e   = properties->energy.getAverage();

  hamiltonian(0,0) = e;

  for(int ai=0; ai<numAI; ai++)
    {
      // < \frac{ \Psi_i }{ \Psi } E_L >
      double pi_e = fwProperties->der(ai,1).getAverage();
      // < \frac{ \Psi_i }{ \Psi } >
      double pi  = fwProperties->der(ai,0).getAverage();
      // < E_{L,i} > (= 0 in limit)
      double ei  = fwProperties->der(ai,3).getAverage();
      
      hamiltonian(ai+1,0) = pi_e - pi*e;
      hamiltonian(0,ai+1) = pi_e - pi*e + ei;

      for(int aj=0; aj<numAI; aj++)
	{
	  // < \frac{ \Psi_i }{ \Psi } \frac{ \Psi_j }{ \Psi } E_L >
	  double ppe  = (fwProperties->hess(0))(ai,aj).getAverage();

	  // < \frac{ \Psi_i }{ \Psi } E_{L,j} >
	  double pi_ej = (fwProperties->hess(2))(ai,aj).getAverage();

	  // < \frac{ \Psi_j }{ \Psi } E_L >
	  double pj_e = fwProperties->der(aj,1).getAverage();

	  // < \frac{ \Psi_j }{ \Psi } >
	  double pj  = fwProperties->der(aj,0).getAverage();

	  // < E_{L,j} > (= 0 in limit) 
	  double ej  = fwProperties->der(aj,3).getAverage();
	  
	  double val = ppe - pi*pj_e - pj*pi_e + pi*pj*e;
	  val += pi_ej - pi*ej;

	  hamiltonian(ai+1,aj+1) = val;
	}
    }

  return hamiltonian;
}

Array2D<double> QMCDerivativeProperties::getParameterOverlap()
{
  Array2D<double> overlap;
  if(fwProperties->hess.dim1() == 0 ||
     fwProperties->hess(0).dim1() == 0 ||
     fwProperties->hess(0).dim2() == 0)
    {
      //clog << "Error: the parameter overlap is not available in this QMCFwProperties object.\n";
      return overlap;
    }
  
  if(globalInput.flags.optimize_Psi_criteria != "generalized_eigenvector")
    return overlap;

  int numAI = fwProperties->der.dim1();
  overlap.allocate(numAI+1,numAI+1);
  overlap = 0.0;

  overlap(0,0) = 1.0;

  for(int ai=0; ai<numAI; ai++)
    {
      overlap(0,ai+1) = 0.0;
      overlap(ai+1,0) = 0.0;

      for(int aj=0; aj<=ai; aj++)
	{
	  // < \frac{ \Psi_i }{ \Psi } \frac{ \Psi_j }{ \Psi } >
	  double pi_pj = (fwProperties->hess(1))(ai,aj).getAverage();

	  // < \frac{ \Psi_i }{ \Psi } >
	  double pi  = fwProperties->der(ai,0).getAverage();
	  // < \frac{ \Psi_j }{ \Psi } >
	  double pj  = fwProperties->der(aj,0).getAverage();

	  double val = pi_pj - pi*pj;
	  overlap(ai+1,aj+1) = val;
	  overlap(aj+1,ai+1) = val;
	}
    }

  return overlap;
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


