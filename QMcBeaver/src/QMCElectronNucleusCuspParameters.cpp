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

#include "QMCElectronNucleusCuspParameters.h"

QMCElectronNucleusCuspParameters::QMCElectronNucleusCuspParameters()
{
  rc = -2.0;
  sgn0 = -2;
  sigma_sq = -1.0;
  Zeff = -1.0;
  Z = -1;
  C = 0.0;
  noC = true;
}

void QMCElectronNucleusCuspParameters::operator=(const QMCElectronNucleusCuspParameters& rhs)
{
  rc = rhs.rc;
  alpha = rhs.alpha;
  idealCurve = rhs.idealCurve;
  orbitalCoefficients = rhs.orbitalCoefficients;

  C = rhs.C;
  noC = rhs.noC;

  sgn0 = rhs.sgn0;
  sigma_sq = rhs.sigma_sq;

  phi0 = rhs.phi0;
  n0 = rhs.n0;
  
  Z = rhs.Z;
  Zeff = rhs.Zeff;
}

void QMCElectronNucleusCuspParameters::initialize
            (const Array2D<qmcfloat>& temp_OrbitalCoefficients, double temp_n0,
	     double temp_rc, const Polynomial& temp_idealCurve, int temp_Z)
{
  orbitalCoefficients = temp_OrbitalCoefficients;
  idealCurve = temp_idealCurve;

  if (temp_Z > 0)
    Z = temp_Z;
  else
    {
      cerr << "ERROR in QMCElectronNucleusCuspParameters: temp_Z = " << temp_Z;
      cerr << endl;
      exit(0);
    }

  n0 = temp_n0;

  if (temp_rc > 0)
    rc = temp_rc;
  else
    {
      cerr << "ERROR in QMCElectronNucleusCuspParameters: temp_rc = ";
      cerr << temp_rc << endl;
      exit(0);
    }
}

double QMCElectronNucleusCuspParameters::get_rc()
{
  if (rc < 0.0)
    {
      cerr << "ERROR in QMCElectronNucleusCuspParameters: rc = " << rc << endl;
      exit(0);
    }
  return rc;
}

double QMCElectronNucleusCuspParameters::getSigmaSq()
{
  if (sigma_sq < 0)
    {
      cerr << "ERROR in QMCElectronNucleusCuspParameters: sigma_sq = ";
      cerr << sigma_sq << endl;
      exit(0);
    }
  return sigma_sq;
}

void QMCElectronNucleusCuspParameters::set_rc(double temp_rc)
{
  if (temp_rc >= 0.0)
    rc = temp_rc;
  else
    {
      cerr << "ERROR in QMCElectronNucleusCuspParameters: temp_rc = ";
      cerr << temp_rc << endl;
      exit(0);
    }
}

void QMCElectronNucleusCuspParameters::replaceOrbitalValues(double x, double y,
   double z, double r, double& orb_value, double& orb_gradx, double& orb_grady,
                                      double& orb_gradz, double& orb_laplacian)
{
  if (r > rc)
    {
      // do nothing
    }

  else if (r < rc)
    {
      // First we evaluate the part of the original orbital due to the s type 
      // Gaussians centered on this nucleus.

      double temp_value;
      double temp_gradx;
      double temp_grady;
      double temp_gradz;
      double temp_laplacian;
  
      evaluateOriginalOrbital(x,y,z,r,temp_value,temp_gradx,temp_grady,
			      temp_gradz,temp_laplacian);

      // Then we subtract these values and replace them with those of the new 
      // orbital.

      orb_value -= temp_value;
      orb_gradx -= temp_gradx;
      orb_grady -= temp_grady;
      orb_gradz -= temp_gradz;
      orb_laplacian -= temp_laplacian;

      evaluateReplacementOrbital(x,y,z,r,temp_value,temp_gradx,temp_grady,
				 temp_gradz,temp_laplacian);

      orb_value += temp_value;
      orb_gradx += temp_gradx;
      orb_grady += temp_grady;
      orb_gradz += temp_gradz;
      orb_laplacian += temp_laplacian;
    }
}

void QMCElectronNucleusCuspParameters::evaluateOriginalOrbital(double x, 
          double y, double z, double r, double& orig_value, double& orig_gradx,
                double& orig_grady, double& orig_gradz, double& orig_laplacian)
{
  // This function evaluates the value, gradient, and laplacian of the part of
  // the original orbital due to s type Gaussian basis functions centered on 
  // this nucleus.

  double r_sq = r*r;

  int nGaussians = orbitalCoefficients.dim1();

  double term_value = 0.0;
  double p0 = 0.0;
  double p1 = 0.0;
  double temp = 0.0;

  orig_value = 0.0;
  orig_gradx = 0.0;
  orig_grady = 0.0;
  orig_gradz = 0.0;
  orig_laplacian = 0.0;

  for (int i=0; i<nGaussians; i++)
    {
      p0 = orbitalCoefficients(i,0);
      p1 = orbitalCoefficients(i,1);
      term_value = p1*exp(-p0*r_sq);

      orig_value += term_value;

      temp = -2*p0*term_value;

      orig_gradx += x*temp;
      orig_grady += y*temp;
      orig_gradz += z*temp;

      orig_laplacian += 2*p0*(2*p0*r_sq-3)*term_value;
    }
}

void QMCElectronNucleusCuspParameters::evaluateReplacementOrbital(double x,
            double y, double z, double r, double& rep_value, double& rep_gradx,
                   double& rep_grady, double& rep_gradz, double& rep_laplacian)
{
  // This function evaluates the value, gradient, and laplacian of the 
  // replacement orbital.

  double exp_value = 0.0;
  double temp_gradient = 0.0;

  alpha.evaluate(r);

  exp_value = sgn0*exp(alpha.getFunctionValue());

  rep_value = C + exp_value;

  temp_gradient = (exp_value*alpha.getFirstDerivativeValue())/r;

  rep_gradx = temp_gradient*x;
  rep_grady = temp_gradient*y;
  rep_gradz = temp_gradient*z;

  rep_laplacian = exp_value*(alpha.getSecondDerivativeValue() + 
      alpha.getFirstDerivativeValue()*(2.0/r+alpha.getFirstDerivativeValue()));
}
      
double QMCElectronNucleusCuspParameters::calculateLocalEnergy(double r, 
							      bool rIsZero)
{
  if (rIsZero == true)
    {
      alpha.evaluate(0.0);
      if (noC == true)
	return -0.5*(alpha.getSecondDerivativeValue() + 
	      alpha.getFirstDerivativeValue()*alpha.getFirstDerivativeValue());
      else
	{
	  double exp_value = sgn0*exp(alpha.getFunctionValue());
	  return (-0.5*exp_value/(C+exp_value))*
	    (alpha.getSecondDerivativeValue() + 
	     alpha.getFirstDerivativeValue()*alpha.getFirstDerivativeValue());
	}
    }

  else
    {
      if (r < 0)
	{
	  cerr << "ERROR in ";
	  cerr << "QMCElectronNucleusCuspParameters::calculateLocalEnergy: ";
	  cerr << "r = " << r << endl;
	  exit(0);
	}

      alpha.evaluate(r);
      
      double localKE = -0.5*(alpha.getSecondDerivativeValue() + 
      alpha.getFirstDerivativeValue()*(2.0/r+alpha.getFirstDerivativeValue()));

      if (noC == true)
	return localKE - Zeff/r;

      else
	{
	  double exp_value = sgn0*exp(alpha.getFunctionValue());
	  return localKE*exp_value/(C+exp_value) - Zeff/r;
	}
    }
}
  
void QMCElectronNucleusCuspParameters::fitReplacementOrbital(double temp_phi0)
{
  // This function fits the exponential replacement orbital to the original 
  // orbital at rc.

  phi0 = temp_phi0;

  double rc_sq = rc*rc;

  int nGaussians = orbitalCoefficients.dim1();

  double term_value = 0.0;
  double orig_value = 0.0;
  double orig_firstDerivative = 0.0;
  double orig_secondDerivative = 0.0;

  double p0,p1;

  for (int i=0; i<nGaussians; i++)
    {
      p0 = orbitalCoefficients(i,0);
      p1 = orbitalCoefficients(i,1);
      term_value = p1*exp(-p0*rc_sq);

      orig_value += term_value;

      orig_firstDerivative += -2*p0*rc*term_value;

      orig_secondDerivative += 2*p0*(2*p0*rc_sq-1)*term_value;
    }

  // We determine the sign of the exponential.
  if (phi0-orig_value > 0)
    sgn0 = 1;
  else
    sgn0 = -1;

  // We check if the replacement orbital has to change sign between 0 and rc.
  noC = true;
  C = 0.0;
  if (orig_value*phi0 < 0)
    {
      noC = false;
      if (orig_value > 0.0 && orig_value < 0.1)
	C = orig_value + 0.1;
      else if (orig_value < 0.0 && orig_value > -0.1)
	C = orig_value - 0.1;
      else
	C = 2*orig_value;
    }
  else
    {
      // If the magnitude of the value of the original orbital at the 
      // correction radius is less than 0.1, the fitting function does not have
      // enough space to work with, and the replacement orbitals end up being
      // very strange.  We use the constant shift to put the orbital where it
      // has greater magnitude.

      if (orig_value > 0.0 && orig_value < 0.1)
	{
	  noC = false;
	  C = orig_value - 0.1;
	}
      else if (orig_value < 0.0 && orig_value > -0.1)
	{
	  noC = false;
	  C = orig_value + 0.1;
	}
    }
 
  if (sgn0 == 1 && phi0 < 0.0)
    {
      noC = false;
      if (orig_value > -0.1)
	C = orig_value - 0.1;
      else
	C = 2*orig_value;
    }
  else if (sgn0 == -1 && phi0 > 0.0)
    {
      noC = false;
      if (orig_value < 0.1)
	C = orig_value + 0.1;
      else
	C = 2*orig_value;
    }

  // Set Zeff for calculating the local energy of the replacement orbital.
  Zeff = Z*fabs(1+(n0/phi0));

  // Match the value at rc.
  double x1 = log(fabs(orig_value-C));

  // Match the first derivative at rc.
  double x2 = orig_firstDerivative/(orig_value-C);

  // Match the second derivative at rc.
  double x3 = orig_secondDerivative/(orig_value-C);

  // Cusp condition at the nucleus.  There are two possibilities to prevent 
  // Zeff from being negative.
  double x4;

  if (n0/phi0 >= -1.0)
    x4 = -Z*(phi0+n0)/(phi0-C);
  else 
    x4 = Z*(phi0+n0)/(phi0-C);

  // If the orbital is negative at the origin and has a negative first
  // derivative, we change the sign of the cusp condition.
  if ( (sgn0 == 1 && phi0<0.0) || (sgn0 == -1 && phi0>0.0) )
    x4 *= -1;

  // Value of the replacement orbital at the nucleus.
  double x5 = log(fabs(phi0-C));

  Array1D<double> temp_cs;
  temp_cs.allocate(5);

  double rc_inv = 1/rc;

  temp_cs(0) = x5;
  temp_cs(1) = x4;
  temp_cs(2) = (6*(x1-x5)*rc_inv-3*(x2+x4))*rc_inv+0.5*(x3-x2*x2);
  temp_cs(3) = ((8*(x5-x1)*rc_inv+5*x2+3*x4)*rc_inv+x2*x2-x3)*rc_inv;
  temp_cs(4) =((3*(x1-x5)*rc_inv-2*x2-x4)*rc_inv+0.5*(x3-x2*x2))*rc_inv*rc_inv;

  alpha.initialize(temp_cs);

  temp_cs.deallocate();

  calculateSigmaSq();
}

void QMCElectronNucleusCuspParameters::calculateSigmaSq()
{
  // This function calculates the maximum square deviation of the local energy
  // of the replacement orbital from the ideal curve.

  if (noC == true)
    // If C is zero, then the replacement orbital has no node between 0 and
    // rc.  The local energy will be finite everywhere on the interval, and
    // the search for the maximum deviation will be simple.
      
    findMaxDeviation(0.0,rc,true,sigma_sq);

  else
    {
      // If C is not zero, then the replacement orbital has a node between 0 
      // and rc.  Because the local energy will blow up at this point, we 
      // exclude a small interval centered on the node from the search for the
      // maximum deviation.

      // The polynomial alpha(r)-ln(C) will have a root where the replacement 
      // orbital has a root.

      Array1D<double> temp_coeffs;
      temp_coeffs = alpha.getCoefficients();
      temp_coeffs(0) -= log(fabs(C));

      Polynomial temp_poly;
      temp_poly.initialize(temp_coeffs);

      Array1D<Complex> roots = temp_poly.getRoots();

      double root = -2;
      const double tol = 1e-8;
      for (int i=0; i<roots.dim1(); i++)
	{
	  double re = roots(i).real();
	  double im = roots(i).imaginary();

	  if (re > 0 && re < rc && fabs(im) < tol)
	    {
	      root = re;
	      break;
	    }
	}

      if (root < 0)
	findMaxDeviation(0.0,rc,true,sigma_sq);

      else
	{
	  double max1 = 0.0;
	  double max2 = 0.0;

	  if (root-rc/5 > 0.0)
	    findMaxDeviation(0.0,root-rc/5,true,max1);
      
	  if (root+rc/5 < rc)
	    findMaxDeviation(root+rc/5,rc,false,max2);

	  sigma_sq = max1 > max2 ? max1 : max2;
	}
    }
}

void QMCElectronNucleusCuspParameters::findMaxDeviation(double lower, 
		                double upper, bool lowerBoundZero, double& max)
{
#define GOLD 0.61803399

  // Finds the max deviation of the replacement local energy from the ideal 
  // curve on the interval lower to upper.

  double lowR = lower;
  double repEnergy = calculateLocalEnergy(lowR,lowerBoundZero);
  idealCurve.evaluate(lowR);
  double idealEnergy = Z*Z*idealCurve.getFunctionValue();
  double lowValue = (repEnergy-idealEnergy)*(repEnergy-idealEnergy);

  double highR = upper;
  repEnergy = calculateLocalEnergy(highR,false);
  idealCurve.evaluate(highR);
  idealEnergy = Z*Z*idealCurve.getFunctionValue();  
  double highValue = (repEnergy-idealEnergy)*(repEnergy-idealEnergy);

  double midR = GOLD*(upper-lower)+lower;
  repEnergy = calculateLocalEnergy(midR,false);
  idealCurve.evaluate(midR);
  idealEnergy = Z*Z*idealCurve.getFunctionValue();
  double midValue = (repEnergy-idealEnergy)*(repEnergy-idealEnergy);

  while (highR-lowR > 1e-4)
    {
      if (highValue > lowValue)
	{
	  lowR = midR;
	  lowValue = midValue;
	  if (highValue > midValue)
	    midR = lowR+GOLD*(highR-lowR);
	  else if (highValue <= midValue)
	    midR = highR-GOLD*(highR-lowR);
	}
      else if (highValue < lowValue)
	{
	  highR = midR;
	  highValue = midValue;
	  if (lowValue > midValue)
	    midR = highR-GOLD*(highR-lowR);
	  else if (lowValue <= midValue)
	    midR = lowR+GOLD*(highR-lowR);
	}
      repEnergy = calculateLocalEnergy(midR,false);
      idealCurve.evaluate(midR);
      idealEnergy = Z*Z*idealCurve.getFunctionValue();      
      midValue = (repEnergy-idealEnergy)*(repEnergy-idealEnergy);
    }

  if (highValue > midValue)
    max = highValue > lowValue ? highValue : lowValue;
  else if (midValue >= highValue)
    max = midValue > lowValue ? midValue : lowValue;

#undef GOLD
}

void QMCElectronNucleusCuspParameters::printParameters()
{
  cout << "***ElectronNucleusCuspParameters***" << endl;
  cout << "rc = " << rc << endl;

  Array1D<double> temp_coeffs;
  temp_coeffs = alpha.getCoefficients();
  cout << "alpha coefficients:" << endl;
  for (int i=0; i<temp_coeffs.dim1(); i++)
    cout << "\t" << temp_coeffs(i);
  cout << endl;
 
  temp_coeffs = idealCurve.getCoefficients();
  cout << "ideal curve coefficients:" << endl;
  for (int i=0; i<temp_coeffs.dim1(); i++)
    cout << "\t" << temp_coeffs(i);
  cout << endl;

  cout << "n0 = " << n0 << endl;
  cout << "C = " << C << endl;
  cout << "sgn0 = " << sgn0 << endl;
  cout << "sigmaSq = " << sigma_sq << endl;
  cout << "phi0 = " << phi0 << endl;
  cout << "Z = " << Z << endl;
  cout << "Zeff = " << Zeff << endl;
  cout << "Orbital coefficients:" << endl;
  for (int i=0; i<orbitalCoefficients.dim1(); i++)
    {
      cout << "\t" << orbitalCoefficients(i,0) << "\t";
      cout << orbitalCoefficients(i,1) << endl;
    }
  cout << "The local energy of the replacement orbital at rc = ";
  cout << calculateLocalEnergy(rc,false) << endl;
  idealCurve.evaluate(rc);
  cout << "The ideal local energy at rc = ";
  cout << Z*Z*idealCurve.getFunctionValue() << endl;

  double ovalue,ogradx,ogrady,ogradz,olap,ox,oy,oz,nor;

  ovalue = 0.0;
  ogradx = 0.0;
  ogrady = 0.0;
  ogradz = 0.0;
  olap = 0.0;

  ox = rc;
  oy = 0.0;
  oz = 0.0;
  nor = rc;

  evaluateOriginalOrbital(ox,oy,oz,nor,ovalue,ogradx,ogrady,ogradz,olap);
  cout << "Original orbital value at rc = " << ovalue << endl;
  cout << "Original gradient at rc = " << ogradx << endl;
  cout << "Original laplacian at rc = " << olap << endl;

  ovalue = 0.0;
  ogradx = 0.0;
  ogrady = 0.0;
  ogradz = 0.0;
  olap = 0.0;

  evaluateReplacementOrbital(ox,oy,oz,nor,ovalue,ogradx,ogrady,ogradz,olap);
  cout << "Replacement orbital value at rc = " << ovalue << endl;
  cout << "Replacement gradient at rc = " << ogradx << endl;
  cout << "Replacement laplacian at rc = " << olap << endl;
  cout << endl << endl;
}

ostream& operator <<(ostream& strm, QMCElectronNucleusCuspParameters &rhs)
{
  // This function is never used.  It was added to get rid of some compiler 
  // errors.
  return strm;
}

