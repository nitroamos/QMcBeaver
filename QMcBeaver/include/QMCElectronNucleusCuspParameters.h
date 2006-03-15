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

#ifndef QMCElectronNucleusCuspParameters_H
#define QMCElectronNucleusCuspParameters_H

#include <iostream>

#include "Polynomial.h"
#include "Complex.h"
#include "Array2D.h"
#include <math.h>

using namespace std;

/**
   A container for the parameters for replacing a Gaussian orbital with an 
   exponential function in the region of a nucleus.  The new function is 
   defined as \( \phi=C+sgn0*exp(p(r)) \), where \(C\) is a constant, 
   \(sgn0=+/-1\), and \(p(r)\) is a polynomial in \(r\). 
*/

class QMCElectronNucleusCuspParameters
{
 private:
  /**
    The radius of correction.
  */
  double rc;

  /**
    The polynomial in r.
  */
  Polynomial alpha;

  /**
    The ideal local energy of the orbital.
  */
  Polynomial idealCurve;

  /**
    This is a constant shift so that the exponential always has the same sign.
  */
  double C;

  /**
    This variable is true if C == 0.0 (Parts of some formulas cancel out if
    this is true, and equality conditions involving doubles are bad form)
  */
  bool noC;

  /**
    sgn0=+/-1 depending on the sign of the orbital at the nucleus.
  */
  int sgn0;

  /**
    The square of the max deviation of the local energy of the corrected 
    orbital from the ideal curve.
  */
  double sigma_sq;

  /**
    The value of the corrected orbital at the nucleus.
  */
  double phi0;

  /**
    The value of the part of the orbital due to basis functions centered on 
    other nuclei evaluated at this nucleus.
  */
  double n0;

  /**
    The charge of the nucleus.
  */
  int Z;

  /**
    The effective nuclear charge used to calculate the local energy of the 
    replacement orbital.
  */
  double Zeff;

  /**
    The exponents and expansion coefficients of the part of the orbital due to 
    s type Gaussians centered on this nucleus.  This array is used to 
    calculate the contribution of this part of the original orbital when an 
    electron is within the radius of correction.  It is not used to fit the 
    parameters of the replacement orbital.
  */
  Array2D<qmcfloat> orbitalCoefficients;

  /**
    Calculates the value, gradient, and laplacian of the original orbital.  
    These values must be subtracted from the appropriate elements of the Slater
    determinant and replaced with those of the new orbital.

    @param x the x coordinate of the electron
    @param y the y coordinate of the electron
    @param z the z coordinate of the electron
    @param r the distance of the electron from the nucleus
    @param orig_value the value of the original orbital
    @param orig_gradx the x gradient of the original orbital
    @param orig_grady the y gradient of the original orbital
    @param orig_gradz the z gradient of the original orbital
    @param orig_laplacian the laplacian of the original orbital
  */
  void evaluateOriginalOrbital(double x, double y, double z, double r,
double& orig_value, double& orig_gradx, double& orig_grady, double& orig_gradz,
			                               double& orig_laplacian);

  /**
    Calculates the value, gradient, and laplacian of the replacement orbital.
    These values replace those calculated from the s type Gaussian basis 
    functions centered on this nucleus.

    @param x the x coordinate of the electron
    @param y the y coordinate of the electron
    @param z the z coordinate of the electron
    @param r the distance of the electron from the nucleus
    @param rep_value the value of the replacement orbital
    @param rep_gradx the x gradient of the replacement orbital
    @param rep_grady the y gradient of the replacement orbital
    @param rep_gradz the z gradient of the replacement orbital
    @param rep_laplacian the laplacian of the replacement orbital
  */
  void evaluateReplacementOrbital(double x, double y, double z, double r,
    double& rep_value, double& rep_gradx, double& rep_grady, double& rep_gradz,
				                        double& rep_laplacian);

  /**
    Calculates the local energy of the replacement orbital at r. 
  */
  double calculateLocalEnergy(double r, bool rIsZero);

  /**
    Calculates the square deviation of the local energy of the replacement 
    orbital from the ideal curve.
  */
  void calculateSigmaSq();

  /**
    Finds the maximum deviation of the local energy of the replacement orbital
    from the ideal curve on the interval lower to upper.

    @param lower the lower bound of the interval
    @param upper the upper bound of the interval
    @param lowerBoundZero is true if the lower bound is zero
    @param max the maximum deviation on the interval
  */
  void findMaxDeviation(double lower, double upper, bool lowerBoundZero, 
			double& max);

 public:

  /**
    Initializes the parameters with default values.
  */
  QMCElectronNucleusCuspParameters();

  /**
    Replaces the original orbital values with those of the replacement orbital
    in the Slater determinant.
    
    @param x the x coordinate of the electron
    @param y the y coordinate of the electron
    @param z the z coordinate of the electron
    @param r the distance from the electron to the nucleus
    @param orb_value the value of the orbital
    @param orb_gradx the x gradient of the orbital
    @param orb_grady the y gradient of the orbital
    @param orb_gradz the z gradient of the orbital
    @param orb_laplacian the laplacian of the orbital
  */
  void replaceOrbitalValues(double x, double y, double z, double r,
    double& orb_value, double& orb_gradx, double& orb_grady, double& orb_gradz,
                                                        double& orb_laplacian);

  /**
    This function fits the parameters of the replacement orbital based on the 
    value and derivatives of the original orbital at the radius of correction.

    @param temp_phi0 the value of the replacement orbital at the nucleus.
  */
  void fitReplacementOrbital(double temp_phi0);

  /**
    Sets the coefficients of the original orbital, the nuclear charge, and the
    value of the part of the orbital due to basis functions centered on other
    nuclei at this nucleus.

    @param temp_OrbitalCoefficients the part of the orbital arising from s type
    Gaussian basis functions centered on this nucleus.
    @param temp_idealCurve the ideal local energy curve
    @param temp_Z the charge of the nucleus
    @param temp_n0 the value of the part of the orbital due to basis functions
    centered on other nuclei at this nucleus.
  */
  void initialize(const Array2D<qmcfloat>& temp_OrbitalCoefficients, 
double temp_n0, double temp_rc, const Polynomial& temp_idealCurve, int temp_Z);

  /**
    Gets the radius of correction for this replacement orbital.

    @return rc the radius of correction
  */
  double get_rc();

  /**
    Sets the radius of correction for this replacement orbital.
    
    @param temp_rc the radius of correction
  */
  void set_rc(double temp_rc);

  /**
    Gets the square of the maximum deviation of the replacement local energy 
    from the ideal curve.

    @return sigma_sq the square of the maximum deviation of the replacement 
    local energy from the ideal curve
  */
  double getSigmaSq();

  /**
    Sets this object equal to another one.

    @param rhs the set of parameters to set this object equal to
  */
  void operator=(const QMCElectronNucleusCuspParameters& rhs);

  /**
    Prints out the parameters of this replacement orbital.
  */
  void printParameters();

  /**
    Formats and prints the replacement orbital to a stream
  */
  friend ostream& operator <<(ostream& strm, QMCElectronNucleusCuspParameters &rhs);
};

#endif
