
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

#ifndef QMCElectronNucleusCusp_H
#define QMCElectronNucleusCusp_H

#include <iostream>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCInput.h"
#include "QMCElectronNucleusCuspParameters.h"
#include "QMCBasisFunction.h"
#include "QMCBasisFunctionCoefficients.h"
#include "QMCMolecule.h"
#include "Polynomial.h"
#include <math.h>

using namespace std;

/**
   A class that replaces Gaussian orbitals with exponential functions near the 
   nuclei.  The new function is defined as \( \phi=C+sgn0*exp(p(r)) \), where 
   \(C\) is a constant, \(sgn0=+/-1\), and \(p(r)\) is a polynomial in \(r\).  
   The value and first and second derivatives of the exponential orbital are 
   fit to the original orbital at the radius of correction.  
*/

class QMCElectronNucleusCusp
{
 private:

  QMCInput* Input;

  /**
    This array holds the parameters of the exponential functions for each
    nucleus and orbital.  The indexes are (nucleus,orbital).
  */
  Array2D<QMCElectronNucleusCuspParameters> ORParams;

  /**
    This is the array of wavefunction coefficients.
  */
  Array2D<qmcfloat> ORWF_coeffs;

  QMCBasisFunction* BF;
  QMCBasisFunctionCoefficients* BF_coeffs;
  QMCMolecule* Molecule;

  /**
    This array holds the exponential and expansion coefficients for the part of
    an orbital arising from the s type Gaussian basis functions centered on a
    nucleus.
  */
  Array2D<qmcfloat> sTypeCoeffs;

  /**
    The ideal curve is a benchmark for judging the behavior of the local 
    energy of the original and replacement orbitals.
  */
  Polynomial idealCurve;

  /**
    The array of ideal curve parameters.
  */
  Array1D<double> idealCurveParams;

  /**
    The charge of the nucleus we are considering.
  */
  int Z;

  int norbitals,natoms;

  /**
    The radius of correction for a certain nucleus and orbital.
  */
  double rc;

  /**
    The value of the orbital due to basis functions centered on other nuclei
    at this nucleus.
  */
  double n0;

  /**
    Determines the radius of correction for a certain nucleus and orbital.  
    The radius of correction is where the local energy of the orbital deviates
    from the ideal curve by more than a tolerance amount.  Fits the ideal curve
    at rc.
  */
  void determineRc();
  
  /**
    Gets the exponential and expansion coefficients of the part of the orbital
    due to the s type basis functions centered on the nucleus.
    Also calculates n0, which is the value of the orbital due to basis
    functions centered on other nuclei at this nucleus.

    @param nuc the index of the nucleus
    @param orb the index of the orbital
  */
  void setSTypeCoeffs(int nuc, int orb);

  /**
    Calculates the local energy of the original orbital at r_orig.

    @param r_orig the distance from the nucleus.
    @return the local energy of the original orbital at r_orig
  */
  double getOrigLocalEnergy(double r_orig);

  /**
    Evaluates the original orbital at r_orig and returns its value.

    @param r_orig the distance from the nucleus
    @return the value of the original orbital at r_orig
  */
  double evaluateOrigOrbital(double r_orig);

  /**
    Fits the replacement orbital to the original orbital at rc.  Minimizes the 
    deviation of the local energy of the replacement orbital from the ideal 
    curve with respect to its value at the nucleus.

    @return ElectronNucleusCuspParameters the replacement orbital for a nucleus
    and orbital
  */
  QMCElectronNucleusCuspParameters fitOrbitalParameters();

 public:

  QMCElectronNucleusCusp();

  void initialize(QMCInput* input, const Array2D<qmcfloat>& WFCoeffs);
  
  /**
    Replaces the entries in the Slater, gradient, and laplacian matrices 
    arising from electrons that are too close to nuclei.

    @param X the array of positions of the electrons
    @param D the Slater matrix
    @param Grad_D the gradient matrix
    @param Laplacian_D the laplacian matrix
  */
  void replaceCusps(Array2D<double> & X,
		    int Start, int Stop,
		    Array2D<qmcfloat> & D, 
		    Array2D<qmcfloat> & GradX,
		    Array2D<qmcfloat> & GradY,
		    Array2D<qmcfloat> & GradZ,
		    Array2D<qmcfloat> & Laplacian_D);

  /**
    Fits exponential functions for each orbital in the region of each nucleus.
  */
  void fitReplacementOrbitals();

  /**
    Sets this object equal to another QMCElectronNucleusCusp
    
    @param rhs the object to set this one equal to
  */
  void operator=(const QMCElectronNucleusCusp& rhs);

  /**
    Formats and prints the replacement orbitals to a stream
  */
  friend ostream& operator <<(ostream& strm, QMCElectronNucleusCusp &rhs);
};
#endif
