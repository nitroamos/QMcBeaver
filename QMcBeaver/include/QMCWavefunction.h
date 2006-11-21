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

#ifndef QMCWavefunction_H
#define QMCWavefunction_H

#include <string>
#include <fstream>
#include <iostream>

#include "Array1D.h"
#include "Array2D.h"

using namespace std;

/**
  The coefficients and parameters describing the trial wavefunction for 
  the system.  These are the coefficients for a wavefunction obtained 
  through standard means (RHF, UHF, MCSCF, DFT, etc).
*/

class QMCWavefunction
{
 private:

  int Norbitals;
  int Nbasisfunc;
  int Nalpha;
  int Nbeta;
  int Nelectrons;
  int Ndeterminants;
  double factor;
  string trialFunctionType;

 public:

  /**
    Creates an instance of the class.
  */
  QMCWavefunction();

  /**
    Gets the number of orbitals.
    @return number of orbitals.
  */
  int getNumberOrbitals();

  /**
    Gets the number of basis functions.
    @return number of basis functions.
  */
  int getNumberBasisFunctions();

  /**
    Gets the number of \f$\alpha\f$ spin electrons.
    @return number of \f$\alpha\f$ spin electrons.
  */
  int getNumberAlphaElectrons();

  /**
    Gets the number of \f$\beta\f$ spin electrons.
    @return number of \f$\beta\f$ spin electrons.
  */
  int getNumberBetaElectrons();

  /**
    Gets the total number of electrons.
    @return total number of electrons.
  */
  int getNumberElectrons();
  
  /**
    Gets the number of determinants in the SCF wavefunction.
    @return number of determinants.
  */
  int getNumberDeterminants();

  /**
     This scales all the values in the orbitals by a constant factor.
     Since the wavefunction does not have to be normalized in a QMC
     wavefunction, this will not change anything.

     Numerically however, the calculation may perform better if
     scaled such that the value of the wavefunction is near 1.0.
  */
  void scaleCoeffs(double factor);

  /**
    Array containing the coefficients used to construct the alpha orbitals.

    For example, orbitals are constructed so that
    \f[
    Orbital_{i}(x,y,z) = \sum_{j=0}^{NumberBasisFunctions-1} Coeffs_{i,j} 
    BasisFunction_{j}(x,y,z)
    \f]
    where the the \f$BasisFunction_{j}(x,y,z)\f$ are from QMCBasisFunction.
    It is assumed that the ordering of the coefficients is the same as
    the basisfunctions in the input file.
    If the trial wavefunction is restricted, the alpha and beta orbitals will
    be identical.
  */
  Array2D<qmcfloat> AlphaCoeffs;

  /**
    Array containing the coefficients used to construct the beta orbitals.
    If the trial wavefunction is restricted, the alpha and beta orbitals will
    be identical.
  */
  Array2D<qmcfloat> BetaCoeffs;

  /**
    Array containing the CI coefficients for a multideterminant wavefunction.
  */
  Array1D<double> CI_coeffs;

  /**
    Array which indicates how many \f$\alpha\f$ spin electron are in each
    orbital for each determinant.  
  */
  Array2D<int> AlphaOccupation;

  /**
    Array which indicates how many \f$\beta\f$ spin electron are in each
    orbital for each determinant.
  */
  Array2D<int> BetaOccupation;

  /**
    Sets two QMCWavefunction objects equal.
    @param rhs object to set this object equal to.
  */
  QMCWavefunction operator=( const QMCWavefunction & rhs );

  /**
    Loads the state of the object from an input stream.
  */
  friend istream& operator >>(istream& strm, QMCWavefunction &rhs);

  /**
    Writes the state of the object to an output stream.
  */
  friend ostream& operator <<(ostream& strm, QMCWavefunction& rhs);

  /**
    Loads the state of the object from a file.
    @param numberOrbitals number of orbitals.
    @param numberBasisFunctions number of basis functions.
    @param numberDeterminants number of determinants.
    @param runfile file to load the object state from.
  */
  void read(int numberOrbitals, int numberBasisFunctions, 
	    int numberDeterminants, string functionType, string runfile);
};

#endif




