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
  int Ncharge;
  int Nelectrons;
  int Ndeterminants;
  double factor;
  string trialFunctionType;

  /**
     For the matrix multiplications, we'll want to have
     the coefficients separated into matrices, one for
     each determinant.
  */
  Array1D< Array2D<qmcfloat> > AlphaCoeffs;
  Array1D< Array2D<qmcfloat> > BetaCoeffs;

  /**
     This function will convert the orbital occupation
     arrays into one of two formats.

     If unordered, the original setting, then a determinant
     will indicate that it uses a particular orbital with a 1,
     and 0 otherwise.

     If ordered, then a determinant will indicate that it uses
     a particular orbital with the index that it will show up in
     for that particular determinant, and -1 otherwise.
     
     The ordered format is required for optimizing the orbitals.

     When the wavefunction is printed, it will be autoamtically
     converted to the unordered format.
  */
  void sortOccupations(bool ordered);

  /**
     This function will convert a restricted wavefunction
     into an unrestricted one by duplicating any orbitals
     that are shared between alpha and beta electrons.

     It can not be undone.

     Print out the occupation arrays to clarify what this
     is doing.
  */
  void unlinkOrbitals();

  /**
     This function will duplicate orbitals so that
     no same spin determinants share orbitals.
     
     It can not be undone.

     Print out the occupation arrays to clarify what this
     is doing.
  */
  void unlinkDeterminants();

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
     The number of active orbitals
     is the number actually used. This
     number might be useful to prevent
     calculating unneeded data.
     
     @return number of active orbitals.
  */
  int getNumberActiveOrbitals();

  /**
     The number of orbitals used by
     alpha electrons.
     @return number of active orbitals.
  */
  int getNumberActiveAlphaOrbitals();

  /**
     The number of orbitals used by beta
     electrons.
     @return number of active orbitals.
  */
  int getNumberActiveBetaOrbitals();

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
     This function indicates whether the occupation arrays
     use -1 or 0 to indicate an unused orbital.
  */
  int getUnusedIndicator();

  /**
     Break apart the OrbitalCoeffs array into
     individual AlphaCoeffs and BetaCoeffs matrices,
     one for each determinant.
  */
  void makeCoefficients();

  /**
     Return a matrix with all the coefficients
     Since we broke apart the matrices, this is
     how we request them.
  */
  Array2D<double> * getCoeff(int ci, bool isAlpha);

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
  Array2D<qmcfloat> OrbitalCoeffs;

  /**
    Array containing the CI coefficients for a multideterminant wavefunction.
  */
  Array1D<double> CI_coeffs;

  /**
    Array which indicates how many \f$\alpha\f$ spin electron are in each
    orbital for each determinant.  

    The dimmensions are Ndeterminants x Norbitals
  */
  Array2D<int> AlphaOccupation;

  /**
    Array which indicates how many \f$\beta\f$ spin electron are in each
    orbital for each determinant.

    The dimmensions are Ndeterminants x Norbitals
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
    @param charge the overall electronic charge of the molecule
    @param numberOrbitals number of orbitals.
    @param numberBasisFunctions number of basis functions.
    @param numberDeterminants number of determinants.
    @param runfile file to load the object state from.
  */
  void read(int charge, int numberOrbitals, int numberBasisFunctions, 
	    int numberDeterminants, string functionType, string runfile);

  /**
     The number of optimizable CI parameters.
  */
  int getNumberCIParameters();

  /**
     The number of optimizable orbital coefficients.
  */
  int getNumberORParameters();

  /**
     Fill the params array with the current CI parameters,
     starting at the shift position.
  */
  void getCIParameters(Array1D<double> & params, int shift);

  /**
     Fill the params array with the current orbitals parameters,
     starting at the shift position.
  */
  void getORParameters(Array1D<double> & params, int shift);

  /**
     Replace the current CI parameters with our new values.
     The shift variable tells us where in the array to look.
  */
  void setCIParameters(Array1D<double> & params, int shift);

  /**
     Replace the current orbital parameters with our new values.
     The shift variable tells us where in the array to look.
  */
  void setORParameters(Array1D<double> & params, int shift);

};

#endif




