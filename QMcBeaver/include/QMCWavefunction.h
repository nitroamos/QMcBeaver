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
  through standard means (HF, DFT, etc.).
*/


class QMCWavefunction
{
 private:

  int Norbitals;
  int Nbasisfunc;
  int Nalpha;
  int Nbeta;
  int Nelectrons;

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
     Array containing the coefficients used to construct the orbitals.

     For example, orbitals are constructed so that
     \f[
     Orbital_{i}(x,y,z) = \sum_{j=0}^{NumberBasisFunctions-1} Coeffs_{i,j} 
     BasisFunction_{j}(x,y,z)
     \f]
     where the the \f$BasisFunction_{j}(x,y,z)\f$ are from QMCBasisFunction.
     It is assumed that the ordering of the coefficients is the same as
     the basisfunctions in the input file.
  */

  Array2D <double> Coeffs;


  /**
     Array which indicates how many \f$\alpha\f$ spin electron are in each
     orbital for the wavefunction.
  */

  Array1D<int> AlphaOccupation;


  /**
     Array which indicates how many \f$\beta\f$ spin electron are in each
     orbital for the wavefunction.
  */

  Array1D<int> BetaOccupation;


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
     @param runfile file to load the object state from.
  */

  void read(int numberOrbitals, int numberBasisFunctions, string runfile);
};

#endif




