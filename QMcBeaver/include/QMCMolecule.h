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

#ifndef QMCMolecule_H
#define QMCMolecule_H

#include <string>
#include <fstream>
#include <iostream>
#include <list>

#include "Array1D.h"
#include "Array2D.h"
#include "StringManipulation.h"

using namespace std;


/**
   Describes a particular molecular geometry.  The geometry is defined by
   3-dimensional cartesian coordinates for each atom, with specified charges 
   and types. 
*/

class QMCMolecule
{
 public:
  /**
     Creates an instance of the class.
  */

  QMCMolecule();


  /**
     Initializes the object.

     @param nAtoms number of atoms in the molecule.
  */

  void initialize(int nAtoms);


  /**
     Gets the number of atoms in the molecule.

     @return number of atoms in the molecule.
  */

  int getNumberAtoms();


  /**
     Array containing the labels for the atoms.  The ith element is the
     label for the ith atom.
  */
  Array1D <string> Atom_Labels;


  /**
     Array containing the 3-dimensional cartesian positions for the atoms.
     The ith element is the position for the ith atom.
  */

  Array2D <double> Atom_Positions;


  /**
     Array containing the nuclear charges for the atoms.  The ith element is 
     the charge for the ith atom.
  */

  Array1D <int>  Z;


  /**
     Array containing all of the different atom labels used in the molecule.
  */

  Array1D <string> NucleiTypes;


  /**
     Sets two QMCMolecule objects equal.

     @param rhs object to set this object equal to.
  */

  QMCMolecule operator=( const QMCMolecule & rhs );

  /**
     Finds the index of the nucleus closest to a specified particle.

     @param x a 2D array containing particle coordinates.
     @param index index of the particle for which to find the closest nucleus.

     @return index of the nucleus closest to the specified particle.
  */
  int findClosestNucleusIndex(Array2D<double> & x, int index);

  /**
     Loads the state of the object from an input stream.
  */

  friend istream& operator >>(istream& strm, QMCMolecule &rhs);


  /**
     Writes the state of the object to an output stream.
  */

  friend ostream& operator <<(ostream& strm, QMCMolecule &rhs);


  /**
     Loads the state of the object from a file.

     @param runfile file to load the object state from.
  */

  void read(string runfile);

 private:
  int Natoms;
};

#endif
