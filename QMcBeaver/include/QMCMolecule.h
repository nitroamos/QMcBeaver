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
  void initialize(int nAtoms, int GridLevel);

  /**
     Gets the number of atoms in the molecule.

     @return number of atoms in the molecule.
  */
  int getNumberAtoms();

  /**
     Gets the total charge from all the nuclei.

     @return total charge
  */
  int getNuclearCharge();

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
     If we used pseudopotentials, then we removed Z-Zeff electrons.
  */
  Array1D <int>  Zeff;

  /**
     Each element in the array corresponds to a nucleus index. If the
     element is true, then that nucleus uses an ecp.
  */
  Array1D<bool> usesPseudo;

  /**
     Save some of the text description of the ecp so that we
     can write it in the restart file.
  */
  Array2D< string > pseudoTitle;

  /**
     The local part of the ecp.
  */
  Array1D< Array2D<double> > Vlocal;

  /**
     The nonlocal part of the ecp.
  */
  Array1D< Array1D< Array2D<double> > > Vnonlocal;

  /**
     This function helps to evaluate Vlocal and Vnonlocal
     at distance r.

     @param V a slice of Vlocal or Vnonlocal
     @param r the value at which to evaluate
     @return the value of the potential
  */
  double evaluatePotential(Array2D<double> & V, double r);

  /**
     Each grid point from the Lebedev-Laikov functions
     has an associated weight, stored here.
  */
  Array1D< Array1D<double> > gridWeights;

  /**
     We store the grid as points on the unit sphere centered
     at the origin. This function is called to provide scales
     and translated grid points.

     @param nuc the index of the nucleus
     @param R the distance of the points from the nucleus
     @param translate should be false if you want to keep
     it centered at the origin
  */
  Array2D<double> getGrid(int nuc, int rg, double R, bool translate);

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
  void readGeometry(string runfile);

  /**
     Read the pseudopotentials, if any, in an input file. The format
     should be the same as GAMESS.

     @return 1 if there were some pseudopotentials read in
  */
  int readPseudoPotential(string runfile);
  int getNumL(int nuc){ return numL(nuc); };
  const static int numRandomGrids = 50;
 private:
  int Natoms;
  
  Array1D<int> numL;

  /**
     This will be initialized to the same value as:
     globalInput.flags.pseudo_gridLevel
     And indicates how many grid points to use.
  */
  int gridLevel;

  /**
     We access the Lebedev-Laikov functions once, and
     save the grid points in these arrays. indexed:
     num nuclei x num grids : num grid points x 3
  */
  Array2D< Array2D<double> > grid;
};

#endif
