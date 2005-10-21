//            Qmcbeaver
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

/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   

Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/


#ifndef QMCBasisFunction_H
#define QMCBasisFunction_H

#include <math.h>
#include <iostream>

#include "CubicSplineWithGeometricProgressionGrid.h"
#include "QMCBasisFunctionCoefficients.h"
#include "QMCMolecule.h"
#include "QMCFlags.h"
#include "Array1D.h"
#include "Array2D.h"

using namespace std;


/**
   This class  stores all of the parameters that a gaussian basis set is
   constructed from for a MOLECULE.  This contains a QMCBasisFunctionCoefficent
   for each atom type.  
*/

class QMCBasisFunction
{
public:
  friend class QMCPsiPotential;
  /**
     Creates an instance of the class.
  */
  QMCBasisFunction();
  
  /**
     Initializes the class with data input to control the calculation and
     provide the molecular geometry.
     
     @param flags input control information
     @param molecule information about the specific molecule
  */
  void initialize(QMCFlags * flags, QMCMolecule * molecule);
  
  /**
     Calculates the value, gradient, and Laplacian of a basis funcion.
  */
  void evaluateBasisFunctions(Array2D<double>& X, int start, int stop,
			      Array2D<qmcfloat>& chi_value,
			      Array2D<qmcfloat>& chi_grx,
			      Array2D<qmcfloat>& chi_gry,
			      Array2D<qmcfloat>& chi_grz,
			      Array2D<qmcfloat>& chi_laplacian);
  
  
  /**
     Sets two QMCBasisFunctions objects equal.
     
     @param rhs object to set this object equal to
  */
  void operator=( const QMCBasisFunction & rhs );
  
  /** 
      Loads the state of the object from a file.
      
      @param runfile file to load
  */
  void read(string runfile);
  
  /**
     Loads the state of the object from an input stream.
     read in coefficients
  */
  friend istream& operator >>(istream& strm,  
			      QMCBasisFunction &rhs);
  
  /**
     Writes the state of the object to an output stream.
  */
  friend ostream& operator <<(ostream& strm,QMCBasisFunction& rhs);
  
  /**
     Returns how many basis functions are located on a specific atom.  
     This can probably be depricated once we have a good initialization
     scheme and not MikesJacked one.
     
     @param i index of atom
     @return number of basis functions on the atom
  */
  int getNumberBasisFunctions(int i);
  
protected:  
  QMCFlags *flags;
  QMCMolecule *Molecule;
  
  int N_BasisFunctions;
  
  /**
     Array containing the QMCBasisFunctionCoefficients for all the atoms.
     Each element is for a different atom.
  */
  Array1D<QMCBasisFunctionCoefficients> BFCoeffs;  // Container for Coeffs

  Array1D <double> Xcalc;      // e position relative to nucleus 
  
  Array2D<int> BFLookupTable;   // Lookup Table to select a BF
                                // BFLookupTable[BF_number][0=atom#,1=orb#]
private:  
  /**
     Radial part of a basis function.  A basis function
     can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where \f$\theta(r)\f$
     is the radial portion of the basis function.  The function is in
     terms of \f$r^2\f$ instead of \f$r\f$ to make it's evaluation faster.
     To evaluate \f$\theta(r)\f$, evaluate the function at \f$r^2\f$.
  */
  static double radialFunction(QMCBasisFunctionCoefficients& BFC,
			       int orbital, double r_sq);
  
  /**
     First derivative of the radial part of a basis function.  A basis
     function can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where
     \f$\theta(r)\f$ is the radial portion of the basis function.
     The function is in terms of \f$r^2\f$ instead of \f$r\f$ to make it's
     evaluation faster.  To evaluate \f$\theta'(r)\f$, evaluate the
     function at \f$r^2\f$.
  */
  static double radialFunctionFirstDerivative(QMCBasisFunctionCoefficients& BFC,
					      int orbital, double r_sq);
  
  /**
     Second derivative of the radial part of a basis function.  A basis
     function can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where
     \f$\theta(r)\f$ is the radial portion of the basis function.
     The function is in terms of \f$r^2\f$ instead of \f$r\f$ to make it's
     evaluation faster.  To evaluate \f$\theta'(r)\f$, evaluate the
     function at \f$r^2\f$.
  */
  static double radialFunctionSecondDerivative(QMCBasisFunctionCoefficients& BFC,
					       int orbital, double r_sq);
  
  /**
     Interpolation of the radial part of a basis function.  A basis function 
     can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where \f$\theta(r)\f$
     is the radial portion of the basis function.  The interpolation is in
     terms of \f$r^2\f$ instead of \f$r\f$ to make it's evaluation faster.  
     To evaluate \f$\theta(r)\f$, evaluate the interpolation at \f$r^2\f$.
  */
  Array2D<CubicSplineWithGeometricProgressionGrid> RadialFunctionInterpolation;
  
  /**
     Interpolation of the first derivative of the radial part of a basis 
     function.  A basis function can be factored into 
     \f$x^{i}y^{j}z^{k}\theta(r)\f$ where \f$\theta(r)\f$
     is the radial portion of the basis function.  The interpolation is in
     terms of \f$r^2\f$ instead of \f$r\f$ to make it's evaluation faster.  
     To evaluate \f$\theta'(r)\f$, evaluate the interpolation at \f$r^2\f$.
  */
  Array2D<CubicSplineWithGeometricProgressionGrid> 
    RadialFunctionFirstDerivativeInterpolation;
  
  /**
     Interpolation of the second derivative of the radial part of a basis 
     function.  A basis function can be factored into 
     \f$x^{i}y^{j}z^{k}\theta(r)\f$ where \f$\theta(r)\f$
     is the radial portion of the basis function.  The interpolation is in
     terms of \f$r^2\f$ instead of \f$r\f$ to make it's evaluation faster.  
     To evaluate \f$\theta''(r)\f$, evaluate the interpolation at \f$r^2\f$.
  */
  Array2D<CubicSplineWithGeometricProgressionGrid> 
    RadialFunctionSecondDerivativeInterpolation;
  
  /**
     Boolean flag indicating if the interpolated radial component of the
     basis functions should be used.
  */
  bool use_radial_interpolation;
  
  /**
     Initializes the interpolations for the radial components of all basis
     functions.
  */
  void initializeInterpolations();
  
  /**
     Initializes the interpolation for the radial component of a basis 
     function or it's derivatives.  
     
     @param S class used to interpolate the function.
     @param whichDerivative the derivative of the radial function to 
     calculate the interpolation for.
  */
  void initializeInterpolation(int bfc_number,int orbital, 
			       CubicSplineWithGeometricProgressionGrid & S,
			       int whichDerivative);
  
};

#endif

