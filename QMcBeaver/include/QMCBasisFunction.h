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

#ifndef QMCBasisFunction_H
#define QMCBasisFunction_H

#include <math.h>
#include <iostream>

#include "CubicSplineWithGeometricProgressionGrid.h"
#include "QMCBasisFunctionCoefficients.h"
#include "QMCMolecule.h"
#include "QMCFlags.h"
#include "fastfunctions.h"
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
    Calculates the value of a basis function.

    @param whichBF which basis function to evaluate
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    @param elNumber which electron in X to calculate the basis function for
    @return basis function value
    */

  double getPsi(int whichBF, Array2D <double>& X, int elNumber);


  /**
    Calculates the gradient of a basis function.

    @param whichBF which basis function to evaluate
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    @param elNumber which electron in X to calculate the basis function for
    @return basis function gradient value
    */

  Array1D <double> getGradPsi(int whichBF, Array2D <double>& X, int elNumber);


  /**
    Calculates the laplacian of a basis function.

    @param whichBF which basis function to evaluate
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    @param elNumber which electron in X to calculate the basis function for
    @return basis function laplacian value
    */
  double getLaplacianPsi(int whichBF, Array2D <double>& X, int elNumber);


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
    */

  friend istream& operator >>(istream& strm,  
      QMCBasisFunction &rhs);    //read in coefficients


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
  
private:

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

  /**
     Radial part of a basis function.  A basis function 
     can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where \f$\theta(r)\f$
     is the radial portion of the basis function.  The function is in
     terms of \f$r^2\f$ instead of \f$r\f$ to make it's evaluation faster.  
     To evaluate \f$\theta(r)\f$, evaluate the function at \f$r^2\f$.
  */
  static double radialFunction(QMCBasisFunctionCoefficients& BFC, int orbital, 
			double r_sq);

  /**
     First derivative of the radial part of a basis function.  A basis 
     function can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where 
     \f$\theta(r)\f$ is the radial portion of the basis function.  
     The function is in terms of \f$r^2\f$ instead of \f$r\f$ to make it's 
     evaluation faster.  To evaluate \f$\theta'(r)\f$, evaluate the 
     function at \f$r^2\f$.
  */
  static double radialFunctionFirstDerivative(
				       QMCBasisFunctionCoefficients& BFC, 
				       int orbital, double r_sq);

  /**
     Second derivative of the radial part of a basis function.  A basis 
     function can be factored into \f$x^{i}y^{j}z^{k}\theta(r)\f$ where 
     \f$\theta(r)\f$ is the radial portion of the basis function.  
     The function is in terms of \f$r^2\f$ instead of \f$r\f$ to make it's 
     evaluation faster.  To evaluate \f$\theta'(r)\f$, evaluate the 
     function at \f$r^2\f$.
  */
  static double radialFunctionSecondDerivative(
					QMCBasisFunctionCoefficients& BFC, 
					int orbital, double r_sq);

  double basis_function(int which_BFC, 
			QMCBasisFunctionCoefficients& BFC, 
			int orbital, Array1D <double>& X);

  Array1D <double> grad_basis_function(int which_BFC, 
					    QMCBasisFunctionCoefficients& BFC, 
					    int orbital, Array1D <double>& X);

  double laplacian_basis_function(int which_BFC,
				  QMCBasisFunctionCoefficients& BFC, 
				  int orbital, Array1D <double>& X);


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
  
