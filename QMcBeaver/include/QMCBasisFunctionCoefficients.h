
#ifndef QMCBasisFunctionCoefficients_H
#define QMCBasisFunctionCoefficients_H

#include <string>
#include <iostream>
#include <fstream>
#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "StringManipulation.h"

using namespace std;

/**
  This class  stores all of the parameters that a gaussian basis set is
  constructed from for an ATOM.  
  For example, a gaussian basis function is 
  \f[
  Gbf(x,y,z)=x^{k}y^{l}z^{m}\sum_{i=0}^{Ngaussians-1}a_{i}e^{-b_{i}r^{2}}
  \f]
  where k,l,m are determined by the type of basis function, \f$a_{i}\f$
  is the contraction parameter, and \f$b_{i}\f$ is the exponential parameter.
  The particular contraction parameter is chosen so that the basis function
  is normalized.  This is slightly different than what is common with linear
  algebra quantum mechanics programs.  The contraction parameters used here 
  can be obtained using the contraction and exponential parameters and 
  k,l,m from a linear algebra basis file.  You will have to look up the
  formula for doing this.

  This reads in basis function coefficients in the following format...

  \verbatim
  AtomLabel  Number_of_orbitals Maximum_Gaussians

  Ngaussians  Type
  exp_param   contraction_param
  ...         ...

  Ngaussians  Type
  exp_param   contraction_param
  ...         ...

  etc...
  \endverbatim

*/

class QMCBasisFunctionCoefficients
{
 public:
  /**
    Creates an instance of the class.
    */

  QMCBasisFunctionCoefficients();


  /**
    Gets the number of basis functions.

    @return number of basis functions
    */

  int getNumberBasisFunctions();


  /**
    Array containing the parameters for the basis functions where
    Coeffs[bf #][Gaussian #][0=exp,1=contract]
    */

  Array3D <qmcfloat> Coeffs;


  /**
    Array containing the k,l,m parameters which indicate the 
    "angular momentum state" of the basis function 
    (\f$bf=x^{k}y^{l}z^{m}*RadialFunction(r)\f$) where xyz[bf #][0=k,1=l,2=m].
    For example, a "px" orbital would have \f$(k,l,m)=(1,0,0)\f$.
    */

  Array2D <int> xyz_powers;
                           
                           
  /**
    Array containing the number of gaussians that need to be contracted
    for the radial portion of the basis function 
    (\f$bf=x^{k}y^{l}z^{m}*RadialFunction(r)\f$) where N_Gauss[bf #].
    */

  Array1D <int> N_Gauss;      
                              
                              
  /**
    Array containing the type of the basis function where Type[bf #].
    The type is a string representation of the "angular momentum state."
    For example, "px", "dxy", and "fxxx" are all types of basis functions.
    */
 
  Array1D <string> Type;      

                              
  /**
    Sets two QMCBasisFunctionCoefficients objects equal.

    @param rhs object to set this object equal to
    */

  void operator=( const QMCBasisFunctionCoefficients & rhs);


  /**
    Loads the state of the object from an input stream.
    */

  friend istream& operator >>(istream& strm,
      QMCBasisFunctionCoefficients &rhs);    


  /**
    Writes the state of the object to an output stream.
    */

  friend ostream& operator <<(ostream& strm,
      QMCBasisFunctionCoefficients& rhs);


  /**
    Loads the state of the object from a file.
    */

  void read(string runfile);

private:
  int Max_Gaussians;
  string Label;
  int N_Orbitals;
};

#endif

