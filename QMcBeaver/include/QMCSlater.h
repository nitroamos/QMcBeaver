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

#ifndef QMCSLATER_H
#define QMCSLATER_H

#include <iostream>

#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "Array4D.h"
#include "LU.h"
#include "QMCInput.h"

using namespace std;

/** 
  An array of Slater determinants describing like spin electrons from a 3N 
  dimensional wavefunction.  This class allows the functions, their gradients, 
  and their laplacians to be calculated.
*/

class QMCSlater
{
public:
  /**
    Initializes the class and sets which region of the \f$3N\f$ dimensional 
    electronic configuration corresponds to electrons in these Slater 
    determinants.  It is assumed that all electrons in a determinant are 
    grouped together in the configuration.  
    
    @param input input data for the calculation
    @param startEl first particle in this determinant.
    @param stopEl last particle in this determinant.
  */
  void initialize(QMCInput *input, int startEl, int stopEl, Array2D<int> occ);

  /**
    Evaluates the slater determinants and their first two derivatives at X.
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
  */
  void evaluate( Array2D<double> &X);

  /**
    Gets an array of values of the Slater determinants for the last evaluated 
    electronic configuration.  The returned values are not normalized to one.  
    Assuming the basis functions ued to make the determinant are normalized, 
    the values can be normalized by dividing by \f$\sqrt{M!}\f$, where 
    \f$M\f$ is the number of electrons in the determinants.
  */
  Array1D<double>* getPsi();

  /**
    Gets an array where each element is a ratio of the Slater determinant 
    gradient over the Slater determinant for the last evaluated electronic 
    configuration.  These values do not depend on the normalization of the 
    Slater determinant.
  */
  Array3D<double>* getGradPsiRatio();

  /**
    Gets an array where each element is a ratio of the Slater determinant 
    laplacian over the Slater determinant for the last evaluated electronic 
    configuration.  These values do not depend on the normalization of the 
    Slater determinant.
  */
  Array1D<double>* getLaplacianPsiRatio();

  /**
    Returns true if the Slater determinant is singular and false otherwise.
  */
  bool isSingular();

  /**
    Sets two QMCSlater objects equal.
    @param rhs object to set this object equal to
    */
  void operator=(const QMCSlater & rhs );

 private:
  QMCInput *Input;
  QMCBasisFunction *BF;
  QMCWavefunction  *WF;
  
  Array1D<double> Psi;
  Array1D<double> Laplacian_PsiRatio;
  Array3D<double> Grad_PsiRatio;

  Array1D<bool> Singular;

  int Start;
  int Stop;
  Array2D<int> occupation;

  double PsiRatio_1electron;

  Array1D< Array2D<qmcfloat> > D;
  Array1D< Array2D<qmcfloat> > D_inv;
  Array1D< Array2D<qmcfloat> > Laplacian_D;
  Array2D< Array2D<qmcfloat> > Grad_D;

  // Scratch Space
  Array2D<qmcfloat> Chi;
  Array2D<qmcfloat> Chi_laplacian;
  Array1D< Array2D<qmcfloat> > Chi_gradient;

  Array2D<qmcfloat> WF_coeffs;

  void allocate(int N);

  /**
    Sets which region of the \f$3N\f$ dimensional electronic configuration 
    corresponds to electrons in this Slater determinant.  It is assumed 
    that all electrons in a determinant are grouped together in the 
    configuration.  
    
    @param startEl first particle in this determinant.
    @param stopEl last particle in this determinant.
  */
  void setStartAndStopElectronPositions(int startEl, int stopEl); 

  void update_Ds(Array2D<double> &X);

  void calculate_DerivativeRatios();

  void update_D_inverse_and_Psi();
};

#endif

