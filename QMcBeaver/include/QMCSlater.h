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
#include "LU.h"
#include "QMCInput.h"

using namespace std;


/** 
  A Slater determinant describing like spin electrons from a 3N dimensional 
  wavefunction.  This class allows the function, it's gradient, and it's 
  laplacian to be calculated.
  */

class QMCSlater
{
public:
  /**
    Initializes the class and 
    sets which region of the \f$3N\f$ dimensional electronic configuration 
    corresponds to electrons in this Slater determinant.  It is assumed 
    that all electrons in a determinant are grouped together in the 
    configuration.  
    
    @param input input data for the calculation
    @param startEl first particle in this determinant.
    @param stopEl last particle in this determinant.
    */

  void initialize(QMCInput *input, int startEl, int stopEl);

  /**
    Evaluates the slater determinant and it's first two derivatives at
    X.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */

  void evaluate( Array2D<double> &X);


  /**
    Gets the value of the Slater determinant for the last evaluated 
    electronic configuration.  The returned value is not normalized to one.  
    Assuming the basis functions ued to make the determinant are normalized, 
    this value can be normalized by dividing it by \f$\sqrt{M!}\f$, where 
    \f$M\f$ is the number of electrons in this determinant.
    */

  double getPsi();


  /**
    Gets the ratio of the Slater determinant gradient over the Slater 
    determinant for the last evaluated electronic configuration.  This value 
    does not depend on the normalization of the Slater determinant.
    */

  Array2D<double> * getGradPsiRatio();


  /**
    Gets the ratio of the Slater determinant laplacian over the Slater 
    determinant for the last evaluated electronic configuration.  This value 
    does not depend on the normalization of the Slater determinant.
    */

  double getLaplacianPsiRatio();


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

  double Psi;
  double Laplacian_PsiRatio;
  Array2D<double> Grad_PsiRatio;
  bool Singular;

  int Start;
  int Stop;

  double PsiRatio_1electron;

  Array2D <double> D;
  Array2D <double> D_inv;
  Array2D <double> Laplacian_D;
  Array3D <double> Grad_D;

  // Scratch Space
  Array1D <double> Chi1D;
  Array2D <double> Chi2D;
  Array1D <double> Grad1e;

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


  void initialize_D(Array2D<double> &X);
  void initialize_Laplacian_D(Array2D<double> &X);
  void initialize_Grad_D(Array2D<double> &X);

  void update_Ds(int electron, Array2D<double> &X);
  void update_D(int electron, Array2D<double> &X);
  void update_Laplacian_D(int electron, Array2D<double> &X);
  void update_Grad_D(int electron, Array2D<double> &X);

  void calculate_DerivativeRatios();
  void calculate_Laplacian_PsiRatio();
  void calculate_Grad_PsiRatio();

  void update_D_inverse_and_Psi();

};

#endif

