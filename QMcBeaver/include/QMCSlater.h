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
#include "Stopwatch.h"

#ifdef QMC_GPU
#include "GPUQMCMatrix.h"
#include "GPUQMCBasisFunction.h"
#endif

using namespace std;

/** 
  An array of Slater determinants describing like spin electrons from a 3N 
  dimensional wavefunction.  This class allows the functions, their gradients, 
  and their laplacians to be calculated.

  This class has now been modified to handle several electron configurations
  simultaneously. It will process (at most) WALKERS_PER_PASS walkers.
*/

class QMCSlater
{
public:
  /**
    Deallocates all memory used by the object.
  */
  ~QMCSlater();

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
    Evaluates the Slater determinants and their first two derivatives at X.
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    @param num how many configurations to process in the X array.
  */
  void evaluate(int num);

  /**
    An array of configurations are processed simultaneously. Call this function
    to set up the basis function and matrix multiplication work, and then
    when the results are actually needed, THEN call evaluate for psi et al
    to be filled in.
  */
  void update_Ds(Array1D<Array2D<double>*> &X, int num);

  /**
    Gets an array of values of the Slater determinants for the last evaluated 
    electronic configuration.  The returned values are not normalized to one.  
    Assuming the basis functions ued to make the determinant are normalized, 
    the values can be normalized by dividing by \f$\sqrt{M!}\f$, where 
    \f$M\f$ is the number of electrons in the determinants.

    @param i of which walker we are requesting the information
  */
  Array1D<double>* getPsi(int i);

  /**
    Gets an array where each element is a ratio of the Slater determinant 
    gradient over the Slater determinant for the last evaluated electronic 
    configuration.  These values do not depend on the normalization of the 
    Slater determinant.

    @param i of which walker we are requesting the information
  */
  Array3D<double>* getGradPsiRatio(int i);

  /**
    Gets an array where each element is a ratio of the Slater determinant 
    laplacian over the Slater determinant for the last evaluated electronic 
    configuration.  These values do not depend on the normalization of the 
    Slater determinant.

    @param i of which walker we are requesting the information
  */
  Array1D<double>* getLaplacianPsiRatio(int i);

  /**
    Gets an array of the densities for the basis functions for the last 
    evaluated electronic configuration.

    @param i of which walker we are requesting the information
  */
  Array1D<double>* getChiDensity(int i);

  /**
     Returns true if the Slater determinant is singular and false otherwise.
     
     @param i of which walker we are requesting the information
  */
  bool isSingular(int i);

  /**
     Sets two QMCSlater objects equal.
     @param rhs object to set this object equal to
  */
  void operator=(const QMCSlater & rhs );

 private:
  QMCInput *Input;
  QMCBasisFunction *BF;
  QMCWavefunction  *WF;
  
  /* The dimensions of these data are numWalkers x numDeterminants*/
  Array1D< Array1D<double> > Psi;
  Array1D< Array1D<double> > Laplacian_PsiRatio;
  Array1D< Array3D<double> > Grad_PsiRatio;
  Array1D< Array1D<double> > Chi_Density;

  Array1D< Array1D<bool> > Singular;

  int Start;
  int Stop;
  Array2D<int> occupation;

  double PsiRatio_1electron;

  /* The dimensions of these data are numWalkers x numDeterminants then numElec x numOrb*/
  Array2D< Array2D<qmcfloat> > D;
  Array2D< Array2D<qmcfloat> > D_inv;
  Array2D< Array2D<qmcfloat> > Laplacian_D;
  Array3D< Array2D<qmcfloat> > Grad_D;
#ifdef QMC_GPU
  GPUQMCBasisFunction gpuBF;
  GPUQMCMatrix gpuMatMult;
  // Scratch Space
  Array1D< Array2D<qmcfloat> > bfData;
  Array2D< Array2D<qmcfloat>* > resultsCollector;
#else
  // Scratch Space
  Array2D<qmcfloat> Chi;
  Array2D<qmcfloat> Chi_laplacian;
  Array1D< Array2D<qmcfloat> > Chi_gradient;
#endif

  /* The dimensions of WF_coeffs is numDeterminants and then numBasisFunction x numOrb */
  Array1D< Array2D<qmcfloat> > WF_coeffs;

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

  void calculate_DerivativeRatios(int walker);
};

#endif

