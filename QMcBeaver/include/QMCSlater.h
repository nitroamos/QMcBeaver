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
#include "QMCElectronNucleusCusp.h"
#include "svdcmp.h"

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
    Creates an uninitialized instance of the object.
  */
  QMCSlater();

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
  void initialize(QMCInput *input, int startEl, int stopEl, bool isAlpha);

  /**
    Evaluates the Slater determinants and their first two derivatives at X.
    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix

    This function only processes on the CPU.

    @param num how many configurations to process in the X array.
    @param start which index in X to start at
    @param X the array of electronic positions indexed by their walker
  */
  void evaluate(Array1D<Array2D<double>*> &X, int num, int start);

  /**
    Calling this function will wrap up the evaluation by unloading results
    from the GPU (if used) and then calculating the inverse.
  */
  void update_Ds();

  /**
    This method will have the inverse and derivative ratios calculated. This function
    was templated to allow flexibility in allowing both double and float without recompiling.

    @param start useful to indicate whether the GPU results have been filled in yet or not
    @param inv, grad, psi, lap: these are the collections of data
  */
  template <class T>
  void processInverse(int start,
      Array2D< Array2D<T> > & psi, Array2D< Array2D<T> > & inv, 
      Array2D< Array2D<T> > & lap, Array3D< Array2D<T> > & grad);

#ifdef QMC_GPU
  /*
    This function is pure GPU calculations. 

    @param num how many configurations to process in the X array.
    @param X the array of electronic positions indexed by their walker
  */
  void gpuEvaluate(Array1D<Array2D<double>*> &X, int num);
#endif


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
     Partial derivatives of the Slater determinants with respect to
     orbital coefficients.
  */
  Array2D<double> * get_p_a(int walker, int ci)
    {
      return  & (p_a(walker))(ci);
    }
  
  /**
     Partial derivative once with respect to orbital coefficient,
     and once with respect to position.
  */  
  Array2D<double> * get_p2_xa(int walker, int ci, int el, int xyz)
    {
      return & (p2_xa(walker))(ci,el,xyz);
    }

  /**
     Partial derivative once with respect to orbital coefficient,
     and twice with respect to position.
  */    
  Array2D<double> * get_p3_xxa(int walker, int ci)
    {
      return & (p3_xxa(walker))(ci);
    }

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

#ifdef QMC_GPU
  GPUQMCBasisFunction gpuBF;
  GPUQMCMatrix gpuMatMult;
#endif

 private:
  /**
     Whether this set of Slater determinants refers to
     alpha for beta electrons.
  */
  bool isAlpha;

  QMCInput *Input;
  QMCBasisFunction *BF;
  QMCWavefunction  *WF;
  
  /** 
    The dimensions of these arrays are numWalkers x numDeterminants
  */
  Array1D< Array1D<double> > Psi;
  Array1D< Array1D<double> > Laplacian_PsiRatio;
  Array1D< Array3D<double> > Grad_PsiRatio;
  Array1D< Array1D<double> > Chi_Density;

  /**
     Data structures to store the partial derivatives
     with respect to orbital coefficients.
  */
  Array1D< Array1D< Array2D<double> > > p_a;
  Array1D< Array3D< Array2D<double> > > p2_xa;
  Array1D< Array1D< Array2D<double> > > p3_xxa;

  Array1D< Array1D<bool> > Singular;

  /**
     The starting and stopping indices in the position
     array for the electrons this Slater is responsible for.
  */
  int Start;
  int Stop;

  Array1D< QMCElectronNucleusCusp > ElectronNucleusCusp;

  double PsiRatio_1electron;

  /** 
    The dimensions of these data are numWalkers x numDeterminants then 
    numElec x numOrb

    These data: D, D_inv, Laplacian_D, and Grad_D are meant to hold
    only the results that were calculated on the CPU
  */

  Array2D< Array2D<qmcfloat> > D;
  Array2D< Array2D<qmcfloat> > D_inv;
  Array2D< Array2D<qmcfloat> > Laplacian_D;
  Array3D< Array2D<qmcfloat> > Grad_D;

#ifdef QMC_GPU
  /** 
    The dimensions of these data are numWalkers x numDeterminants then 
    numElec x numOrb

    These data: D, D_inv, Laplacian_D, and Grad_D are meant to hold
    only the results that were calculated on the GPU
  */

  Array2D< Array2D<float> > gpu_D;
  Array2D< Array2D<float> > gpu_D_inv;
  Array2D< Array2D<float> > gpu_Laplacian_D;
  Array3D< Array2D<float> > gpu_Grad_D;

  /**
    This holds pointers to the GPU data. It is currently
    (unless i'm forgetting something i did) only useful for
    getting data from GPUQMCMatrix.
  */
  Array2D< Array2D<float>* > resultsCollector;
#endif

  /*
    The first dimension is the number of determinants.
    If we're optimizing orbitals, then we need to save all of these.

    If we are not optimizing, then they don't need to be kept after they're
    used in matrix multiplication.
  */
  Array1D< Array2D<qmcfloat> > Chi;
  Array1D< Array2D<qmcfloat> > Chi_laplacian;
  Array1D< Array1D< Array2D<qmcfloat> > > Chi_gradient;

  void allocate();

  /**
     Start and Stop are the indices in the electron coordinate
     array that this Slater determinant will include. Therefore,
     the number of electrons is Stop - Start + 1.
     @return the number of electrons for this Slater determinant
  */
  int getNumberElectrons()
    {
      return Stop-Start+1;
    }

  template <class T>
  void calculate_DerivativeRatios(int walker, int start, Array2D< Array2D<T> > & inv, 
      Array2D< Array2D<T> > & lap, Array3D< Array2D<T> > & grad);
};

#endif

