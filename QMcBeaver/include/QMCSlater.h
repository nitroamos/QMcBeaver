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
#include "QMCInput.h"
#include "Stopwatch.h"
#include "QMCElectronNucleusCusp.h"
#include "svdcmp.h"
#include "QMCWalkerData.h"

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
    @param whichE the index of the electron to move (-1 if all)
  */
  void evaluate(Array1D<Array2D<double>*> & X,
		int num, int start, int whichE);

  /**
    Calling this function will wrap up the evaluation by unloading results
    from the GPU (if used) and then calculating the inverse.
    
    The new idea here is to split all the work that the evaluate function used
    to do.  This enables QMCSCFJastrow to control when the multiplication
    happens relative to the inverse.
  */
  void update_Ds(Array1D<QMCWalkerData *> & walkerData);

  /**
    This method will have the inverse and derivative ratios calculated. 

    @param ci which determinant we are working on
    @param row if moving one electron per iteration, this is the row
    the current electron belongs to
    
    @param psi (input) Slater matrix
    @param inv (output) inverse of psi
    @param lap (input) laplacian of Slater
    @param grad (input) gradient of Slater
    @param det (output) determinant of psi
    @param gradPR (output) gradient of psi ratio
    @param lapPR (output) laplacian of psi ratio
    @return whether the matrix was invertible
  */
  template <class T>
    bool calculate_DerivativeRatios(int ci, int row,
				    Array2D<T> & psi,
				    Array2D<double> & inv,
				    Array2D<T> & lap,
				    Array2D<T> & gradx,
				    Array2D<T> & grady,
				    Array2D<T> & gradz,
				    Array1D<double> & det,
				    Array3D<double> & gradPR,
				    Array1D<double> & lapPR);    
  
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
     These arrays are the local storage data that
     we use when we update all electrons at once.
  */
  Array1D< Array1D<double> > Psi;
  Array1D< Array1D<double> > Laplacian_PsiRatio;
  Array1D< Array3D<double> > Grad_PsiRatio;
  Array1D< Array1D<double> > Chi_Density;
  Array1D< Array1D<bool> >   Singular;

  /**
     These arrays are where we store the basis functions
     and their derivatives, evaluated for a particular
     electronic configuration.
  */
  Array1D<          Array2D<qmcfloat>   > Chi;
  Array1D<          Array2D<qmcfloat>   > Chi_laplacian;
  Array1D< Array1D< Array2D<qmcfloat> > > Chi_gradient;

  /**
    If we are updating one electron at a time, then
    these pointers need to be set to the data stored
    in the associated QMCWalkerData.

    If we are updating all together, then these pointers
    will be set to our local QMCSlater data.

    The reason we have two choices for storage location is
    because it will change how much memory we need. If we
    update all at once, then we don't need to save
    intermediate data between iterations, so we can just
    save everything right in QMCSlater, a class for which
    there are only 2 instances for the entire calculation.

    If we only update one at a time, then we need a
    per walker storage location.
  */
  Array1D< Array1D<double> * > pointer_det;
  Array1D< Array1D<double> * > pointer_lapPR;
  Array1D< Array3D<double> * > pointer_gradPR;

  /**
     Data structures to store the partial derivatives
     with respect to orbital coefficients.
  */
  Array1D< Array1D< Array2D<double> > > p_a;
  Array1D< Array3D< Array2D<double> > > p2_xa;
  Array1D< Array1D< Array2D<double> > > p3_xxa;

  /**
     The starting and stopping indices in the position
     array for the electrons this Slater is responsible for.
  */
  int Start;
  int Stop;

  QMCElectronNucleusCusp ElectronNucleusCusp;

  /** 
    The dimensions of these data are numWalkers x numDeterminants then 
    numElec x numOrb

    These data: D, D_inv, Laplacian_D, and Grad_D are meant to hold
    only the results that were calculated on the CPU

    D_inv is the only one that can be in double since we'll explicitly
    typecast D before inversion.
  */
  Array2D< double > * ciDet;
  Array1D< Array2D<qmcfloat> > D;
  Array2D< Array2D<qmcfloat> > Grad_D;
  Array1D< Array2D<qmcfloat> > Laplacian_D;
  Array2D< Array2D<double>   > D_inv;

#ifdef QMC_GPU
  /** 
    The dimensions of these data are numWalkers x numDeterminants then 
    numElec x numOrb

    These data: D, D_inv, Laplacian_D, and Grad_D are meant to hold
    only the results that were calculated on the GPU
  */

  Array1D< Array2D<float> > gpu_D;
  Array2D< Array2D<float> > gpu_D_inv;
  Array1D< Array2D<float> > gpu_Laplacian_D;
  Array2D< Array2D<float> > gpu_Grad_D;

  /**
    This holds pointers to the GPU data. It is currently
    (unless i'm forgetting something i did) only useful for
    getting data from GPUQMCMatrix.
  */
  Array2D< Array2D<float>* > resultsCollector;
#endif

  /*
    This function will make sure that
    all the data structures have the right dimensions
    for either matrix-matrix multiplication or
    vector-matrix multiplication.

    @param whichE the whichE index provided to evaluate
    @param start will be assigned which electron to start the update
    @param stop will be assigned which electron to stop the update
  */
  void allocateIteration(int whichE, int & start, int & stop);

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
};

#endif

