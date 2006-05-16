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

#ifndef QMCFUNCTIONS_H
#define QMCFUNCTIONS_H

#include <iostream>
#include <sstream>

#include "QMCSlater.h"
#include "QMCInput.h"
#include "QMCJastrow.h"
#include "QMCPotential_Energy.h"
#include "QMCGreensRatioComponent.h"
#include "QMCHartreeFock.h"
#include "QMCConfigIO.h"
#include "IeeeMath.h"
#include "QMCNuclearForces.h"
#include "QMCWalkerData.h"

using namespace std;

/**
  This class calculates the value of the wavefunction, its first two 
  derivatives, and any other properties which are calculated from the 
  wavefunction (local energy, etc.).  

  The wavefunction is assumed to be of the form
  \f[
  \Psi_{QMC} = J\sum_{i}D_{\uparrow,i}D_{\downarrow,i}
  \f]
  where
  \f[
  J=exp(\sum{u_{i,j}(r_{i,j})})
  \f]
  is a Jastrow type correlation function,
  and \f$D_{\uparrow}\f$ and \f$D_{\downarrow}\f$ are Slater determinants
  for the up and down electrons respectively.

  It has now been modified to treat several configurations simultaneously.
  The point of this is so that a QMCSlater has more of a SIMD nature than
  it had previously, and that if benefits can be derived by multiplying
  several matricies simultaneously, it is now much easier to do so. This
  change anticipates GPU (and similiar architecture...) modifications.
*/

class QMCFunctions
{
public:
  /**
    Creates a new instance of the class.
  */
  QMCFunctions();

  /**
    Creates a new instance of the class and initializes it with the data 
    controlling the QMC calculation.

    @param input input data for the calculation
    @param HF object for calculating mean field potential
  */
  QMCFunctions(QMCInput *input, QMCHartreeFock* HF);

  /**
    Creates a new instance of the class and initializes it with the data 
    controlling the QMC calculation.

    @param input input data for the calculation
  */
  QMCFunctions(QMCInput *input);

  /**
    Creates a new instance of the class that is identical to another 
    instance of QMCFunctions.

    @param rhs object to make a copy of
  */
  QMCFunctions(const QMCFunctions & rhs );

  /**
    Deallocates all memory used by the object.
  */
  ~QMCFunctions();

  /**
    Initializes the object with the data controlling the QMC calculation.

    @param input input data for the calculation
  */
  void initialize(QMCInput *input, QMCHartreeFock *HF); 

  /**
    Evaluates all of the calculated properties at X and places the calculated
    data into the QMCWalkerData struct provided. Two overloaded functions are
    provided, one of them processes a array of parameters, the other processes
    just one (useful during a QMCWalker's initialization)

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    @param data all the data that a QMCWalker should require
    @param writeConfig if the program is writing configs, we need to know here.
    if true, the walkerData.configOutput will be given its info
  */
  void evaluate(Array2D<double> &X, QMCWalkerData & data);
  void evaluate(Array1D<QMCWalkerData *> &walkerData, 
		Array1D<Array2D<double> * > &xData, int num, bool writeConfig);

  /**
    Sets two QMCFunctions objects equal.

    @param rhs object to set this object equal to
  */
  void operator=(const QMCFunctions & rhs );

 private:
  QMCInput *Input; 
 
  int nalpha, nbeta;

  /**
     Corresponding to QMCFunction's ability to process several walkers 
     simultaneously, each QMCSlater object is able to do the same in analogous
     fashion.
  */
  QMCSlater Alpha, Beta; 
  QMCJastrow Jastrow;  
  QMCPotential_Energy PE;
  QMCNuclearForces nf;

  QMCGreensRatioComponent Psi;
  double Laplacian_PsiRatio;

  double SCF_Laplacian_PsiRatio;
  Array2D<double> SCF_Grad_PsiRatio;

  double E_Local;

  /**
     This function processes the results from a QMCSlater. The input parameters
     are the data structures this function is meant to fill and are from
     a QMCWalkerData struct.
  */
  void calculate_Psi_quantities(Array2D<double> & Grad_PsiRatio,
				Array1D<double> & Chi_Density, int walker);

  /**
     This function continues to process results from a QMCSlater calculation.
     Grad_PsiRatio should be passed in as const, but doing so while passing
     it in with an & requires modification of several Array2D functions...
  */
  void calculate_Modified_Grad_PsiRatio(Array2D<double> & X, 
					Array2D<double> & Modified_Grad_PsiRatio, 
					Array2D<double> & Grad_PsiRatio);
          
  /**
    Calculate the local energy.
    @param which indicates which walker
    */
  void calculate_E_Local(int which);

  /**
    Gets the value of the wavefunction at the last evaluated electronic
    configuration.  The returned value is not normalized to one.

    @return wavefunction value
  */
  QMCGreensRatioComponent getPsi();

  /**
    Gets the local energy at the last evaluated electronic configuration.

    @return local energy
  */
  double getLocalEnergy();

  /**
    Gets the kinetic energy at the last evaluated electronic configuration.

    @return kinetic energy.
  */
  double getKineticEnergy();

  /**
    Gets the potential energy at the last evaluated electronic configuration.
    @param which indicates which walker
    @return potential energy.
  */
  double getPotentialEnergy(int which);

  double getEnergyEE(int which);
  double getEnergyNE(int which);

  /**
    Returns true if the last evaluated electronic configuration gives a
    singular Slater matrix and false otherwise.

    @param walker indicate which QMCWalker we are interested in
    @return true if the Slater matrix is singular and false otherwise
  */
  bool isSingular(int walker);
};

#endif





