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

#include "QMCSlater.h"
#include "QMCInput.h"
#include "QMCJastrow.h"
#include "QMCPotential_Energy.h"
#include "QMCGreensRatioComponent.h"

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
    controling the QMC calculation.

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
    Initializes the object with the data controling the QMC calculation.

    @param input input data for the calculation
  */
  void initialize(QMCInput *input); 

  /**
    Evaluates all of the calculated properties at X.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
  */
  void evaluate(Array2D<double> &X);

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

    @return potential energy.
  */
  double getPotentialEnergy();

  /**
    Gets the densities for the basis functions for the last evaluated
    electronic configuration.

    @return chi density.
  */
  Array1D<double>* getChiDensity();

  /**
    Gets the ratio of the wavefunction gradient to the wavefunction value at
    the last evaluated electronic configuration.  This is also known as the
    quantum force.

    @return wavefunction gradient ratio (quantum force)
  */
  Array2D<double>* getGradPsiRatio();

  /**
    Gets a modified version of the ratio of the wavefunction gradient to the 
    wavefunction value at the last evaluated electronic configuration.  
    The modifications typically help deal with singularities near nodes,
    and the particular type of modification can be selected.  
    This is also known as the modified quantum force.

    @return modified wavefunction gradient ratio (modified quantum force)
  */
  Array2D<double>* getModifiedGradPsiRatio();

  /**
    Returns true if the last evaluated electronic configuration gives a
    singular Slater matrix and false otherwise.

    @return true if the Slater matrix is singular and false otherwise
  */
  bool isSingular();

  /**
    Sets two QMCFunctions objects equal.

    @param rhs object to set this object equal to
  */
  void operator=(const QMCFunctions & rhs );

  /**
    Writes the state of this object to a stream for use in correlated
    sampling calculations.

    @param strm output stream
  */
  void writeCorrelatedSamplingConfiguration(ostream& strm);

 private:
  QMCInput *Input; 
 
  QMCSlater Alpha, Beta; 
  QMCJastrow Jastrow;  
  QMCPotential_Energy PE;

  QMCGreensRatioComponent Psi;
  double Laplacian_PsiRatio;
  Array2D<double> Grad_PsiRatio;
  Array2D<double> Modified_Grad_PsiRatio;

  Array1D<double> Chi_Density;

  double SCF_Laplacian_PsiRatio;
  Array2D<double> SCF_Grad_PsiRatio;

  double E_Local;

  void calculate_Psi_quantities();
  void calculate_Modified_Grad_PsiRatio(Array2D<double> &X);
  void calculate_E_Local();
};

#endif





