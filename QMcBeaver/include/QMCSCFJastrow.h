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

#ifndef QMCSCFJastrow_H
#define QMCSCFJastrow_H

#include <iostream>
#include <sstream>

#include "QMCFunctions.h"
#include "QMCSlater.h"
#include "QMCInput.h"
#include "QMCJastrow.h"
#include "QMCPotential_Energy.h"
#include "QMCDouble.h"
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

class QMCSCFJastrow : public QMCFunctions
{
public:
  /**
    Creates a new instance of the class.
  */
  QMCSCFJastrow();

  /**
    Creates a new instance of the class and initializes it with the data 
    controlling the QMC calculation.

    @param input input data for the calculation
    @param HF object for calculating mean field potential
  */
  QMCSCFJastrow(QMCInput *input, QMCHartreeFock* HF);

  /**
    Creates a new instance of the class and initializes it with the data 
    controlling the QMC calculation.

    @param input input data for the calculation
  */
  QMCSCFJastrow(QMCInput *input);

  /**
    Creates a new instance of the class that is identical to another 
    instance of QMCSCFJastrow.

    @param rhs object to make a copy of
  */
  QMCSCFJastrow(const QMCSCFJastrow & rhs );

  /**
    Deallocates all memory used by the object.
  */
  ~QMCSCFJastrow();

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
    if true, the walkerData.configOutput will be given its info
  */
  void evaluate(Array2D<double> &X, QMCWalkerData & data);
  void evaluate(Array1D<QMCWalkerData *> &walkerData, 
		Array1D<Array2D<double> * > &xData, int num);

  /**
    Sets two QMCSCFJastrow objects equal.

    @param rhs object to set this object equal to
  */
  void operator=(const QMCSCFJastrow & rhs );

  int getNumTimers();
  void aggregateTimers(Array1D<Stopwatch> & timers,
		       int & idx);
 private:
  Array1D<Stopwatch> swTimers;

  QMCWalkerData * wd;
  Array2D<double> * x;
  int iWalker;

  /**
     Corresponding to QMCFunction's ability to process several walkers 
     simultaneously, each QMCSlater object is able to do the same in analogous
     fashion.
  */
  QMCSlater Alpha, Beta; 
  QMCJastrow Jastrow;  
  QMCPotential_Energy PE;
  QMCNuclearForces nf;

  QMCDouble Psi;

  double Laplacian_PsiRatio;

  double E_Local;

  /**
     These next several peices of data are the intermediary
     data used to calculate psi and the energy.

     Since they are likely to be needed by the parameter
     derivative functions, we need to make them available
     outside the function that assigns them.
  */
  Array1D<double> termPR;
  double** term_AlphaGrad;
  double** term_BetaGrad;
  
  Array1D<double>* alphaPsi;
  Array3D<double>* alphaGrad;
  Array1D<double>* alphaLaplacian;
  
  Array1D<double>* betaPsi;
  Array3D<double>* betaGrad;
  Array1D<double>* betaLaplacian;

  Array2D<double>* JastrowGrad;

  /**
     This function processes the results from a QMCSlater. The input parameters
     are the data structures this function is meant to fill and are from
     a QMCWalkerData struct.
  */
  void calculate_Psi_quantities();

  void update_SCF();

  /**
     This must be called after calculate_Psi_quantities.
     It assumes that we're only interested Jastrow parameter modifications,
     recalculate the energy for each set of parameters.
  */
  void calculate_CorrelatedSampling(Array1D<QMCWalkerData *> &walkerData,
				    Array1D<Array2D<double> * > &xData,
				    int num);

  /**
     This function continues to process results from a QMCSlater calculation.
     Grad_PsiRatio should be passed in as const, but doing so while passing
     it in with an & requires modification of several Array2D functions...
  */
  void calculate_Modified_Grad_PsiRatio();

  /**
     Calculate the derivative of the energy with respect to
     parameters in the Jastrow.

     @param ai the index where to start putting derivatives for
     this type
  */
  void calculate_JastrowDerivatives(int & ai);

  /**
     Calculate the derivative of the energy with respect to
     the CI coefficients.

     @param ai the index where to start putting derivatives for
     this type
  */
  void calculate_CIDerivatives(int & ai);

  /**
     Calculate the derivative of the energy with respect to
     parameters in the orbitals.

     @param ai the index where to start putting derivatives for
     this type
  */
  void calculate_OrbitalDerivatives(int & ai);

  /**
     Call this function to print out the parameter derivatives
     and then exit. The derivatives will be formatted in a way that it's easy to
     check them in a spreadsheet.

     To use this, uncomment the function call.
     Run the program and copy/paste the output into a spreadsheet.
     Change one of the parameters in the input file by a small
     amount "h".
     Rerun the program with the modified input file, and copy this
     data into your spreadsheet.
     Use the spreadsheet to calculate ( f(x+h) - f(x) ) / h
     and verify that this quantity matches the derivative the program
     calculated.

     The derivatives should match to a relative error of at least 1e-5,
     depending on the magnitude of h. If h is too small, then the quality
     of your derivative is going to suffer from poor output precision.
     If h is too large, then quality is going to suffer from the derivative
     approximation formula.
  */
  void checkParameterDerivatives();
            
  /**
     Gets the value of the wavefunction at the last evaluated electronic
     configuration.  The returned value is not normalized to one.
     
     @return wavefunction value
  */
  QMCDouble getPsi()
    {
      return Psi;
    }
  
  /**
     Gets the local energy at the last evaluated electronic configuration.

     E_L &=& -\frac{1}{2}\frac{\nabla^2 \left(D e^{U}\right)}{D e^{U}} + V \\
     &=&
     - \frac{1}{2} \left[\frac{\nabla^2 D}{D}
     +2 \left( \frac{\nabla D}{D} \right ) \cdot \nabla U
     + \nabla U \cdot \nabla U
     + \nabla^2 U \right]
     + V

     @return local energy
  */
  double getLocalEnergy()
    {
      return E_Local;
    }
  
  /**
     Gets the kinetic energy at the last evaluated electronic configuration.
     
     @return kinetic energy.
  */
  double getKineticEnergy()
    {
      return -0.5 * Laplacian_PsiRatio; 
    }
  
  /**
     Gets the potential energy at the last evaluated electronic configuration.
     @param which indicates which walker
     @return potential energy.
  */
  double getPotentialEnergy(int whichWalker)
    {
      return PE.getEnergy(whichWalker); 
    }
  
  double getEnergyEE(int whichWalker)
    {
      return PE.getEnergyEE(whichWalker);
    }
  
  double getEnergyNE(int whichWalker)
    {
      return PE.getEnergyNE(whichWalker);
    }
  
  /**
     Calculate the local energy.
     must calculate Laplacian_PsiRatio before calling this
     @param which indicates which walker
  */
  void calculate_E_Local(int whichWalker)
    {
      E_Local = -0.5 * Laplacian_PsiRatio + PE.getEnergy(whichWalker);
    }
  
  /**
     Returns true if the last evaluated electronic configuration gives a
     singular Slater matrix and false otherwise.
     
     @param walker indicate which QMCWalker we are interested in
     @return true if the Slater matrix is singular and false otherwise
  */
  bool isSingular(int whichWalker);
};

#endif





