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

#ifndef QMCPotential_Energy_H
#define QMCPotential_Energy_H

#include <iostream>

#include "Stopwatch.h"
#include "Array1D.h"
#include "Array2D.h"
#include "QMCInput.h"
#include "QMCHartreeFock.h"
#include "QMCWalkerData.h"
#include "QMCElectronNucleusCusp.h"

using namespace std;


/**
  The potential energy of the system.
  */

class QMCPotential_Energy
{
public:
  /**
    Creates an instance of the class.
    */
  QMCPotential_Energy();


  /**
    Initialize the object.

    @param input data input to control the calculation
    */
  void initialize(QMCInput *input, QMCHartreeFock *HF);

  /**
     If pseudopotentials are used, call this function to evaluate it.

     @param R the list of electron positions
     @param elec the electron for which the ecp is evaluated
     @param nuc the nucleus with the ecp
  */
  double evaluatePseudoPotential(Array2D<double> & R, int elec, int nuc);

  /**
    Evaluates the potential energy for the given electronic configuration.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */
  void evaluate(Array1D<Array2D<double>*> &X,
		Array1D<QMCWalkerData *> &walkerData,
		int num);

  /**
    Gets the potential energy of the last configuration evaluated.
  */
  double getEnergy(int which);

  /**
    Gets the nuc-elec potential energy of the last config evaluated.
  */
  double getEnergyNE(int which);

  /**
    Gets the elec-elec potential energy of the last config evaluated.
  */
  double getEnergyEE(int which);

  /**
    Sets two QMCPotential_Energy objects equal.

    @param rhs object to set this object equal to
  */
  void operator=( const QMCPotential_Energy & rhs );

  int getNumTimers();
  void aggregateTimers(Array1D<Stopwatch> & timers,
		       int & idx);

  /**
     This variable is used by QMCRun to control printing of debug data.
  */
  static int printElec;
 private:
  Array1D<Stopwatch> swTimers;

  QMCInput *Input;
  QMCHartreeFock *HartreeFock;
  QMCWalkerData * wd;
  QMCWavefunction * WF;
  QMCMolecule * MOL;

  Array1D<double> Energy_total;
  Array1D<double> Energy_ne;
  Array1D<double> Energy_ee;

  //A bunch of data for ecp evaluation
  Array2D< Array2D< Array1D<double> > > angularGrids;
  Array1D<double> integrand;
  Array2D<qmcfloat> * coeffs;
  Array2D<double> X;
  Array2D<double> D;
  Array2D<double> ciDet;
  Array1D< Array2D<double> >   Dc_inv;
  QMCElectronNucleusCusp ElectronNucleusCuspA;
  QMCElectronNucleusCusp ElectronNucleusCuspB;

  double P_nn;
  double P_en;
  double P_ee;

  void calc_P_nn();
  void calc_P_en(Array2D<double> &R);
  void calc_P_ee(Array2D<double> &R);
  
  /**
     This function prints out the values of the local & non local
     pseudopotential energies
  */
  void printPseudoPotential(double max, int num, int nuc);
};

#endif
