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

#include "Array1D.h"
#include "Array2D.h"
#include "QMCInput.h"
#include "QMCHartreeFock.h"
#include "QMCWalkerData.h"

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

  /**
     Calculates distance between electrons i and j
     based on the coordinates found in R
  */
  static double rij(Array2D<double> &R, int i, int j);

  /**
     Calculates the distance between particles i and j,
     where i is found in positioni and
     j is found in positionj.
  */
  static double rij(Array2D<double> &positioni,Array2D<double> &positionj, 
		    int i, int j);
  
 private:
  QMCInput *Input;
  QMCHartreeFock *HartreeFock;
  QMCWalkerData * wd;

  Array1D<double> Energy_total;
  Array1D<double> Energy_ne;
  Array1D<double> Energy_ee;

  double P_nn;
  double P_en;
  double P_ee;

  void calc_P_nn();
  void calc_P_en(Array2D<double> &R);
  void calc_P_ee(Array2D<double> &R);
};

#endif
