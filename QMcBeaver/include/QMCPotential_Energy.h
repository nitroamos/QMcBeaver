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
  void initialize(QMCInput *input);


  /**
    Evaluates the potential energy for the given electronic configuration.

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */

  void evaluate(Array2D<double> &X);


  /**
    Gets the potential energy of the last configuration evaluated.
    */
  
  double getEnergy();


  /**
    Sets two QMCPotential_Energy objects equal.

    @param rhs object to set this object equal to
    */

  void operator=( const QMCPotential_Energy & rhs );

 private:
  QMCInput *Input;

  double Energy_total;

  double P_nn;
  double P_en;
  double P_ee;

  void calc_P_nn();
  void calc_P_en(Array2D<double> &R);
  void calc_P_ee(Array2D<double> &R);
  double rij(Array2D<double> &R, int i, int j);
  double rij(Array2D<double> &positioni,Array2D<double> &positionj, 
	     int i, int j);

};

#endif
