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

#ifndef QMCWALKERDATA_H
#define QMCWALKERDATA_H

#include <iostream>
#include <sstream>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCGreensRatioComponent.h"

using namespace std;

/**
   The QMCWalkerData data type is meant to hold all the information
   calculated by QMCFunction that is useful to QMCWalker. This should
   in effect decouple QMCFunction and QMCWalker from each other enabling
   QMCFunction to be treated a little bit differently without significant
   modifications to QMCWalker.
*/
struct QMCWalkerData {
  double localEnergy, kineticEnergy, potentialEnergy;
  double neEnergy, eeEnergy;

  QMCGreensRatioComponent psi;
  bool isSingular;
  
  /**
    Gets the ratio of the wavefunction gradient to the wavefunction value at
    the last evaluated electronic configuration.  This is also known as the
    quantum force.
  */
  Array2D<double> gradPsiRatio;

  /**
    Gets a modified version of the ratio of the wavefunction gradient to the 
    wavefunction value at the last evaluated electronic configuration.  
    The modifications typically help deal with singularities near nodes,
    and the particular type of modification can be selected.  
    This is also known as the modified quantum force.
  */
  Array2D<double> modifiedGradPsiRatio;

  /**
    This holds derivatives of the Hamiltonian (Hellman-Feynmann theorem)
    evaluated with respect to each nucleus.
  */  
  Array2D<double> nuclearDerivatives;
  
  Array1D<double> chiDensity;

  stringstream * configOutput;
};

#endif





