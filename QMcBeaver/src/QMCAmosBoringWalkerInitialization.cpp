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

#include "QMCAmosBoringWalkerInitialization.h"

Array2D<double> QMCAmosBoringWalkerInitialization::initializeWalkerPosition()
{
  int nelectrons = Input->WF.getNumberElectrons();
  int nbeta  = Input->WF.getNumberBetaElectrons();
  int nalpha = Input->WF.getNumberAlphaElectrons();
  int numDimensions = Input->WF.getNumberBasisFunctions();

  /**
     Why restrain ourselves to Cartesian coordinates?
  */
  Array2D<double> R(nelectrons,numDimensions);

  for(int i=0; i<R.dim1(); i++){
    for(int j=0; j<R.dim2(); j++){
      //R(i,j) = (ran1(&Input->flags.iseed) - 0.5)/3.0;
      R(i,j) = gasdev(&Input->flags.iseed);
    }
  }
  return R;
}















