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

#ifndef CambridgeThreeBodyCorrelationFunction_H
#define CambridgeThreeBodyCorrelationFunction_H

#include "QMCThreeBodyCorrelationFunction.h"
#include "Array3D.h"

/** 
   This is the three body correlation function described in 
   Phys Rev B 70, 235119 (2004)
*/

class CambridgeThreeBodyCorrelationFunction : public QMCThreeBodyCorrelationFunction
{
 protected:
  double FunctionValue;
  Array1D<double> grad1;
  Array1D<double> grad2;
  double LaplacianValue;
  
 private:

  Array3D<double> coeffs;

  int Nen, Nee, C;
  double cutoff;

 public:

  void initializeParameters(int electron_nucleus, int electron_electron, 
			    Array1D<double> &Parameters, int power, 
			    double max_dist);

  void evaluate(Array1D<double> &xyz1, double dist1, Array1D<double> &xyz2,
		double dist2);

  double getFunctionValue();
  double get_p_a(int ai);
  
  Array1D<double> * getElectron1Gradient();
  Array1D<double> * getElectron2Gradient();

  double get_p2_xa(int ai);
  
  double getLaplacianValue();
  double get_p3_xxa(int ai);

  double getCutoffDist();
};

#endif
