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

#ifndef ZeroThreeBodyCorrelationFunction_H
#define ZeroThreeBodyCorrelationFunction_H

#include "QMCThreeBodyCorrelationFunction.h"

/** 
   This is the three body correlation function that returns zeros- there is no
   interaction between the particles.
*/

class ZeroThreeBodyCorrelationFunction : public QMCThreeBodyCorrelationFunction
{
 public:

  void initializeParameters(int electron_nucleus, int electron_electron, 
			    Array1D<double> &Parameters, int C, double max);

  bool setElectron(bool first, Array1D<double> &xyz, double dist);

  void evaluate(Array1D<double> &xyz12, double r12);

  double getFunctionValue();
  double getFunctionValue(double r12, double r1, double r2);
  double get_p_a(int ai);
  
  Array1D<double> * getElectron1Gradient();
  Array1D<double> * getElectron2Gradient();

  double get_p2_xa(bool e1, int xyz, int ai);
  
  double getLaplacianValue();
  double get_p3_xxa(int ai);

  double getCutoffDist();

 private:

  Array1D<double> grad1;
  Array1D<double> grad2;
};

#endif
