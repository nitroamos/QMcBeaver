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

  Array1D<double> p_a;
  Array2D<double> p2_x1a;
  Array2D<double> p2_x2a;
  Array1D<double> p3_xxa;

 private:

  Array3D<double> coeffs;

  int Nen, Nee, C;
  double cutoff;

  double pre1, pre2;
  double d1pre1, d1pre2;
  double d2pre1, d2pre2;

  double dU_dr1, dU_dr2, dU_dr12;

  double r1, r2, r12;
  Array1D<double> r1v, r2v, r12v;

 public:

  void initializeParameters(int electron_nucleus, int electron_electron, 
			    Array1D<double> &Parameters, int power, 
			    double max_dist);

  void evaluate(Array1D<double> &xyz1, double dist1,
		Array1D<double> &xyz2, double dist2,
		Array1D<double> &xyz12, double r12);

  double getFunctionValue();
  double get_p_a(int ai);

  double getLapPoly(double sign, double lterm, double nterm,
		    double l2term, double lnterm, double n2term,
		    double dpre, double pre, double dist);

  Array1D<double> * getElectron1Gradient();
  Array1D<double> * getElectron2Gradient();

  double get_p2_xa(bool e1, int xyz, int ai);
  
  double getLaplacianValue();
  double get_p3_xxa(int ai);

  double getCutoffDist();

  void print(ostream& strm);
};

#endif
