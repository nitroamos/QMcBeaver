#ifndef Umrigar2CorrelationFunction_H
#define Umrigar2CorrelationFunction_H

#include "FixedCuspPadeCorrelationFunction.h"

/** 
    This 2 body Jastrow function is based on the version
    described by Filippi and Umrigar in
    "Multiconfiguration wave functions for quantum Monte Carlo
    calculations of first-row diatomic molecules" in
    JCP 105, pg 213 1996
*/

class Umrigar2CorrelationFunction : public FixedCuspPadeCorrelationFunction
{
 private: 
  double U;
  double dU_dr, dU_drr;
  double dU_dk, dU_dkr, dU_dkrr;

  double k;
  double g;

  double a;
  double b, B, dB;

  double r;

  double den;
  double iden, iden2;

  bool optimizeG;
public:

  void initializeParameters(Array1D<int> & BeginningIndexOfParameterType, 
			    Array1D<double> &Parameters,
			    Array1D<int> & BeginningIndexOfConstantType, 
			    Array1D<double> & Constants);

  void evaluate( double r );
  double get_p_a(int ai);
  double get_p2_xa(int ai);
  double get_p3_xxa(int ai);

  bool isSingular();

  void print(ostream& strm);
};

#endif
