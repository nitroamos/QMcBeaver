#ifndef Williamson2CorrelationFunction_H
#define Williamson2CorrelationFunction_H

#include "FixedCuspPadeCorrelationFunction.h"

/** 

*/

class Williamson2CorrelationFunction : public FixedCuspPadeCorrelationFunction
{
 private: 
  double g, A, A2, F, A2F, s2g;

  double B;
  double r,x,L;

  //Some temporary variables
  double yuk, dyuk, d2yuk;
  double dyuk_a, d2yuk_a, d3yuk_a;
  double ir, t1, t2;

  double temper, dtemper, d2temper;

  double dG_a, dG_xa, dG_xxa;
  double pre, dpre, d2pre;
  
  int n;
  double xbar, dxbar, d2xbar;
  double sum_Tn, sum_dTn, sum_d2Tn;
  Array1D<double>   Tn;
  Array1D<double>  dTn;
  Array1D<double> d2Tn;
  Array1D<double>   co;
public:

  void initializeParameters(Array1D<int> & BeginningIndexOfParameterType, 
			    Array1D<double> &Parameters,
			    Array1D<int> & BeginningIndexOfConstantType, 
			    Array1D<double> & Constants);

  void evaluate( double r );
  double get_p_a(int ai);
  double get_p2_xa(int ai);
  double get_p3_xxa(int ai);

  void ChebyshevT_old(int n, double x, double &y, double &dy, double &d2y);
  void ChebyshevT(int n, double x, double &y, double &dy, double &d2y);
  bool isSingular();

  void print(ostream& strm);
};

#endif
