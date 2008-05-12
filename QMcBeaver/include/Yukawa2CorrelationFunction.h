#ifndef Yukawa2CorrelationFunction_H
#define Yukawa2CorrelationFunction_H

#include "FixedCuspPadeCorrelationFunction.h"

/** 

*/

class Yukawa2CorrelationFunction : public FixedCuspPadeCorrelationFunction
{
 private: 
  double g, A, A2, F, A2F, s2g;
  
  double ir, t1, t2;

  double r;

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
