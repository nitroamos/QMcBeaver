
#ifndef QMCValueStepLength_H
#define QMCValueStepLength_H

#include <iostream>

#include "IeeeMath.h"
#include "QMCLineSearchStepLengthSelectionAlgorithm.h"

using namespace std;

class QMCValueStepLength : 
public QMCLineSearchStepLengthSelectionAlgorithm
{
public:
  QMCValueStepLength(double _value): value(_value)
    {

    }
  
  double stepLength(QMCObjectiveFunction *function, 
		    Array1D<double> & position,
		    Array1D<double> & searchDirection,
		    Array1D<double> & gradient,
		    double functionValue)
    {
      return value;
    }
  
 private:
  const double value;
};

#endif
