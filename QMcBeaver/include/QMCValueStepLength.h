
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
		    Array1D<double> & unused1,
		    Array1D<double> & unused2,
		    Array1D<double> & unused3,
		    Array2D<double> & unused4,
		    double functionValue)
    {
      return value;
    }
  
 private:
  const double value;
};

#endif
