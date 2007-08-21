
#ifndef QMCLinearizeStepLength_H
#define QMCLinearizeStepLength_H

#include <iostream>

#include "IeeeMath.h"
#include "QMCLineSearchStepLengthSelectionAlgorithm.h"

using namespace std;

class QMCLinearizeStepLength : 
public QMCLineSearchStepLengthSelectionAlgorithm
{
public:
  QMCLinearizeStepLength()
    {

    }

  bool isLinear(int ai);

  /**
     This method is from the
     JCP 126 084102 (2007)
     paper. There are a couple variations on this method,
     specifically on how to choose ksi. I have not fully
     studied all the possibilities, and I'm not even 100%
     sure that i programmed this correctly.

     This is designed to go with the "generalized_eigenvector"
     optimization criteria.
  */
  double stepLength(QMCObjectiveFunction *function, 
		    Array1D<double> & delta_x,
		    Array1D<double> & unused1,
		    Array1D<double> & unused2,
		    Array2D<double> & overlap,
		    double ksi);  
 private:

};

#endif
