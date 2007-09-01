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

#ifndef QMCThreeBodyCorrelationFunctionFactory_H
#define QMCThreeBodyCorrelationFunctionFactory_H

#include <string>

#include "QMCThreeBodyCorrelationFunction.h"

using namespace std;

/** 
  Object factory which returns the correct QMCThreeBodyCorrelationFunction when
  a string keyword describing the correlation function is provided.
*/

class QMCThreeBodyCorrelationFunctionFactory
{
public:
  /** 
    Returns the correct QMCThreeBodyCorrelationFunction when a string keyword 
    describing the correlation function is provided.
  */

  static QMCThreeBodyCorrelationFunction * threeBodyCorrelationFunctionFactory(string & Type);
};

#endif





