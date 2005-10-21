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

#ifndef QMCCorrelationFunctionFactory_H
#define QMCCorrelationFunctionFactory_H

#include <string>

#include "QMCCorrelationFunction.h"
#include "ZeroCorrelationFunction.h"
#include "PadeCorrelationFunction.h"
#include "FixedCuspPadeCorrelationFunction.h"
#include "JuliusCorrelationFunction.h"

using namespace std;


/** 
  Object factory which returns the correct QMCCorrelationFunction when a 
  string keyword describing the correlation function is provided.
  */

class QMCCorrelationFunctionFactory
{
public:
  /** 
    Returns the correct QMCCorrelationFunction when a string keyword 
    describing the correlation function is provided.
    */

  static QMCCorrelationFunction * correlationFunctionFactory(string & Type);
};

#endif





