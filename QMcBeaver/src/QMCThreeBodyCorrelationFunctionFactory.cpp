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

#include "QMCThreeBodyCorrelationFunctionFactory.h"
#include "CambridgeThreeBodyCorrelationFunction.h"
#include "ZeroThreeBodyCorrelationFunction.h"

QMCThreeBodyCorrelationFunction * 
  QMCThreeBodyCorrelationFunctionFactory::threeBodyCorrelationFunctionFactory(string & Type)
{
  QMCThreeBodyCorrelationFunction * ThreeBodyCorrelationFunction = 0;

  if( Type == "None" )
    ThreeBodyCorrelationFunction = new ZeroThreeBodyCorrelationFunction();
  else if( Type == "Cambridge" )
    ThreeBodyCorrelationFunction = new CambridgeThreeBodyCorrelationFunction();
  else
    {
      cerr << "ERROR: Unknown correlation function type (" << Type
           << ") being assigned in QMCThreeBodyCorrelationFunctionFactory!" 
	   << endl;
      exit(0);
    }

  return ThreeBodyCorrelationFunction;
}
