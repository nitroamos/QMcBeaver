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

#include "QMCCorrelationFunctionFactory.h"
#include "Cambridge2CorrelationFunction.h"
#include "Umrigar2CorrelationFunction.h"
#include "Yukawa2CorrelationFunction.h"
#include "Williamson2CorrelationFunction.h"
#include "Anderson2CorrelationFunction.h"

QMCCorrelationFunction * 
  QMCCorrelationFunctionFactory::correlationFunctionFactory(string & Type)
{
  QMCCorrelationFunction * CorrelationFunction = 0;
  if( Type == "None" )
    {
      CorrelationFunction = new ZeroCorrelationFunction();
    }
  else if( Type == "Pade" )
    {
      CorrelationFunction = new PadeCorrelationFunction();
    }
  else if( Type == "FixedCuspPade" )
    {
      CorrelationFunction = new FixedCuspPadeCorrelationFunction();
    }
  else if( Type == "Cambridge2" )
    {
      CorrelationFunction = new Cambridge2CorrelationFunction();
    }
  else if( Type == "Umrigar2" )
    {
      CorrelationFunction = new Umrigar2CorrelationFunction();
    }
  else if( Type == "Yukawa" )
    {
      CorrelationFunction = new Yukawa2CorrelationFunction();
    }
  else if( Type == "Williamson" )
    {
      CorrelationFunction = new Williamson2CorrelationFunction();
    }
  else if( Type == "Anderson" )
    {
      CorrelationFunction = new Anderson2CorrelationFunction();
    }
  else if( Type == "Julius" )
    {
      CorrelationFunction = new JuliusCorrelationFunction();
    }
  else
    {
      cerr << "ERROR: Unknown correlation function type (" 
           << Type
           << ") being assigned in QMCCorrelationFunctionFactory!" << endl;
      exit(0);
    }

  return CorrelationFunction;
}
