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

#ifndef QMCFunctionsFactory_H
#define QMCFunctionsFactory_H

#include <string>

#include "QMCFunctions.h"
#include "QMCInput.h"
#include "QMCSCFJastrow.h"
#include "QMCHarmonicOscillator.h"

using namespace std;


/** 
    Object factory which returns the correct QMCFunctions when a 
    string keyword describing the correlation function is provided.
*/

class QMCFunctionsFactory
{
public:
  /** 
      Returns the correct QMCFunctions when a string keyword 
      describing the correlation function is provided.
  */
  
  static QMCFunctions * functionsFactory(QMCInput * input, string & Type);
};

#endif





