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

#include "QMCFunctionsFactory.h"

QMCFunctions * 
  QMCFunctionsFactory::functionsFactory(QMCInput * input, string & Type)
{
  QMCFunctions * Functions = 0;

  if( Type == "restricted" || Type == "unrestricted")
    {
      Functions = new QMCSCFJastrow();
    }
  else if( Type == "harmonicoscillator" )
    {
      Functions = new QMCHarmonicOscillator(input);
    }
  else
    {
      cerr << "ERROR: Unknown wavefunction type (" 
           << Type
           << ") being assigned in QMCFunctionsFactory!" << endl;
      exit(0);
    }

  return Functions;
}
