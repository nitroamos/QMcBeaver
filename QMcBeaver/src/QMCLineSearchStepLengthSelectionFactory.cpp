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

#include "QMCLineSearchStepLengthSelectionFactory.h"

QMCLineSearchStepLengthSelectionAlgorithm * 
  QMCLineSearchStepLengthSelectionFactory::factory(string & Type)
{
  QMCLineSearchStepLengthSelectionAlgorithm * algorithm = 0;

  if( Type == "MikesBracketing" )
    {
      algorithm = new QMCMikesBracketingStepLengthSelector();
    }
  else if( Type == "Wolfe" )
    {
      algorithm = new QMCWolfeStepLengthSelector();
    }
  else
    {
      cerr << "ERROR: Unknown line search step length selection algorithm (" 
           << Type
           << ") being assigned in QMCLineSearchStepLengthSelectionFactory!" 
	   << endl;
      exit(0);
    }

  return algorithm;
}
