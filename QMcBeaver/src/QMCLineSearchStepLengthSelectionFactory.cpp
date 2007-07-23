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
#include "QMCMikesBracketingStepLengthSelector.h"
#include "QMCWolfeStepLengthSelector.h"
#include "QMCValueStepLength.h"
#include "QMCLinearizeStepLength.h"

QMCLineSearchStepLengthSelectionAlgorithm * 
  QMCLineSearchStepLengthSelectionFactory::factory(string & Type)
{
  QMCLineSearchStepLengthSelectionAlgorithm * algorithm = 0;

  bool isValue = false;
  double value = 0;

  char ** ptr = 0;
  value = strtod(Type.c_str(), ptr);

  if( Type == "MikesBracketing" )
    {
      algorithm = new QMCMikesBracketingStepLengthSelector();
    }
  else if( Type == "Wolfe" )
    {
      algorithm = new QMCWolfeStepLengthSelector();
    }
  else if( Type == "Linearize" )
    {
      algorithm = new QMCLinearizeStepLength();
    }
  else if( Type == "None" )
    {
      algorithm = new QMCValueStepLength(1.0);
    }
  else if( fabs(value) > 1.0e-10 )
    {
      algorithm = new QMCValueStepLength(value);
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
