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

#ifndef QMCLineSearchStepLengthSelectionFactory_H
#define QMCLineSearchStepLengthSelectionFactory_H

#include <string>

#include "QMCLineSearchStepLengthSelectionAlgorithm.h"
#include "QMCMikesBracketingStepLengthSelector.h"
#include "QMCWolfeStepLengthSelector.h"
#include "QMCValueStepLength.h"

using namespace std;

/** 
  Object factory which returns the correct 
  QMCLineSearchStepLengthSelectionAlgorithm when a 
  string keyword describing the correlation function is provided.
  */

class QMCLineSearchStepLengthSelectionFactory
{
public:
  /** 
    Returns the correct QMCLineSearchStepLengthSelectionAlgorithm when a 
    string keyword describing the correlation function is provided.
    */

  static QMCLineSearchStepLengthSelectionAlgorithm * factory(string & Type);
};

#endif





