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

#ifndef QMCInitializeWalkerFactory_H
#define QMCInitializeWalkerFactory_H

#include <string>

#include "QMCMikesJackedWalkerInitialization.h"
#include "QMCMikesBetterWalkerInitialization.h"
#include "QMCDansWalkerInitialization.h"
#include "QMCAmosBoringWalkerInitialization.h"

using namespace std;

/** 
  Object factory which returns the correct QMCInitialize walker when a 
  string keyword describing the correlation function is provided.
  */

class QMCInitializeWalkerFactory
{
public:

  /** 
    Returns the correct QMCInitializeWalker when a string keyword 
    describing the initialization method is provided.
    @param input input input data for the calculation
    @param type string describing which initialization algorithm to choose
    @return the selected QMCInitializeWalker method.
    */

  static QMCInitializeWalker * initializeWalkerFactory(QMCInput *input, 
						       string & type);
};

#endif





