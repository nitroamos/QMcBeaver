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

#include "QMCInitializeWalkerFactory.h"

QMCInitializeWalker * QMCInitializeWalkerFactory::
  initializeWalkerFactory(QMCInput *input, string & Type)
{
  QMCInitializeWalker * initializeWalker = 0;

  if( Type == "mikes_jacked_initialization" )
    {
      initializeWalker = new QMCMikesJackedWalkerInitialization(input);
    }
   else if( Type == "mikes_better_initialization" )
    {
      initializeWalker = new QMCMikesBetterWalkerInitialization(input);
    }
  else
    {
      cerr << "ERROR: Unknown walker initialization type (" 
           << Type
           << ") being assigned in QMCInitializeWalkerFactory!" << endl;
      exit(0);
    }

  return initializeWalker;
}
