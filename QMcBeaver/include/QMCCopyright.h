//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef QMCCopyright_H
#define QMCCopyright_H

#include <string>
#include <iostream>

#include "StringManipulation.h"

using namespace std;

/**
  Central localtion for all copyright information relevant to QMcBeaver.
*/

class QMCCopyright
{
 public:
  /**
     Writes the copyright information to a stream in a human readable format.
  */
  friend ostream & operator<<(ostream & strm, QMCCopyright & rhs);
};

#endif
