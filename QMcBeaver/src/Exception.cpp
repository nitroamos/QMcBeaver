//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2002.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "Exception.h"

Exception::Exception()
{
}

Exception::Exception(string message)
{
  setMessage(message);
}

void Exception::setMessage(string message)
{
  this->message = message;
}

string Exception::getMessage()
{
  return message;
}
