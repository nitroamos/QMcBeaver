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

#ifndef Exception_H
#define Exception_H

#include <string>

using namespace std;

/**
  An Exception is thrown when an error occurs.  This can be extended to
  deal with special types of errors.
  */

class Exception
{
private: string message;

/**
  Creates an exception.
 */
public: Exception();

/**
  Creates an exception.
 
  @param message A message describing what went wrong.
 */
public: Exception(string message);


/**
  Sets the error message for the exception.
  */
public: void setMessage(string message);

/**
  Gets the error message for the exception.
  */
public: string getMessage();
};


#endif
