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


#ifndef XMLParseException_H
#define XMLParseException_H

#include <string>

#include "StringManipulation.h"
#include "Exception.h"


/**
  An XMLParseException is thrown when an error occures while parsing an XML
  stream.
 */

class XMLParseException : public Exception
{
/**
  Indicates that no line number has been associated with this exception.
  */
public: static int NO_LINE;


/**
  The line number in the source code where the error occurred, or
  <code>NO_LINE</code> if the line number is unknown.
  */
private: int lineNr;

/**
  Creates an exception.
  
  @param name The name of the element where the error is located.
  @param message A message describing what went wrong.
  */
public: XMLParseException(string name, string message);

/**
  Creates an exception.
  
  @param name The name of the element where the error is located.
  @param lineNr The number of the line in the input.
  @param message A message describing what went wrong.
  */
public: XMLParseException(string name, int lineNr, string message);

/**
  Where the error occurred, or <code>NO_LINE</code> if the line number is
  unknown.

  @return Line number where the error occurred.
  */
public: int getLineNr();
};

#endif
