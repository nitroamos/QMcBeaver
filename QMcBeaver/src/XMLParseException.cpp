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

#include "XMLParseException.h"

int XMLParseException::NO_LINE = -1;


XMLParseException::XMLParseException(string name, string message)
{
  string tempMessage = "XML Parse Exception during parsing of "
    + ( name.empty() ? "the XML definition" 
	: ("a " + name + " element")) + ": " + message;
  setMessage(tempMessage);

  this->lineNr = XMLParseException::NO_LINE;
}


XMLParseException::XMLParseException(string name, int lineNr, string message)
{
  string tempMessage = "XML Parse Exception during parsing of "
    + ( name.empty() ? "the XML definition"
	: ("a " + name + " element")) + " at line " 
    + StringManipulation::intToString(lineNr) + ": " + message;
  setMessage(tempMessage);

  this->lineNr = lineNr;
}


int XMLParseException::getLineNr()
{
  return this->lineNr;
}

