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

#ifndef StringManipulation_H
#define StringManipulation_H

#include <stdio.h>
#include <string>
#include <strstream>

using namespace std;

/**
  A set of functions to manipulate strings.
  */

class StringManipulation
{
public:
  /**
    Converts a string to all upper case.

    @param s a string
    */

  static string toAllUpper(string & s);


  /** 
    Converts a string to all lower case. 

    @param s a string
   */

  static string toAllLower(string & s);

  
  /**
    Capitalizes the first letter and lowers all others in a string.

    @param s a string
    */
  
  static string toFirstUpperRestLower(string & s);

  /**
    Makes a character upper case.

    @param c a character
    */
  static char toUpperChar(char c);

  /**
    Makes a character lower case.

    @param c a character
    */
  static char toLowerChar(char c);
  

  /**
    Returns a string representation of an integer.

    @param i an integer.
    */
  static string intToString(int i);

  /**
    Returns a hexadecimal string representation of an integer.

    @param i an integer.
    */
  static string intToHexString(int i);

  /**
    Returns a string representation of a double.

    @param d a double.
    */
  static string doubleToString(double d);

  /**
    Returns a string representation of a long.

    @param l a long.
    */
  static string longToString(long l);

  /**
    Returns an int representation of a string.

    @param s a string.
    */
  static int stringToInt(string & s);

  /**
    Returns a long representation of a string.

    @param s a string.
    */
  static long stringToLong(string & s);

  /**
    Returns an representation of a hexadecimal string.

    @param s a string.
    */
  static int hexstringToInt(string & s);

  /**
    Returns an double representation of a string.

    @param s a string.
    */
  static double stringToDouble(string & s);
};

#endif


