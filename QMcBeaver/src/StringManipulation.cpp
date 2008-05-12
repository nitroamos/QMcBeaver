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

#include "StringManipulation.h"
#include "math.h"

using namespace std;

// to upper case
string StringManipulation::toAllUpper(string & s)
{
  int size = (int)s.size();
  string result;

  for(int i=0; i<size; i++)
    {
      char upper = StringManipulation::toUpperChar(s[i]);
      result += upper;
    }

  return result;
}

// to lower case
string StringManipulation::toAllLower(string & s)
{
  int size = (int)s.size();
  string result;

  for(int i=0; i<size; i++)
    {
      char lower = StringManipulation::toLowerChar(s[i]);
      result += lower;
    }

  return result;
}


// capitialize first letter and lower all others
string StringManipulation::toFirstUpperRestLower(string & s)
{
  int size = (int)s.size();
  string result;

  if( size < 1 )
    {
      return result;
    }

  char upper = StringManipulation::toUpperChar(s[0]);
  result += upper;

  for(int i=1; i<size; i++)
    {
      char lower = StringManipulation::toLowerChar(s[i]);
      result += lower;
    }

  return result;
}

// Make a character upper case
char StringManipulation::toUpperChar(char c)
{
  switch( c )
    {
    case 'a': return 'A';
    case 'b': return 'B';
    case 'c': return 'C';
    case 'd': return 'D';
    case 'e': return 'E';
    case 'f': return 'F';
    case 'g': return 'G';
    case 'h': return 'H';
    case 'i': return 'I';
    case 'j': return 'J';
    case 'k': return 'K';
    case 'l': return 'L';
    case 'm': return 'M';
    case 'n': return 'N';
    case 'o': return 'O';
    case 'p': return 'P';
    case 'q': return 'Q';
    case 'r': return 'R';
    case 's': return 'S';
    case 't': return 'T';
    case 'u': return 'U';
    case 'v': return 'V';
    case 'w': return 'W';
    case 'x': return 'X';
    case 'y': return 'Y';
    case 'z': return 'Z';
    default: return c;
    }
}

// Make a character lower case
char StringManipulation::toLowerChar(char c)
{
  switch( c )
    {
    case 'A': return 'a';
    case 'B': return 'b';
    case 'C': return 'c';
    case 'D': return 'd';
    case 'E': return 'e';
    case 'F': return 'f';
    case 'G': return 'g';
    case 'H': return 'h';
    case 'I': return 'i';
    case 'J': return 'j';
    case 'K': return 'k';
    case 'L': return 'l';
    case 'M': return 'm';
    case 'N': return 'n';
    case 'O': return 'o';
    case 'P': return 'p';
    case 'Q': return 'q';
    case 'R': return 'r';
    case 'S': return 's';
    case 'T': return 't';
    case 'U': return 'u';
    case 'V': return 'v';
    case 'W': return 'w';
    case 'X': return 'x';
    case 'Y': return 'y';
    case 'Z': return 'z';
    default: return c;
    }
}

string StringManipulation::intToString(int i)
{
  //  ostringstream ostr;
  //  ostr << i << ends;
  //  return ostr.str();

  char charArray[64];
  snprintf(charArray,64,"%d",i);
  string result(charArray);
  return result;
}

string StringManipulation::intToHexString(int i)
{
  ostringstream ostr;
  ostr << hex << i << ends;
  return ostr.str();

  //  char charArray[64];
  //  snprintf(charArray,64,"%d",i);
  //  string result(charArray);
  //  return result;
}

string StringManipulation::doubleToString(double d)
{
  ostringstream ostr;
  ostr << d << ends;
  return ostr.str();

  //  char charArray[64];
  //  snprintf(charArray,64,"%f",d);
  //  string result(charArray);
  //  return result;
}

string StringManipulation::fancyDoubleToString(int sigFig, int width, double d)
{
  ostringstream ostr;
  ostr.precision(sigFig);
  ostr << showpos;
  string output;

  if(d == 0.0){
    ostr.precision(1);
    ostr.setf(ios::fixed);
    ostr << 0.0;
    output = ostr.str();
  } else if(fabs(d) > 1e-2)
    {
      ostr.setf(ios::fixed);
      ostr << d;
      output = ostr.str();
    } else {
      ostr.setf(ios::scientific);
      ostr << d;
      string temp = ostr.str();      
      string::size_type loc = temp.find("e", 0 );
      temp.replace(loc,1,"*^");
      output = temp;
    }

  int len = output.length();
  if(len < width)
    output.insert(0,width-len,' ');
  len = output.length();
  
  return output;
}

string StringManipulation::longToString(long l)
{
  ostringstream ostr;
  ostr << l << ends;
  return ostr.str();

  //  char charArray[64];
  //  snprintf(charArray,64,"%f",d);
  //  string result(charArray);
  //  return result;
}

int StringManipulation::stringToInt(string & s)
{
  istringstream istr(s, istringstream::in);
  int val;
  istr >> val;
  return val;
}

long StringManipulation::stringToLong(string & s)
{  
  istringstream istr(s, istringstream::in);
  long val;
  istr >> val;
  return val;
}

int StringManipulation::hexstringToInt(string & s)
{
  istringstream istr(s, istringstream::in);
  int val;
  istr >> hex >> val;
  return val;
}

double StringManipulation::stringToDouble(string & s)
{
  istringstream istr(s, istringstream::in);
  double val;
  istr >> val;
  return val;
}

