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

#include "fastfunctions.h"

double fastPower(double x, int n)
{
  double temp = 0;
  switch(n)
    {
    case 0: 
      return 1.0;
    case 1:
      return x;
    case 2:
      return x*x;
    case 3:
      return x*x*x;
    case 4:
      temp = x*x;
      return temp * temp;
    case 5:
      temp = x*x;
      return temp * temp * x;
    case 6:
      temp = x*x*x;
      return temp * temp;
    default:
      return pow(x,n);
    }
}
