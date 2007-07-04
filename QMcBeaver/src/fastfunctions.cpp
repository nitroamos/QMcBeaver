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

double pythag(double a, double b)
{
  /*
    I copied this routine from:
    http://www3.telus.net/thothworks/EigRSvaloCPP0.html

    Thanks!
  */

  // Returns the square root of (a*a + b*b)
  // without overflow or destructive underflow

  double p, r, s, t, u;

  t = fabs(a);
  u = fabs(b);

  p = ((t >= u) ? t : u);
  if (p > 0){
    r = ((t <= u) ? t : u);
    r /= p;
    r *= r;
    t = 4.0 + r;

    while (t > 4.0){
      s = r/t;
      u = 1.0 + 2.0*s;
      p = u*p;
      t = s/u;
      r *= t*t;
      t = 4.0 + r;
    } // while (t > 4.0);
  } // End if (p > 0)

  return p;
} // End pythag
