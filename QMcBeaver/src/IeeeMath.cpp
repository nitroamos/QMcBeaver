//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2004.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   

Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/

#include "IeeeMath.h"

#define USING_QSC

template <class T>
bool IeeeMath::isNaN(T x)
{
  /**
    It seems difficult to get a portable version of this function.
    your mileage may vary!
    
    Hopefully this function won't give
    any more false positives.
  */
#ifdef USING_QSC
  if( x != x )
    {
      return true;
    }
  else
    {
      return false;
    }
#else
  if( x != x || !isfinite(x) )
    {
      return true;
    }
  else
    {
      return false;
    }
#endif
}

//QSC compiler quits if this line is uncommented...
#if ! defined USING_QSC
template bool IeeeMath::isNaN<double>(double x);
#endif

