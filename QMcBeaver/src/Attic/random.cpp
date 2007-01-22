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

#include <math.h>
#include "random.h"

double gasdev(long *idum)
{
  double ran1(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (*idum < 0) iset=0;
  if  (iset == 0) 
  {
    do 
    {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } 
    while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } 
  else 
  {
    iset=0;
    return gset;
  }
}
#define IA 16807
/* note #undef's at end of file */
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
  /*
    Notice that the use of static variables means
    that calling ran1 with the same seed will not result
    in the same random number. It has a MEMORY!!!
  */
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0 || !iy) 
  {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) 
    {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double expdev(long *idum)
{
  return -log(1.0-ran1(idum));
}

double sindev(long *idum)
{
  return acos(1.0-2.0*ran1(idum));
}

double randomDistribution1(long *idum)
{
  while(true)
    {
      double x = expdev(idum)/0.4;

      if( ran1(idum) < ( 0.5*x*x*exp(-x) )/( 0.8*exp(-0.4*x) ) )
	{
	  return x;
	}
    }
}



