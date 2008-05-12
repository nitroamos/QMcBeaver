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

#include "Random.h"

Random::Random()
{
#ifdef USESPRNG
  stream = 0;
#endif
}

Random::Random(long seed)
{
#ifdef USESPRNG
  stream = 0;
#endif
  initialize(seed,0);
}

Random::Random(const Random & rhs)
{
  start = rhs.start;
  current = rhs.current;

#ifdef USESPRNG
  //Would we ever need more than one generator?
  stream = rhs.stream;
#endif

  iy = rhs.iy;
  for(int i=0; i<NTAB; i++)
    iv[i] = rhs.iv[i];
}

Random::~Random()
{
#ifdef USESPRNG
  if(stream != 0)
    stream->free_sprng();
#endif
}

void Random::initialize(long seed, int rank)
{
  if(seed == 0)
    {
      srand( time(NULL) );
      current = -rand();
      seed = intdev();
      reset();
      //intdev returns a postive number, but we'll be switching the sign
      clog << "Using iseed = -" << seed << endl;
    }

#ifdef USESPRNG
  int sprng_seed = (int)(seed);

  //The sign bit probably doesn't matter...
  if(sprng_seed < 0)
    sprng_seed *= -1;
  start = (long)(sprng_seed);

  /*
    case 0: return new LFG;
    Combined Multiple Recursive Generator

    case 1: return new LCG;
    48 Bit Linear Congruential Generator with Prime Addend

    case 2: return new LCG64;
    64 Bit Linear Congruential Generator with Prime Addend

    case 3: return new CMRG;
    Modified Lagged Fibonacci Generator

    case 4: return new MLFG;
    Multiplicative Lagged Fibonacci Generator

    case 5: return new PMLCG;
    Prime Modulus Linear Congruential Generator
    (not installed in SPRNG by default)
  */
  if(stream == 0)
    stream = SelectType(0);

  int numStreams = 1;

#ifdef PARALLEL
  MPI_Comm_size(MPI_COMM_WORLD, &numStreams);
#endif

  clog << "Segment fault next line? Something is wrong with SPRNG\n";
  if(! stream->init_sprng(rank,numStreams,sprng_seed,SPRNG_DEFAULT) )
    {
      cerr << "Error: SPRNG initialization failed\n";
      exit(0);
    }

#else
  // ran1 needs to be initialized with a negative number
  if(seed > 0) seed *= -1;
  start = seed;
  current = seed;

  int my_seed = current;
  for(int i=0; i<rank; i++)
    my_seed = intdev();
  current = -1*abs(my_seed);
  reset();
#endif
}

void Random::reset()
{
#ifdef USESPRNG

#else
  //resetting ran1's internal static variables
  iy = 0;
  for(int i=0; i<NTAB; i++)
    iv[i] = 0;
#endif
}
void Random::printStream(ostream & strm)
{
  strm << "Random stream info:" << endl;
#ifdef USESPRNG
  //Apparently I can't choose a stream
  stream->print_sprng();
#else
  strm << "current ran1 seed = " << current << endl;
#endif
}

void Random::writeXML(ostream & strm)
{
#ifdef USESPRNG
  char * buffer;
  int bytes_required = stream->pack_sprng(&buffer);
  strm << "<sprng>\n" << *buffer << "\n</sprng>" << endl;
#else
  //Note: this doesn't fully save the state of ran1
  strm << "<iseed>\n" << current << "\n</iseed>" << endl;
  strm << "<internal>\n";
  strm << iy << " ";
  for(int i=0; i<NTAB; i++)
    strm << iv[i] << " ";
  strm << "\n</internal>" << endl;
#endif
}

bool Random::readXML(istream & strm)
{
#ifdef USESPRNG
  clog << "readXML in Random.cpp needs to be debugged!\n";
  int size;
  strm >> size;
  char temp[size];

  strm >> temp;  
  if (temp != "<sprng>")
    return false;
  strm >> temp;
  if(! stream->unpack_sprng(temp) )
    {
      cerr << "Error: Unpacking SPRNG stream failed\n";
      exit(0);
    }
  strm >> temp;
  if (temp != "</sprng>")
    return false;

  return true;

#else

  string temp;
  strm >> temp;
  if (temp != "<iseed>")
    return false;
  strm >> temp;
  long iseed = atoi( temp.c_str() );
  current = iseed;
  //initialize(iseed,0);
  strm >> temp;
  if (temp != "</iseed>")
    return false;

  strm >> temp;
  strm >> iy;
  for(int i=0; i<NTAB; i++)
    strm >> iv[i];
  strm >> temp;
  if(temp != "</internal>")
    return false;

  return true;

#endif
}

int Random::intdev()
{
#ifdef USESPRNG
  return stream->isprng();
#else
  double temp = ran1(&current);
  int val = (int)(INT_MAX*temp);
  return val;
#endif
}
double Random::unidev()
{
#ifdef USESPRNG
  return stream->sprng();
#else
  return ran1(&current);
#endif
}

double Random::gasdev()
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (current < 0) iset=0;
  if  (iset == 0) 
  {
    do 
    {
      v1=2.0*unidev()-1.0;
      v2=2.0*unidev()-1.0;
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

double Random::expdev()
{
  return -log(1.0-unidev());
}

double Random::sindev()
{
  return acos(1.0-2.0*unidev());
}

double Random::sindev(double u)
{
  double scale = 1.0;
  if(u < 3.14159265359 * 0.5)
    scale = sin(u);
  
  while(true)
    {
      double x = u*unidev();
      if( unidev() < sin(x)/scale )
	return x;
    }
}

double Random::randomDistribution1()
{
  while(true)
    {
      double x = expdev()/0.4;

      if( unidev() < ( 0.5*x*x*exp(-x) )/( 0.8*exp(-0.4*x) ) )
	{
	  return x;
	}
    }
}

double Random::pdf2(double a, double zeta, double x)
{
  return fabs((1.0+a*x)*exp(-zeta * x) * sqrt(x));
}

double Random::randomDistribution2(double a, double zeta, double l, double r)
{
#ifdef QMC_DEBUG
  if(l > r){
    cout << "Error: lower bound (" << l << ") is higher than upper bound (" << r << ")" << endl;
  }
#endif

  double scale = 1.0;
  if(fabs(a) < 1e-50)
    {
      double root1  = 0.5 / zeta;
      double vroot1 = pdf2(a, zeta, root1);
      scale = vroot1;

      /*
      printf("1=%20.10f\n",root1);
      printf("1=%20.10f\n",vroot1);
      */
    } else if(fabs(zeta) < 1e-50){
      double root1  = -1.0 / (3.0*a);
      double vroot1 = pdf2(a, zeta, root1);
      scale = vroot1;

      /*
      printf("1=%20.10f\n",root1);
      printf("1=%20.10f\n",vroot1);
      */
    } else {
      double root1  = (3.0*a - 2.0*zeta + sqrt(9.0*a*a - 4.0*a*zeta + 4.0*zeta*zeta))/(4.0*a*zeta);
      double root2  = (3.0*a - 2.0*zeta - sqrt(9.0*a*a - 4.0*a*zeta + 4.0*zeta*zeta))/(4.0*a*zeta);
      double vroot1 = pdf2(a, zeta, root1);
      double vroot2 = pdf2(a, zeta, root2);
      scale = max(vroot1,vroot2);

      /*
      printf("1=%20.10f 2=%20.10f\n",root1,root2);
      printf("1=%20.10f 2=%20.10f\n",vroot1,vroot2);
      */
    }

  double vl    = pdf2(a, zeta, l);
  double vr    = pdf2(a, zeta, r);
  scale = max(scale,vl);
  scale = max(scale,vr);

  /*
  printf("l=%20.10f r=%20.10f\n",l,r);
  printf("l=%20.10f r=%20.10f scale=%20.10f\n",vl,vr,scale);
  */

  while(true)
    {
      double x = (r-l)*unidev() + l;
      if( unidev() < pdf2(a, zeta, x)/ scale )
	{
	  return x;
	}
    }
}

double Random::randomDistribution3(double F)
{
  //I assume F can never be negative
#ifdef QMC_DEBUG
  if(F <= 0.0){
    cout << "Error: why is F = " << F << " negative?" << endl;
  }
#endif
  double scale = 1.0 + F;
  while(true)
    {
      double x = 2.0*3.14159265359*unidev();
      if( unidev() < (1.0 + F*cos(x))/ scale )
	{
	  return x;
	}
    }
}

/* note #undef's at end of file */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Random::ran1(long *idum)
{
  int j;
  long k;
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


