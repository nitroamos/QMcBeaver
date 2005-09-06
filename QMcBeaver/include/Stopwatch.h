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

#ifndef Stopwatch_H
#define Stopwatch_H

#ifdef _WIN32

#include <sys/timeb.h>
#include <sys/types.h>
//#include <winsock.h>
#include <windows.h>
#define USE_HIGH_PRECISION

struct timezone {
    int tz_minuteswest;
    int tz_dsttime;
};

//timeval is defined in winsock, but cygwin's version seems to have a problem...
//included here just in case winsock doesn't play nice.
#if defined(__CYGWIN__)
struct timeval {
        long    tv_sec;
        long    tv_usec;
};
#endif

/**
  High precision timing is not available to Windows through time.h, and I don't
  know how precise it is on Linux, Unix, etc. This function has been "overloaded"
  to use special Windows API when it is running on Windows/cygwin. This precision
  is probably only important for GPU related timings. When I get to Linux,
  I'll figure out how to make a function that will work there as well.
  
  For Intel Pentiums, there is a command called rdtsc. It's possible this works on
  AMD cpus also. From the Intel software forums, this was recommended as sample code:
  
  unsigned __int64 tick = 0;
  void getTick(void)
  { 
     __asm{
        _emit 0x0f
        _emit 0x31
        lea ecx,tick
        mov dword ptr [ecx],  eax
        mov dword ptr [ecx+4],edx
      }
  }

  main() {
    while (1) {
      getTick();
      printf("%I64u\n", tick);
    } 
  }

  I believe this code would have to be matched with a frequency number to turn the number
  of ticks into an actual time.
*/
static void gettimeofday(struct timeval* t,void* timezone){
#ifdef USE_HIGH_PRECISION
  LARGE_INTEGER timeLI, freqLI;
  QueryPerformanceCounter(&timeLI);
  QueryPerformanceFrequency(&freqLI);
  t->tv_sec=0;
  double temp = 1.0e6*(double)(timeLI.QuadPart)/freqLI.QuadPart;
  t->tv_usec=(long)(temp + 0.5);
#else
  struct _timeb timebuffer;
  _ftime( &timebuffer );
  t->tv_sec=(long)timebuffer.time;
  t->tv_usec=1000*timebuffer.millitm;
#endif
}
#else
#include <sys/time.h>
#endif


#include <string>
#include <sstream>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

using namespace std;

/**
  An accurate software stopwatch.
  */

class Stopwatch
{
 private:
  long stime1, stime2, result_ms, total_ms;
  int milli1, milli2;
  bool running;
  struct timeval tp;
  struct timezone tz;

public:
  /**
    Creates an instance of the stopwatch that is zeroed and not running.
    */

  Stopwatch();


  /**
    Resets and stops the stopwatch.
    */

  void reset();


  /**
    Starts the stopwatch.
    */

  void start();


  /**
    Stops the stopwatch.
    */

  void stop();


  /**
    Gets the time in milliseconds.
    */

  long timeMS();


  /**
    Returns true if the stopwatch is running and false otherwise.
    */

  bool isRunning();


  /**
    Gets the time formatted as a string.
    */

  string toString();
  
  /**
     Sets two objects equal.
    */

  void operator = ( const Stopwatch &rhs);

  /**
    Returns a stopwatch which contains the total time from two
    stopwatch objects.
    */
  Stopwatch operator+(Stopwatch & rhs);

  /**
    Formats and prints the time to a stream.
    */
  friend ostream& operator <<(ostream& strm, Stopwatch &watch);


#ifdef PARALLEL

private:
  /**
    A flag which tells if MPI_TYPE has been generated.
    */

  static bool mpiTypeCreated;

  
  /**
    Build MPI_TYPE.
    */

  static void buildMpiType();

  /**
    Build MPI_REDUCE.
    */

  static void buildMpiReduce();


  /**
    An MPI function which allows MPI_Reduce to be used in adding Stopwatches.
    */

  static void Reduce_Function(Stopwatch *in, Stopwatch *inout, 
                              int *len, MPI_Datatype *dptr);

public:

  /**
    The MPI data type for a Stopwatch.
    */
  static MPI_Datatype MPI_TYPE;


  /** 
    The MPI operation for performing MPI_Reduce on Stopwatch objects.
    */
  
  static MPI_Op MPI_REDUCE;

#endif

};

#endif
