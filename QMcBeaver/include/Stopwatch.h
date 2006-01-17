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

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

#if defined(_WIN32) || defined(__CYGWIN__)

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
typedef double microType;
struct timeval {
        long         tv_sec;
        microType    tv_usec;
};
#else
typedef double microType;
#endif //defined(__CYGWIN__)

static void gettimeofday(struct timeval* t,void* timezone){
#ifdef USE_HIGH_PRECISION
  LARGE_INTEGER timeLI, freqLI;
  QueryPerformanceCounter(&timeLI);
  QueryPerformanceFrequency(&freqLI);
  t->tv_sec=0;
  double temp = 1.0e6*(double)(timeLI.QuadPart)/freqLI.QuadPart;
  t->tv_usec = temp;
#else
  struct _timeb timebuffer;
  _ftime( &timebuffer );
  t->tv_sec=(long)timebuffer.time;
  t->tv_usec=1000*timebuffer.millitm;
#endif
}

#else
#include <sys/time.h>
#endif//defined(_WIN32) || defined(__CYGWIN__)


#include <string>
#include <sstream>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

/**
  An accurate software stopwatch.
  */

class Stopwatch
{
 private:
  long stime1, stime2, result_us, total_us;
  microType micro1, micro2;
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
    Gets the time in microseconds.
    */

  long timeUS();

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
