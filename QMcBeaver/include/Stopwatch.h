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


#ifndef Stopwatch_H
#define Stopwatch_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

/*
  When using longs for microseconds (useful in GPU timings) then
  the longest the datatype will work for is (assuming 32 bits) ~ 1.2 hrs
*/
typedef long long longType;
//typedef long longType;

typedef double microType;

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
struct timeval {
        long         tv_sec;
        microType    tv_usec;
};
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
  longType stime1, stime2, result_us, total_us;
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

  longType timeMS();

  /**
    Gets the time in microseconds.
    */

  longType timeUS();

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
