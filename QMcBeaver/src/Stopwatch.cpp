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

#include "Stopwatch.h"

Stopwatch::Stopwatch()
{
  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;

  reset();


#ifdef PARALLEL
  if( !mpiTypeCreated )
    {
      mpiTypeCreated = true;
      buildMpiType();
      buildMpiReduce();
    }
#endif
}

void Stopwatch::reset()
{
  total_us=0;
  running=false;
}

void Stopwatch::start()
{
  if(running)
    {
      //stopped a stopwatch that was not running
      cerr << "ERROR: You just started a stopwatch that was running" << endl;
      cerr << "check to make sure this is not a bug" << endl;
      cerr << "exit in void Stopwatch::start()" << endl;
      exit(1);
    }
  else
    {
      gettimeofday(&tp,&tz);
      stime1 = tp.tv_sec;
      micro1 = tp.tv_usec;
      running=true;
    }
}

void Stopwatch::stop()
{
  if(running)
    {
      gettimeofday(&tp,&tz);
      stime2 = tp.tv_sec;
      micro2 = tp.tv_usec;
      result_us = ((stime2-stime1)*1e6 + micro2 - micro1);
      total_us += result_us;
      running   = false;
    }
  else
    {
      //stopped a stopwatch that was not running
      cerr << "ERROR: You just stopped a stopwatch that was not running" 
	   << endl;
      cerr << "check to make sure this is not a bug" << endl;
      cerr << "exit in void Stopwatch::stop()" << endl;
      exit(1);
    }
}

long Stopwatch::timeMS()
{
  return total_us/1000;
}

long Stopwatch::timeUS()
{
  return total_us;
}

bool Stopwatch::isRunning()
{
  return running;
}

string Stopwatch::toString()
{
  ostringstream stream;
  stream << timeMS() << " ms " << ends;
  return stream.str();
}

void Stopwatch::operator =( const Stopwatch & rhs )
{
  total_us = rhs.total_us;

  // don't include any members that are not in the mpi type because
  // they can cause seg faults with some MPI implementations.
}  

Stopwatch Stopwatch::operator+(Stopwatch & rhs)
{
  Stopwatch result;
  result.total_us = this->total_us + rhs.total_us;
  result.running = false;

  return result;
}


ostream& operator <<(ostream& strm, Stopwatch &watch)
{
  strm << watch.toString();
  return strm;
}



#ifdef PARALLEL

bool Stopwatch::mpiTypeCreated = false;

void Stopwatch::buildMpiReduce()
{
  MPI_Op_create((MPI_User_function*)Reduce_Function,
                true,&MPI_REDUCE);
}

void Stopwatch::buildMpiType()
{
  Stopwatch indata;

  int          block_lengths[1];
  MPI_Aint     displacements[1];
  MPI_Aint     addresses[2];
  MPI_Datatype typelist[1];

  typelist[0] = MPI_LONG_INT;    // total time
  
  block_lengths[0] = 1;
  
  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.total_us), &addresses[1]);

  displacements[0] = addresses[1] - addresses[0];
  
  MPI_Type_struct(1, block_lengths, displacements, typelist, 
                  &MPI_TYPE);
  MPI_Type_commit(&MPI_TYPE);
}


MPI_Datatype Stopwatch::MPI_TYPE;

MPI_Op Stopwatch::MPI_REDUCE;

void Stopwatch::Reduce_Function(Stopwatch *in, Stopwatch *inout, 
				int *len, MPI_Datatype *dptr)
{
  *inout = *inout + *in;
}


#endif









