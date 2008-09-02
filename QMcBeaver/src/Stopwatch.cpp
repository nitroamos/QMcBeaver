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


#include "Stopwatch.h"

Stopwatch::Stopwatch()
{
  tz.tz_minuteswest = 0;
  tz.tz_dsttime     = 0;
  name              = "";
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

Stopwatch::~Stopwatch()
{
}

void Stopwatch::reset()
{
  stime1   = 0;
  micro1   = 0;
  total_us = 0;
  running  = false;

  average_us = 0;
  numSamples = -5;
}

void Stopwatch::reset(string newName)
{
  reset();
  setName(newName);
}

void Stopwatch::setName(string newName)
{
  name = newName;
}

string Stopwatch::getName()
{
  return name;
}

void Stopwatch::start()
{
  if(running)
    {
#ifdef QMC_DEBUG
      //stopped a stopwatch that was not running
      //cerr << "Warning: You just started a stopwatch that was already running." << endl;
#endif
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
      result_us = (longType)((stime2-stime1)*1e6 + micro2 - micro1);
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

void Stopwatch::lap(int weight)
{
  stop();
  if(numSamples >= 0)
    average_us += total_us;
  numSamples += weight;
  total_us = 0;
}

void Stopwatch::lap()
{
  lap(1);
}

double Stopwatch::getAverage()
{
  if(numSamples <= 0)
    return 0.0;
  
  double avg = (double)average_us;
  return avg/numSamples;
}

void Stopwatch::print(ostream & strm)
{
  strm << setw(30) << name << ": " << (int)(total_us+0.5);
  strm << " ( avg = " << setprecision(2) << setw(10) << getAverage() << " us";
  strm << ", n=" << setw(15) << numSamples << ")";
}

longType Stopwatch::timeMS()
{
  return total_us/1000;
}

longType Stopwatch::timeUS()
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
  longType time = timeMS();
  if(time <= 0 && isRunning())
    {
      gettimeofday(&tp,&tz);
      stime2 = tp.tv_sec;
      micro2 = tp.tv_usec;
      time = (longType)((stime2-stime1)*1e6 + micro2 - micro1);
      time /= 1000;
    }
  double hrs = time/(1000.0*60.0*60.0);
  double mins = time/(1000.0*60.0);
  double secs = time/(1000.0);
  int wdth = 20;
  int prec  = 5;
  stream.width(wdth); stream.precision(prec);
  stream << time << " ms  " << " = ";
  stream.width(wdth); stream.precision(prec);
  stream << secs << " s   " << " = ";
  stream.width(wdth); stream.precision(prec);
  stream << mins << " mins" << " = ";
  stream.width(wdth); stream.precision(prec);
  stream << hrs <<  " hrs ";
  return stream.str();
}

void Stopwatch::operator =( const Stopwatch & rhs )
{
  total_us   = rhs.total_us;
  average_us = rhs.average_us;
  numSamples = rhs.numSamples;

  // don't include any members that are not in the mpi type because
  // they can cause seg faults with some MPI implementations.
}  

Stopwatch Stopwatch::operator+(Stopwatch & rhs)
{
  Stopwatch result;
  result.total_us   = this->total_us + rhs.total_us;
  result.average_us = this->average_us + rhs.average_us;
  result.numSamples = this->numSamples + rhs.numSamples;
  result.running = false;

  return result;
}

void Stopwatch::aggregateTimer(Stopwatch & rhs)
{
  if(name != "" && name != rhs.name){
    cerr << "Warning: name mismatch in Stopwatch. " << name << " != " << rhs.name << endl;
  }
  setName(rhs.name);
  *this = *this + rhs;
  rhs.reset();
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

  int          block_lengths[3];
  MPI_Aint     displacements[3];
  MPI_Aint     addresses[4];
  MPI_Datatype typelist[3];

  //typelist[0] = MPI_LONG_INT;    // total time
  typelist[0] = MPI_LONG_LONG_INT;    // total time
  typelist[1] = MPI_LONG_LONG_INT;    // total time
  typelist[2] = MPI_LONG_LONG_INT;    // total time
  
  block_lengths[0] = 1;
  block_lengths[1] = 1;
  block_lengths[2] = 1;
  
  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.total_us),   &addresses[1]);
  MPI_Address(&(indata.average_us), &addresses[2]);
  MPI_Address(&(indata.numSamples), &addresses[3]);

  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];

  //See discussion in QMCProperties for what this is all about
#ifdef QMC_OLDMPICH
  MPI_Type_struct(3, block_lengths, displacements, typelist, 
                  &MPI_TYPE);
#else
  MPI_Datatype temp;
  MPI_Type_struct(3, block_lengths, displacements, typelist, 
                  &temp);
  MPI_Type_create_resized(temp,0,sizeof(Stopwatch),&MPI_TYPE);
#endif

  MPI_Type_commit(&MPI_TYPE);
}

MPI_Datatype Stopwatch::MPI_TYPE;

MPI_Op Stopwatch::MPI_REDUCE;

void Stopwatch::Reduce_Function(Stopwatch *in, Stopwatch *inout, 
				int *len, MPI_Datatype *dptr)
{
  for(int i=0; i < *len; i++)
    inout[i] = inout[i] + in[i];
}


#endif
