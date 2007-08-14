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

#include "QMCStopwatches.h"

QMCStopwatches::QMCStopwatches()
{
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

void QMCStopwatches::stop()
{
  if(Initialization.isRunning())        Initialization.stop();
  if(Equilibration.isRunning())         Equilibration.stop();
  if(Propagation.isRunning())           Propagation.stop();
  if(Communication_send.isRunning())    Communication_send.stop();
  if(Communication_reduce.isRunning())  Communication_reduce.stop();
  if(Communication_synch.isRunning())   Communication_synch.stop();
  if(Communication_poll.isRunning())    Communication_poll.stop();
  if(Total.isRunning())                 Total.stop();
}

void QMCStopwatches::reset()
{
  Total.reset();
  Initialization.reset();
  Equilibration.reset();
  Propagation.reset();
  Communication_send.reset();
  Communication_reduce.reset();
  Communication_synch.reset();
  Communication_poll.reset();
  Optimization.reset();
}

Stopwatch * QMCStopwatches::getInitializationStopwatch()
{
  return &(this->Initialization);
}

Stopwatch * QMCStopwatches::getEquilibrationStopwatch()
{
  return &(this->Equilibration);
}

Stopwatch * QMCStopwatches::getPropagationStopwatch()
{
  return &(this->Propagation);
}

Stopwatch * QMCStopwatches::getSendCommandStopwatch()
{
  return &(this->Communication_send);
}

Stopwatch * QMCStopwatches::getGatherPropertiesStopwatch()
{
  return &(this->Communication_reduce);
}

Stopwatch * QMCStopwatches::getCommunicationSynchronizationStopwatch()
{
  return &(this->Communication_synch);
}

Stopwatch * QMCStopwatches::getCommandPollingStopwatch()
{
  return &(this->Communication_poll);
}

Stopwatch * QMCStopwatches::getOptimizationStopwatch()
{
  return &(this->Optimization);
}

Stopwatch * QMCStopwatches::getTotalTimeStopwatch()
{
  return &(this->Total);
}

void QMCStopwatches::operator =(const QMCStopwatches & rhs)
{
  Initialization       = rhs.Initialization;
  Equilibration        = rhs.Equilibration;
  Propagation          = rhs.Propagation;
  Communication_send   = rhs.Communication_send;
  Communication_reduce = rhs.Communication_reduce;
  Communication_synch  = rhs.Communication_synch;
  Communication_poll   = rhs.Communication_poll;
  Optimization         = rhs.Optimization;
  Total                = rhs.Total;
}

QMCStopwatches QMCStopwatches::operator+(QMCStopwatches & rhs)
{
  QMCStopwatches result;

  result.Initialization       = this->Initialization + rhs.Initialization;
  result.Equilibration        = this->Equilibration + rhs.Equilibration;
  result.Propagation          = this->Propagation + rhs.Propagation;
  result.Communication_send   = this->Communication_send + 
                                  rhs.Communication_send;
  result.Communication_reduce = this->Communication_reduce + 
                                  rhs.Communication_reduce;
  result.Communication_synch  = this->Communication_synch + 
                                  rhs.Communication_synch;
  result.Communication_poll   = this->Communication_poll + 
                                  rhs.Communication_poll;
  result.Optimization         = this->Optimization + rhs.Optimization;
  result.Total                = this->Total + rhs.Total;

  return result;
}

ostream& operator <<(ostream& strm, QMCStopwatches & rhs)
{
  double total = rhs.getTotalTimeStopwatch()->timeUS() / 100.0;
  strm.setf(ios::fixed);
  strm.unsetf(ios::scientific);
  int prec = 2;
  int width = 10;
  strm << "Total Time:             " << *rhs.getTotalTimeStopwatch() 
       << setprecision(prec) << setw(width) << (double) rhs.getTotalTimeStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Initialization Time:    " << *rhs.getInitializationStopwatch() 
       << setprecision(prec) << setw(width) << (double) rhs.getInitializationStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Equilibration Time:     " << *rhs.getEquilibrationStopwatch()
       << setprecision(prec) << setw(width) << (double) rhs.getEquilibrationStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Propagation Time:       " << *rhs.getPropagationStopwatch() 
       << setprecision(prec) << setw(width) << (double) rhs.getPropagationStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Send Command Time:      " << *rhs.getSendCommandStopwatch()
       << setprecision(prec) << setw(width) << (double) rhs.getSendCommandStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Gather Properties Time: " << *rhs.getGatherPropertiesStopwatch() 
       << setprecision(prec) << setw(width) << (double) rhs.getGatherPropertiesStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Synchronization Time:   " << *rhs.getCommunicationSynchronizationStopwatch()
       << setprecision(prec) << setw(width) << (double) rhs.getCommunicationSynchronizationStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Poll for Command Time:  " << *rhs.getCommandPollingStopwatch() 
       << setprecision(prec) << setw(width) << (double) rhs.getCommandPollingStopwatch()->timeUS() / total
       << " %" << endl;
  strm << "Optimization Time:      " << *rhs.getOptimizationStopwatch() 
       << setprecision(prec) << setw(width) << (double) rhs.getOptimizationStopwatch()->timeUS() / total
       << " %" << endl;

  return strm;
}



#ifdef PARALLEL

bool QMCStopwatches::mpiTypeCreated = false;

void QMCStopwatches::buildMpiReduce()
{
  MPI_Op_create((MPI_User_function*)Reduce_Function,
                true,&MPI_REDUCE);
}

void QMCStopwatches::buildMpiType()
{
  QMCStopwatches indata;

  const int numberDataTypes = 9;

  int          block_lengths[numberDataTypes];
  MPI_Aint     displacements[numberDataTypes];
  MPI_Aint     addresses[numberDataTypes+1];
  MPI_Datatype typelist[numberDataTypes];

  typelist[0] = Stopwatch::MPI_TYPE; 
  typelist[1] = Stopwatch::MPI_TYPE; 
  typelist[2] = Stopwatch::MPI_TYPE; 
  typelist[3] = Stopwatch::MPI_TYPE; 
  typelist[4] = Stopwatch::MPI_TYPE; 
  typelist[5] = Stopwatch::MPI_TYPE; 
  typelist[6] = Stopwatch::MPI_TYPE; 
  typelist[7] = Stopwatch::MPI_TYPE; 
  typelist[8] = Stopwatch::MPI_TYPE;
 
  block_lengths[0] = 1;
  block_lengths[1] = 1;
  block_lengths[2] = 1;
  block_lengths[3] = 1;
  block_lengths[4] = 1;
  block_lengths[5] = 1;
  block_lengths[6] = 1;
  block_lengths[7] = 1;
  block_lengths[8] = 1;

  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.Initialization), &addresses[1]);
  MPI_Address(&(indata.Equilibration), &addresses[2]);
  MPI_Address(&(indata.Propagation), &addresses[3]);
  MPI_Address(&(indata.Communication_send), &addresses[4]);
  MPI_Address(&(indata.Communication_reduce), &addresses[5]);
  MPI_Address(&(indata.Communication_synch), &addresses[6]);
  MPI_Address(&(indata.Communication_poll), &addresses[7]);
  MPI_Address(&(indata.Optimization), &addresses[8]);
  MPI_Address(&(indata.Total), &addresses[9]);

  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];
  displacements[3] = addresses[4] - addresses[0];
  displacements[4] = addresses[5] - addresses[0];
  displacements[5] = addresses[6] - addresses[0];
  displacements[6] = addresses[7] - addresses[0];
  displacements[7] = addresses[8] - addresses[0];
  displacements[8] = addresses[9] - addresses[0];
  
  MPI_Type_struct(numberDataTypes, block_lengths, displacements, typelist, 
                  &MPI_TYPE);
  MPI_Type_commit(&MPI_TYPE);
}


MPI_Datatype QMCStopwatches::MPI_TYPE;

MPI_Op QMCStopwatches::MPI_REDUCE;

void QMCStopwatches::Reduce_Function(QMCStopwatches *in, 
				     QMCStopwatches *inout, 
				     int *len, MPI_Datatype *dptr)
{
  *inout = *inout + *in;
}


#endif

