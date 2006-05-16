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

#include "QMCStatistic.h"

QMCStatistic::QMCStatistic()
{
  zeroOut();

#ifdef PARALLEL
  if( !mpiTypeCreated )
    {
      mpiTypeCreated = true;
      buildMpiType();
      buildMpiReduce();
    }
#endif
}
  

void QMCStatistic::zeroOut()
{
  sum      = 0.0;
  sum2     = 0.0;
  weights  = 0.0;
  nsamples = 0;
}

unsigned long QMCStatistic::getNumberSamples()
{
  return nsamples;
}

double QMCStatistic::getAverage()
{
  if( nsamples == 0) return 0;
  return (sum/weights);
}

double QMCStatistic::getVariance()
{
  if( nsamples == 0) return 0;
  return (sum2/weights-getAverage()*getAverage());
}

double QMCStatistic::getStandardDeviation()
{
  return sqrt(getVariance());
}

void QMCStatistic::newSample(double s, double weight)
{
  weights += weight;
  sum     += s*weight;
  sum2    += s*s*weight;
  nsamples++;
}

void QMCStatistic::operator = ( const QMCStatistic & rhs )
{
  weights = rhs.weights;
  sum = rhs.sum;
  sum2 = rhs.sum2;
  nsamples = rhs.nsamples;
}

QMCStatistic QMCStatistic::operator + (const QMCStatistic &rhs)
{
  QMCStatistic result;
  result.weights  = weights+rhs.weights;
  result.sum      = sum+rhs.sum;
  result.sum2     = sum2+rhs.sum2;
  result.nsamples = nsamples+rhs.nsamples;
  return result;
}

void QMCStatistic::reWeight(double w)
{
  //weights *= w;
  sum     *= w;
  sum2    *= w;
}

void QMCStatistic::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCStatistic>" << endl;
  
  // sum
  strm << "<Sum>\n" << double(sum) << "\n</Sum>" << endl;
  
  // sum sq
  strm << "<SumSquared>\n" << double(sum2) << "\n</SumSquared>" << endl;

  // weights
  strm << "<SumWeights>\n" << double(weights) << "\n</SumWeights>" << endl;

  // nsamples
  strm << "<NumberOfSamples>\n" << nsamples << "\n</NumberOfSamples>" << endl;

  // Close XML
  strm << "</QMCStatistic>" << endl;
}

void QMCStatistic::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;
  
  // sum
  strm >> temp;
  strm >> temp;
  sum = atof(temp.c_str());
  strm >> temp;
  
  // sum sq
  strm >> temp;
  strm >> temp;
  sum2 = atof(temp.c_str());
  strm >> temp;

  // weights
  strm >> temp;
  strm >> temp;
  weights = atof(temp.c_str());
  strm >> temp;

  // nsamples
  strm >> temp;
  strm >> temp;
  nsamples = atoi(temp.c_str());
  strm >> temp;

  // Close XML
  strm >> temp;
}

ostream& operator <<(ostream& strm, QMCStatistic &rhs)
{
  strm << rhs.getAverage() << " +/- " << rhs.getStandardDeviation();
  return strm;
}


#ifdef PARALLEL

bool QMCStatistic::mpiTypeCreated = false;

void QMCStatistic::buildMpiReduce()
{
  MPI_Op_create((MPI_User_function*)Reduce_Function,
                true,&MPI_REDUCE);
}

void QMCStatistic::buildMpiType()
{
  QMCStatistic indata;

  int          block_lengths[4];
  MPI_Aint     displacements[4];
  MPI_Aint     addresses[5];
  MPI_Datatype typelist[4];

  typelist[0] = MPI_LONG_DOUBLE; //this is the weights 
                                 //(roughly the # of samples)
  typelist[1] = MPI_LONG_DOUBLE; //this is the sum
  typelist[2] = MPI_LONG_DOUBLE; //this is the sum2
  typelist[3] = MPI_LONG;        //this is nsamples
  
  block_lengths[0] = 1;
  block_lengths[1] = 1;
  block_lengths[2] = 1;
  block_lengths[3] = 1;
  
  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.weights), &addresses[1]);
  MPI_Address(&(indata.sum), &addresses[2]);
  MPI_Address(&(indata.sum2), &addresses[3]);
  MPI_Address(&(indata.nsamples), &addresses[4]);

  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];
  displacements[3] = addresses[4] - addresses[0];
  
  MPI_Type_struct(4, block_lengths, displacements, typelist, 
                  &MPI_TYPE);
  MPI_Type_commit(&MPI_TYPE);
}


MPI_Datatype QMCStatistic::MPI_TYPE;

MPI_Op QMCStatistic::MPI_REDUCE;

void QMCStatistic::Reduce_Function(QMCStatistic *in, QMCStatistic *inout, 
				   int *len, MPI_Datatype *dptr)
{
  *inout = *inout + *in;
}


#endif








