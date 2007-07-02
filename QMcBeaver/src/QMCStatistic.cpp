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
  //sum      = 0.0;
  //sum2     = 0.0;
  //weights  = 0.0;
  nsamples = 0;
}

unsigned long QMCStatistic::getNumberSamples() const
{
  return nsamples;
}

double QMCStatistic::getAverage() const
{
  if( nsamples == 0) return 0;
  return (sum/weights);
}

double QMCStatistic::getVariance() const
{
  if( nsamples == 0) return 0;
  return (sum2/weights-getAverage()*getAverage());
}

double QMCStatistic::getStandardDeviation() const
{
  return sqrt(getVariance());
}

void QMCStatistic::newSample(long double s, long double weight)
{
  if(nsamples == 0)
    {
      weights = weight;
      sum     = s*weight;
      sum2    = s*s*weight;
    } else {
      weights += weight;
      sum     += s*weight;
      sum2    += s*s*weight;
    }
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
  if(nsamples == 0)
    {
      result = rhs;
    } else {
      result.weights  = weights+rhs.weights;
      result.sum      = sum+rhs.sum;
      result.sum2     = sum2+rhs.sum2;
      result.nsamples = nsamples+rhs.nsamples;
    }
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

  if(nsamples == 0)
    {
      // sum
      strm << "<Sum>\n" << 0.0 << "\n</Sum>" << endl;
      // sum sq
      strm << "<SumSquared>\n" << 0.0 << "\n</SumSquared>" << endl; 
      // weights
      strm << "<SumWeights>\n" << 0.0 << "\n</SumWeights>" << endl;
    } else {
      // sum
      strm << "<Sum>\n" << double(sum) << "\n</Sum>" << endl;
      // sum sq
      strm << "<SumSquared>\n" << double(sum2) << "\n</SumSquared>" << endl; 
      // weights
      strm << "<SumWeights>\n" << double(weights) << "\n</SumWeights>" << endl;
    }
  
  // nsamples
  strm << "<NumberOfSamples>\n" << nsamples << "\n</NumberOfSamples>" << endl;

  // Close XML
  strm << "</QMCStatistic>" << endl;
}

bool QMCStatistic::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;
  if (temp != "<QMCStatistic>")
    return false;

  // sum
  strm >> temp;
  if (temp != "<Sum>")
    return false;
  strm >> temp;
  sum = atof(temp.c_str());
  strm >> temp;
  if (temp != "</Sum>")
    return false;
  
  // sum sq
  strm >> temp;
  if (temp != "<SumSquared>")
    return false;
  strm >> temp;
  sum2 = atof(temp.c_str());
  strm >> temp;
  if (temp != "</SumSquared>")
    return false;

  // weights
  strm >> temp;
  if (temp != "<SumWeights>")
    return false;
  strm >> temp;
  weights = atof(temp.c_str());
  strm >> temp;
  if (temp != "</SumWeights>")
    return false;

  // nsamples
  strm >> temp;
  if (temp != "<NumberOfSamples>")
    return false;
  strm >> temp;
  nsamples = atoi(temp.c_str());
  strm >> temp;
  if (temp != "</NumberOfSamples>")
    return false;

  // Close XML
  strm >> temp;
  if (temp != "</QMCStatistic>")
    return false;

  return true;
}

ostream& operator << (ostream& strm, const QMCStatistic &rhs)
{
  strm.precision(12);
  strm.width(20);
  strm << scientific << rhs.getAverage() << " +/- ";
  if( fabs(rhs.getStandardDeviation()) > 1e-300 )
    if( log( fabs( rhs.getStandardDeviation() )) > 10.0)
      strm << scientific;
  strm.precision(12);
  strm.width(20);
  strm << scientific << rhs.getStandardDeviation() << " (" << rhs.getNumberSamples() << " samples)" << endl;
  
  //strm << rhs.getAverage() << " +/- " << rhs.getStandardDeviation();
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








