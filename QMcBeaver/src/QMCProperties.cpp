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

#include "QMCProperties.h"

QMCProperties::QMCProperties()
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

void QMCProperties::zeroOut()
{
  energy.zeroOut();
  kineticEnergy.zeroOut();
  potentialEnergy.zeroOut();
  logWeights.zeroOut();
  acceptanceProbability.zeroOut();
  distanceMovedAccepted.zeroOut();
  distanceMovedTrial.zeroOut();
}

void QMCProperties::operator = ( const QMCProperties &rhs )
{
  energy                = rhs.energy;
  kineticEnergy         = rhs.kineticEnergy;
  potentialEnergy       = rhs.potentialEnergy;
  logWeights            = rhs.logWeights;
  acceptanceProbability = rhs.acceptanceProbability;
  distanceMovedAccepted = rhs.distanceMovedAccepted;
  distanceMovedTrial    = rhs.distanceMovedTrial;
}

QMCProperties QMCProperties::operator + ( QMCProperties &rhs )
{
  QMCProperties result;

  result.energy = energy+rhs.energy;
  result.kineticEnergy = kineticEnergy + rhs.kineticEnergy;
  result.potentialEnergy = potentialEnergy + rhs.potentialEnergy;
  result.logWeights = logWeights+rhs.logWeights;
  result.acceptanceProbability = acceptanceProbability + 
    rhs.acceptanceProbability;
  result.distanceMovedAccepted = distanceMovedAccepted + 
    rhs.distanceMovedAccepted;
  result.distanceMovedTrial = distanceMovedTrial + 
    rhs.distanceMovedTrial;

  return result;
}

void QMCProperties::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCProperties>" << endl;

  // Energy
  strm << "<Energy>" << endl;
  energy.toXML(strm);
  strm << "</Energy>" << endl;

  // Kinetic energy
  strm << "<KineticEnergy>" << endl;
  kineticEnergy.toXML(strm);
  strm << "</KineticEnergy>" << endl;

  // Potential energy
  strm << "<PotentialEnergy>" << endl;
  potentialEnergy.toXML(strm);
  strm << "</PotentialEnergy>" << endl;

  // log weights
  strm << "<LogWeights>" << endl;
  logWeights.toXML(strm);
  strm << "</LogWeights>" << endl;

  // Acceptance probability
  strm << "<AcceptanceProbability>" << endl;
  acceptanceProbability.toXML(strm);
  strm << "</AcceptanceProbability>" << endl;

  // Distance Moved Accepted
  strm << "<DistanceMovedAccepted>" << endl;
  distanceMovedAccepted.toXML(strm);
  strm << "</DistanceMovedAccepted>" << endl;

  // Distance Moved Trial
  strm << "<DistanceMovedTrial>" << endl;
  distanceMovedTrial.toXML(strm);
  strm << "</DistanceMovedTrial>" << endl;

  // Close XML
  strm << "</QMCProperties>" << endl;
}

void QMCProperties::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;

  // Read energy
  strm >> temp;
  energy.readXML(strm);
  strm >> temp;

  // Read kinetic energy
  strm >> temp;
  kineticEnergy.readXML(strm);
  strm >> temp;

  // Read potential energy
  strm >> temp;
  potentialEnergy.readXML(strm);
  strm >> temp;

  // Read log weights
  strm >> temp;
  logWeights.readXML(strm);
  strm >> temp;

  // Acceptance probability
  strm >> temp;
  acceptanceProbability.readXML(strm);
  strm >> temp;

  // Distance Moved Accepted
  strm >> temp;
  distanceMovedAccepted.readXML(strm);
  strm >> temp;

  // Distance Moved Trial
  strm >> temp;
  distanceMovedTrial.readXML(strm);
  strm >> temp;

  // Close XML
  strm >> temp;
}

ostream& operator <<(ostream& strm, QMCProperties &rhs)
{
  strm << "----------------- Energy --------------------" << endl;
  strm << rhs.energy;

  strm << "-------------- Kinetic Energy ---------------" << endl;
  strm << rhs.kineticEnergy;

  strm << "------------ Potential Energy ---------------" << endl;
  strm << rhs.potentialEnergy;

  strm << "----------- AcceptanceProbability -----------" << endl;
  strm << rhs.acceptanceProbability;


  strm << "----------- DistanceMovedAccepted -----------" << endl;
  strm << rhs.distanceMovedAccepted;


  strm << "------------ DistanceMovedTrial -------------" << endl;
  strm << rhs.distanceMovedTrial;


  strm << "---------------- logWeights -----------------" << endl;
  strm << rhs.logWeights;

  return strm;
}





#ifdef PARALLEL

bool QMCProperties::mpiTypeCreated = false;

void QMCProperties::buildMpiReduce()
{
  MPI_Op_create((MPI_User_function*)Reduce_Function,
                true,&MPI_REDUCE);
}

void QMCProperties::buildMpiType()
{
  QMCProperties indata;

  // The number of properties 
  // ADJUST THIS WHEN ADDING NEW PROPERTIES
  const int NumberOfProperties = 7;

  int          block_lengths[NumberOfProperties];
  MPI_Aint     displacements[NumberOfProperties];
  MPI_Aint     addresses[NumberOfProperties+1];
  MPI_Datatype typelist[NumberOfProperties];

  // Set up the types and lengths of each part of the new MPI data struct
  for(int i=0; i<NumberOfProperties; i++)
    {
      typelist[i] = QMCProperty::MPI_TYPE;
      block_lengths[i] = 1;
    }
  
  // Find the addresses of all of the data elements in the struct
  // ADJUST THIS WHEN ADDING NEW PROPERTIES
  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.energy), &addresses[1]);
  MPI_Address(&(indata.logWeights), &addresses[2]);
  MPI_Address(&(indata.acceptanceProbability), &addresses[3]);
  MPI_Address(&(indata.distanceMovedAccepted), &addresses[4]);
  MPI_Address(&(indata.distanceMovedTrial), &addresses[5]);  
  MPI_Address(&(indata.kineticEnergy), &addresses[6]);  
  MPI_Address(&(indata.potentialEnergy), &addresses[7]);  


  // Find the relative addresses of the data elements to the start of 
  // the struct
  for(int i=0; i<NumberOfProperties; i++)
    {  
      displacements[i] = addresses[i+1] - addresses[0];
    }
  
  MPI_Type_struct(NumberOfProperties, block_lengths, displacements, typelist, 
                  &MPI_TYPE);
  MPI_Type_commit(&MPI_TYPE);
}

MPI_Datatype QMCProperties::MPI_TYPE;

MPI_Op QMCProperties::MPI_REDUCE;

void QMCProperties::Reduce_Function(QMCProperties *in, QMCProperties *inout, 
                                   int *len, MPI_Datatype *dptr)
{
  *inout = *inout + *in;
}


#endif

