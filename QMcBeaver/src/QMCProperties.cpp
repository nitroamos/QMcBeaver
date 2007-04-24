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
  calc_density = false;
  nBasisFunc   = 0;
  zeroOut();

#ifdef PARALLEL
  if (!mpiTypeCreated)
    {
      mpiTypeCreated = true;
      buildMpiType();
      buildMpiReduce();
    }
#endif
}

QMCProperties::~QMCProperties()
{
  if(calc_density)
    {
      chiDensity.deallocate();
    }
}

void QMCProperties::setCalcDensity(bool calcDensity, int nbasisfunctions)
{
  calc_density = calcDensity;

  if(calc_density == true)
    {
      nBasisFunc   = nbasisfunctions;
      chiDensity.allocate(nBasisFunc);

      for (int i=0; i<nBasisFunc; i++)
	chiDensity(i).zeroOut();
    }
  else
    {
      nBasisFunc = 0;
      chiDensity.deallocate();
    }
}

void QMCProperties::zeroOut()
{
  walkerAge.zeroOut();
  weightChange.zeroOut();
  growthRate.zeroOut();

  energy.zeroOut();
  kineticEnergy.zeroOut();
  potentialEnergy.zeroOut();
  neEnergy.zeroOut();
  eeEnergy.zeroOut();
  logWeights.zeroOut();
  acceptanceProbability.zeroOut();
  distanceMovedAccepted.zeroOut();
  distanceMovedTrial.zeroOut();
  
  if (calc_density)
    for (int i=0; i<nBasisFunc; i++)
      chiDensity(i).zeroOut();
}

void QMCProperties::newSample(QMCProperties* newProperties, double weight, 
			      int nwalkers)
{
  if(newProperties->walkerAge.getNumberSamples() > 0)
    walkerAge.newSample(newProperties->walkerAge.getAverage(), 1.0);
  if(newProperties->weightChange.getNumberSamples() > 0)
    weightChange.newSample(newProperties->weightChange.getAverage(), 1.0);;
  if(newProperties->growthRate.getNumberSamples() > 0)
    growthRate.newSample(newProperties->growthRate.getAverage(), 1.0);

  energy.newSample(newProperties->energy.getAverage(), weight);
  kineticEnergy.newSample(newProperties->kineticEnergy.getAverage(), weight);
  potentialEnergy.newSample(newProperties->potentialEnergy.getAverage(),
			    weight);
  neEnergy.newSample(newProperties->neEnergy.getAverage(), weight);
  eeEnergy.newSample(newProperties->eeEnergy.getAverage(), weight);
  acceptanceProbability.newSample
    (newProperties->acceptanceProbability.getAverage(), weight);
  distanceMovedAccepted.newSample
    (newProperties->distanceMovedAccepted.getAverage(), weight);
  distanceMovedTrial.newSample(newProperties->distanceMovedTrial.getAverage(),
			       weight);
  logWeights.newSample(newProperties->logWeights.getAverage(), nwalkers);
  
  if (calc_density)
    for (int i=0; i<nBasisFunc; i++)
      chiDensity(i).newSample(newProperties->chiDensity(i).getAverage(),weight);
}

void QMCProperties::matchParametersTo( const QMCProperties &rhs )
{
  setCalcDensity(rhs.calc_density,rhs.nBasisFunc);
}

void QMCProperties::operator = ( const QMCProperties &rhs )
{
  matchParametersTo(rhs);

  walkerAge             = rhs.walkerAge;
  weightChange          = rhs.weightChange;
  growthRate            = rhs.growthRate;

  energy                = rhs.energy;
  kineticEnergy         = rhs.kineticEnergy;
  potentialEnergy       = rhs.potentialEnergy;
  neEnergy              = rhs.neEnergy;
  eeEnergy              = rhs.eeEnergy;
  logWeights            = rhs.logWeights;
  acceptanceProbability = rhs.acceptanceProbability;
  distanceMovedAccepted = rhs.distanceMovedAccepted;
  distanceMovedTrial    = rhs.distanceMovedTrial;
}

QMCProperties QMCProperties::operator + ( QMCProperties &rhs )
{
  QMCProperties result;
  
  matchParametersTo(rhs);
  result.matchParametersTo(rhs);
  
  result.weightChange          = weightChange + rhs.weightChange;
  result.walkerAge             = walkerAge + rhs.walkerAge;
  result.growthRate            = growthRate + rhs.growthRate;

  result.energy                = energy + rhs.energy;
  result.kineticEnergy         = kineticEnergy + rhs.kineticEnergy;
  result.potentialEnergy       = potentialEnergy + rhs.potentialEnergy;
  result.neEnergy              = neEnergy + rhs.neEnergy;
  result.eeEnergy              = eeEnergy + rhs.eeEnergy;
  result.logWeights            = logWeights + rhs.logWeights;
  result.acceptanceProbability = acceptanceProbability + 
    rhs.acceptanceProbability;
  result.distanceMovedAccepted = distanceMovedAccepted + 
    rhs.distanceMovedAccepted;
  result.distanceMovedTrial    = distanceMovedTrial + rhs.distanceMovedTrial;

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

  // nuc-elec energy
  strm << "<NucElecEnergy>" << endl;
  neEnergy.toXML(strm);
  strm << "</NucElecEnergy>" << endl;

  // elec-elec energy;
  strm << "<ElecElecEnergy>" << endl;
  eeEnergy.toXML(strm);
  strm << "</ElecElecEnergy>" << endl;

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

  strm << "<walkerAge>" << endl;
  walkerAge.toXML(strm);
  strm << "</walkerAge>" << endl;

  strm << "<weightChange>" << endl;
  weightChange.toXML(strm);
  strm << "</weightChange>" << endl;

  strm << "<growthRate>" << endl;
  growthRate.toXML(strm);
  strm << "</growthRate>" << endl;

  // Chi Density
  if (calc_density)
    {
      strm << "<ChiDensity>" << endl;
      for (int i=0; i<nBasisFunc; i++)
        chiDensity(i).toXML(strm);
      strm << "</ChiDensity>" << endl;
    }

  // Close XML
  strm << "</QMCProperties>" << endl;
}

bool QMCProperties::readXML(istream& strm)
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

  // Read nuc-elec energy
  strm >> temp;
  neEnergy.readXML(strm);
  strm >> temp;

  // Read elec-elec energy
  strm >> temp;
  eeEnergy.readXML(strm);
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

  strm >> temp;
  walkerAge.readXML(strm);
  strm >> temp;

  strm >> temp;
  weightChange.readXML(strm);
  strm >> temp;

  strm >> temp;
  growthRate.readXML(strm);
  strm >> temp;

  // Chi Density
  if (calc_density)
    {
      strm >> temp;
      for (int i=0; i<nBasisFunc; i++)
        chiDensity(i).readXML(strm);
      strm >> temp;
    }
    
  // Close XML
  strm >> temp;
  if(temp != "</QMCProperties>")
    {
      clog << "Error: checkpoint read failed in QMCProperties. We expected a \"</QMCProperties>\"" <<
	" tag, but found \"" << temp << "\"." << endl;
      return false;
    }
  return true;
}

ostream& operator <<(ostream& strm, QMCProperties &rhs)
{
  strm << endl << "------------------- Energy -------------------" << endl;
  strm << rhs.energy;

  strm << endl << "--------------- Kinetic Energy ---------------" << endl;
  strm << rhs.kineticEnergy;

  strm << endl << "-------------- Potential Energy --------------" << endl;
  strm << rhs.potentialEnergy;

  strm << endl << "-------------- Nuc-Elec Energy ---------------" << endl;
  strm << rhs.neEnergy;

  strm << endl << "-------------- Elec-Elec Energy --------------" << endl;
  strm << rhs.eeEnergy;

  strm << endl << "------------ AcceptanceProbability -----------" << endl;
  strm << rhs.acceptanceProbability;

  strm << endl << "------------ DistanceMovedAccepted -----------" << endl;
  strm << rhs.distanceMovedAccepted;

  strm << endl << "------------- DistanceMovedTrial -------------" << endl;
  strm << rhs.distanceMovedTrial;

  strm << endl << "----------------- Walker Age -----------------" << endl;
  strm << rhs.walkerAge;

  strm << endl << "--------------- Weight Change ----------------" << endl;
  strm << rhs.weightChange;

  strm << endl << "---------------- Growth Rate -----------------" << endl;
  strm << rhs.growthRate;

  strm << endl << "----------------- logWeights -----------------" << endl;
  strm << rhs.logWeights;

  return strm;
}

#ifdef PARALLEL

bool QMCProperties::mpiTypeCreated = false;

void QMCProperties::buildMpiReduce()
{
  MPI_Op_create((MPI_User_function*)Reduce_Function,true,&MPI_REDUCE);
}

void QMCProperties::buildMpiType()
{
  QMCProperties indata;

  // The number of properties 
  // ADJUST THIS WHEN ADDING NEW PROPERTIES
  const int NumberOfProperties = 12;
  
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
  MPI_Address(&(indata.neEnergy), &addresses[8]);
  MPI_Address(&(indata.eeEnergy), &addresses[9]);
  MPI_Address(&(indata.walkerAge), &addresses[10]);
  MPI_Address(&(indata.weightChange), &addresses[11]);
  MPI_Address(&(indata.growthRate), &addresses[12]);
  
  // Find the relative addresses of the data elements to the start of 
  // the struct
  for(int i=0; i<NumberOfProperties; i++)
    {  
      displacements[i] = addresses[i+1] - addresses[0];
    }

  /*
  MPI_Type_struct(NumberOfProperties, block_lengths, displacements, typelist, 
		  &MPI_TYPE);
  //Because we're including all kinds of junk in QMCProperties, the displacements
  //are not going to find everything. This ensures that MPI knows the exact size
  */
  
  MPI_Datatype temp;
  MPI_Type_struct(NumberOfProperties, block_lengths, displacements, typelist, 
                  &temp);
  MPI_Type_create_resized(temp,0,sizeof(QMCProperties),&MPI_TYPE);


  MPI_Type_commit(&MPI_TYPE);
}

MPI_Datatype QMCProperties::MPI_TYPE;

MPI_Op QMCProperties::MPI_REDUCE;

void QMCProperties::Reduce_Function(QMCProperties *in, QMCProperties *inout, 
                                   int *len, MPI_Datatype *dptr)
{
  /*
    we only want to add some of the elements (only the ones who's size is known
    at compile time). the others will be collected in other ways. we'll
    leave the operator + for a more general purpose.
  */
  //cout << globalInput.flags.my_rank << ":" << __FILE__ << ":" << __FUNCTION__ << ":" << __LINE__ << " len " << *len << endl;
 
  for(int i=0; i < *len; i++)
    {
      //inout[i] = inout[i] + in[i];
      
      inout[i].energy                = inout[i].energy + in[i].energy;
      inout[i].kineticEnergy         = inout[i].kineticEnergy + in[i].kineticEnergy;
      inout[i].potentialEnergy       = inout[i].potentialEnergy + in[i].potentialEnergy;
      inout[i].neEnergy              = inout[i].neEnergy + in[i].neEnergy;
      inout[i].eeEnergy              = inout[i].eeEnergy + in[i].eeEnergy;
      inout[i].logWeights            = inout[i].logWeights + in[i].logWeights;
      inout[i].acceptanceProbability = inout[i].acceptanceProbability + 
	in[i].acceptanceProbability;
      inout[i].distanceMovedAccepted = inout[i].distanceMovedAccepted + 
	in[i].distanceMovedAccepted;
      inout[i].distanceMovedTrial    = inout[i].distanceMovedTrial + in[i].distanceMovedTrial;

      inout[i].walkerAge             = inout[i].walkerAge + in[i].walkerAge;
      inout[i].weightChange          = inout[i].weightChange + in[i].weightChange;
      inout[i].growthRate            = inout[i].growthRate + in[i].growthRate;
    }
}

#endif

