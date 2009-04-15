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

#include "QMCInput.h"
#include "QMCDerivativeProperties.h"

QMCProperties::QMCProperties()
{  
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
}

void QMCProperties::zeroOut()
{
  walkerAge.zeroOut();
  weightChange.zeroOut();
  growthRate.zeroOut();

  energy.zeroOut();
  energy2.zeroOut();
  kineticEnergy.zeroOut();
  potentialEnergy.zeroOut();
  neEnergy.zeroOut();
  eeEnergy.zeroOut();
  x2.zeroOut();
  y2.zeroOut();
  z2.zeroOut();
  logWeights.zeroOut();
  acceptanceProbability.zeroOut();
  distanceMovedAccepted.zeroOut();
  distanceMovedTrial.zeroOut();
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

  if(newProperties->energy.getNumberSamples() > 0)
    {
      energy.newSample(newProperties->energy.getAverage(), weight);

      double e = newProperties->energy.getAverage();
      energy2.newSample(e*e, weight);
      
      kineticEnergy.newSample(newProperties->kineticEnergy.getAverage(), weight);
      potentialEnergy.newSample(newProperties->potentialEnergy.getAverage(),
				weight);
      neEnergy.newSample(newProperties->neEnergy.getAverage(), weight);
      eeEnergy.newSample(newProperties->eeEnergy.getAverage(), weight);
      x2.newSample(newProperties->x2.getAverage(),weight);
      y2.newSample(newProperties->y2.getAverage(),weight);
      z2.newSample(newProperties->z2.getAverage(),weight);
      acceptanceProbability.newSample
	(newProperties->acceptanceProbability.getAverage(), weight);
      distanceMovedAccepted.newSample
	(newProperties->distanceMovedAccepted.getAverage(), weight);
      distanceMovedTrial.newSample(newProperties->distanceMovedTrial.getAverage(),
				   weight);
      logWeights.newSample(newProperties->logWeights.getAverage(), nwalkers);
    }
}

void QMCProperties::operator = ( const QMCProperties &rhs )
{
  walkerAge             = rhs.walkerAge;
  weightChange          = rhs.weightChange;
  growthRate            = rhs.growthRate;

  energy                = rhs.energy;
  energy2               = rhs.energy2;
  kineticEnergy         = rhs.kineticEnergy;
  potentialEnergy       = rhs.potentialEnergy;
  neEnergy              = rhs.neEnergy;
  eeEnergy              = rhs.eeEnergy;
  x2                    = rhs.x2;
  y2                    = rhs.y2;
  z2                    = rhs.z2;
  logWeights            = rhs.logWeights;
  acceptanceProbability = rhs.acceptanceProbability;
  distanceMovedAccepted = rhs.distanceMovedAccepted;
  distanceMovedTrial    = rhs.distanceMovedTrial;
}

QMCProperties QMCProperties::operator + ( QMCProperties &rhs )
{
  QMCProperties result;
    
  result.weightChange          = weightChange + rhs.weightChange;
  result.walkerAge             = walkerAge + rhs.walkerAge;
  result.growthRate            = growthRate + rhs.growthRate;

  result.energy                = energy + rhs.energy;
  result.energy2               = energy2 + rhs.energy2;
  result.kineticEnergy         = kineticEnergy + rhs.kineticEnergy;
  result.potentialEnergy       = potentialEnergy + rhs.potentialEnergy;
  result.neEnergy              = neEnergy + rhs.neEnergy;
  result.eeEnergy              = eeEnergy + rhs.eeEnergy;
  result.x2                    = x2 + rhs.x2;
  result.y2                    = y2 + rhs.y2;
  result.z2                    = z2 + rhs.z2;
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

  if (globalInput.flags.checkpoint_energy_only == 0)
    {
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

      strm << "<x2>" << endl;
      x2.toXML(strm);
      strm << "</x2>" << endl;

      strm << "<y2>" << endl;
      y2.toXML(strm);
      strm << "</y2>" << endl;

      strm << "<z2>" << endl;
      z2.toXML(strm);
      strm << "</z2>" << endl;

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

      // Walker Age
      strm << "<walkerAge>" << endl;
      walkerAge.toXML(strm);
      strm << "</walkerAge>" << endl;
      
      // Weight Change
      strm << "<weightChange>" << endl;
      weightChange.toXML(strm);
      strm << "</weightChange>" << endl;

      // Growth Rate
      strm << "<growthRate>" << endl;
      growthRate.toXML(strm);
      strm << "</growthRate>" << endl;
    }

  // Close XML
  strm << "</QMCProperties>" << endl;
}

bool QMCProperties::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;
  if (temp != "<QMCProperties>")
    return false;

  // Read energy
  strm >> temp;
  if (temp != "<Energy>")
    return false;
  if (!energy.readXML(strm))
    return false;
  strm >> temp;
  if (temp != "</Energy>")
    return false;

  if (globalInput.flags.checkpoint_energy_only == 0)
    {
      // Read kinetic energy
      strm >> temp;
      if (temp != "<KineticEnergy>")
	return false;
      if (!kineticEnergy.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</KineticEnergy>")
	return false;

      // Read potential energy
      strm >> temp;
      if (temp != "<PotentialEnergy>")
	return false;
      if (!potentialEnergy.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</PotentialEnergy>")
	return false;

      // Read nuc-elec energy
      strm >> temp;
      if (temp != "<NucElecEnergy>")
	return false;
      if (!neEnergy.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</NucElecEnergy>")
	return false;

      // Read elec-elec energy
      strm >> temp;
      if (temp != "<ElecElecEnergy>")
	return false;
      if (!eeEnergy.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</ElecElecEnergy>")
	return false;

      strm >> temp;
      if (temp != "<x2>")
	return false;
      if (!x2.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</x2>")
	return false;

      strm >> temp;
      if (temp != "<y2>")
	return false;
      if (!y2.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</y2>")
	return false;

      strm >> temp;
      if (temp != "<z2>")
	return false;
      if (!z2.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</z2>")
	return false;

      // Read log weights
      strm >> temp;
      if (temp != "<LogWeights>")
	return false;
      if (!logWeights.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</LogWeights>")
	return false;

      // Acceptance probability
      strm >> temp;
      if (temp != "<AcceptanceProbability>")
	return false;
      if (!acceptanceProbability.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</AcceptanceProbability>")
	return false;

      // Distance Moved Accepted
      strm >> temp;
      if (temp != "<DistanceMovedAccepted>")
	return false;
      if (!distanceMovedAccepted.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</DistanceMovedAccepted>")
	return false;

      // Distance Moved Trial
      strm >> temp;
      if (temp != "<DistanceMovedTrial>")
	return false;
      if (!distanceMovedTrial.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</DistanceMovedTrial>")
	return false;

      // Walker Age
      strm >> temp;
      if (temp != "<walkerAge>")
	return false;
      if (!walkerAge.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</walkerAge>")
	return false;

      // Weight Change
      strm >> temp;
      if (temp != "<weightChange>")
	return false;
      if (!weightChange.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</weightChange>")
	return false;

      // Growth Rate
      strm >> temp;
      if (temp != "<growthRate>")
	return false;
      if (!growthRate.readXML(strm))
	return false;
      strm >> temp;
      if (temp != "</growthRate>")
	return false;
    }

  // Close XML
  strm >> temp;
  if(temp != "</QMCProperties>")
    return false;

  return true;
}

ostream& operator <<(ostream& strm, QMCProperties &rhs)
{
  strm << endl << "------------------- Energy -------------------" << endl;
  rhs.energy.printAll(strm);

  strm << endl << "----------------- Energy^2 -------------------" << endl;
  strm << rhs.energy2;

  strm << endl << "--------------- Kinetic Energy ---------------" << endl;
  strm << rhs.kineticEnergy;

  strm << endl << "-------------- Potential Energy --------------" << endl;
  strm << rhs.potentialEnergy;

  strm << endl << "-------------- Nuc-Elec Energy ---------------" << endl;
  strm << rhs.neEnergy;

  strm << endl << "-------------- Elec-Elec Energy --------------" << endl;
  strm << rhs.eeEnergy;

  strm << endl << "-------------------- x2 ----------------------" << endl;
  strm << rhs.x2;

  strm << endl << "-------------------- y2 ----------------------" << endl;
  strm << rhs.y2;

  strm << endl << "-------------------- z2 ----------------------" << endl;
  strm << rhs.z2;

  strm << endl << "------------ AcceptanceProbability -----------" << endl;
  strm << rhs.acceptanceProbability;

  strm << endl << "------------ DistanceMovedAccepted -----------" << endl;
  strm << rhs.distanceMovedAccepted;

  strm << endl << "------------- DistanceMovedTrial -------------" << endl;
  strm << rhs.distanceMovedTrial;

  strm << endl << "----------------- Walker Age -----------------" << endl;
  strm << rhs.walkerAge;

  if(globalInput.flags.run_type == "diffusion"){
    strm << endl << "--------------- Weight Change ----------------" << endl;
    strm << rhs.weightChange;
    
    strm << endl << "---------------- Growth Rate -----------------" << endl;
    strm << rhs.growthRate;

    strm << endl << "----------------- logWeights -----------------" << endl;
    strm << rhs.logWeights;
  }

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
  const int NumberOfProperties = 16;
  
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
  MPI_Address(&(indata.x2), &addresses[10]);
  MPI_Address(&(indata.y2), &addresses[11]);
  MPI_Address(&(indata.z2), &addresses[12]);
  MPI_Address(&(indata.walkerAge), &addresses[13]);
  MPI_Address(&(indata.weightChange), &addresses[14]);
  MPI_Address(&(indata.growthRate), &addresses[15]);
  MPI_Address(&(indata.energy2), &addresses[16]);
  
  // Find the relative addresses of the data elements to the start of 
  // the struct
  for(int i=0; i<NumberOfProperties; i++)
    {  
      displacements[i] = addresses[i+1] - addresses[0];
    }


#ifdef QMC_OLDMPICH
  /*
    If you're having trouble with MPICH (e.g. it can't find MPI_Type_create_resized)
    then compile with the USE_OLDMPICH flag.
  */
  MPI_Type_struct(NumberOfProperties, block_lengths, displacements, typelist, 
		  &MPI_TYPE);
#else
  /*
    Because we're including all kinds of junk in QMCProperties, the displacements
    are not going to find everything. This ensures that MPI knows the exact size.

    I've only MPICH version 1.2.* complain about the following lines, so if you want
    to use it, then you'll have to compile this file and QMCProperty with the
    QMC_OLDMPICH macro enabled.
  */
  MPI_Datatype temp;
  MPI_Type_struct(NumberOfProperties, block_lengths, displacements, typelist, 
                  &temp);
  MPI_Type_create_resized(temp,0,sizeof(QMCProperties),&MPI_TYPE);
#endif

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
      inout[i].x2                    = inout[i].x2 + in[i].x2;
      inout[i].y2                    = inout[i].y2 + in[i].y2;
      inout[i].z2                    = inout[i].z2 + in[i].z2;
      inout[i].logWeights            = inout[i].logWeights + in[i].logWeights;
      inout[i].acceptanceProbability = inout[i].acceptanceProbability + 
	in[i].acceptanceProbability;
      inout[i].distanceMovedAccepted = inout[i].distanceMovedAccepted + 
	in[i].distanceMovedAccepted;
      inout[i].distanceMovedTrial    = inout[i].distanceMovedTrial + in[i].distanceMovedTrial;

      inout[i].walkerAge             = inout[i].walkerAge + in[i].walkerAge;
      inout[i].weightChange          = inout[i].weightChange + in[i].weightChange;
      inout[i].growthRate            = inout[i].growthRate + in[i].growthRate;
      inout[i].energy2               = inout[i].energy2 + in[i].energy2;
    }
}

#endif

