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

#include "QMCFutureWalkingProperties.h"

int QMCFutureWalkingProperties::numFutureWalking = -1;

QMCFutureWalkingProperties::QMCFutureWalkingProperties()
{
  /*
    we need to be sure that globalInput has been initialized
    because we need to know how many future walking properties
    we are going to be collecting.
    
    globalInput is a global variable declared in QMCBeaver. i
    did this because it would be extremely tedious and prone
    to error to modify all instances of QMCFutureWalkingProperties to
    require a QMCInput to it's constructor.
    
    lastly, the globalInput needs to be initiated before
    any QMCFutureWalkingProperties are initiated because if we are running
    a parallel calculation, then we want the QMCFutureWalkingProperties on
    all the processors to have the right (static) MPI datatype.
    therefore, this information needs to be available before
    a QMCManager is created.
    
    this information is only available after the input file
    has been read in.
  */
  if(globalInput.flags.input_file_name == "")
  {
    cerr << "Error: The global QMCInput needs to have been initiated by this point!\n";
    exit(1);
  }
  
  //the extra one is to ensure we have at least one
  numFutureWalking = globalInput.flags.future_walking.size();
  
  fwEnergy.allocate(numFutureWalking);
  fwKineticEnergy.allocate(numFutureWalking);
  fwPotentialEnergy.allocate(numFutureWalking);
  r12.allocate(numFutureWalking);
  r2.allocate(numFutureWalking);
  
  calc_forces = false;
  zeroOut();
}

QMCFutureWalkingProperties::~QMCFutureWalkingProperties()
{
  if(calc_forces)
  {
    for (int i=0; i<nuclearForces.dim1(); i++)
      nuclearForces(i).deallocate();
    nuclearForces.deallocate();    
  }

  fwEnergy.deallocate();
  fwKineticEnergy.deallocate();
  fwPotentialEnergy.deallocate();
  r12.deallocate();
  r2.deallocate();
}

void QMCFutureWalkingProperties::setCalcForces(bool calcForces, int dim1, int dim2)
{
  calc_forces = calcForces;
  
  if(calc_forces == true)
    {
      nuclearForces.allocate(numFutureWalking);
      for (int fw=0; fw<nuclearForces.dim1(); fw++)
	{
	  nuclearForces(fw).allocate(dim1,dim2);
	  for(int i=0; i<nuclearForces(fw).dim1(); i++)
	    for(int j=0; j<nuclearForces(fw).dim2(); j++)
	      (nuclearForces(fw))(i,j).zeroOut();
	}
    }
  else
    {
      for (int i=0; i<nuclearForces.dim1(); i++)
	nuclearForces(i).deallocate();
      nuclearForces.deallocate();
    }
}

void QMCFutureWalkingProperties::zeroOut()
{
  for(int fw=0; fw<numFutureWalking; fw++)
  {
    fwEnergy(fw).zeroOut();
    fwKineticEnergy(fw).zeroOut();
    fwPotentialEnergy(fw).zeroOut();
    r12(fw).zeroOut();
    r2(fw).zeroOut();
    
    if (calc_forces)
      for (int i=0; i<nuclearForces(fw).dim1(); i++)
        for (int j=0; j<nuclearForces(fw).dim2(); j++)
          (nuclearForces(fw))(i,j).zeroOut();
  }
}

void QMCFutureWalkingProperties::newSample(QMCFutureWalkingProperties* newProperties, double weight, 
			      int nwalkers)
{  
  for(int fw=0; fw<numFutureWalking; fw++)
  {
    if(newProperties->fwKineticEnergy(fw).getNumberSamples() <= 0)
      continue;
      
    if (calc_forces)
    {
      for (int i=0; i<nuclearForces(fw).dim1(); i++)
        for (int j=0; j<nuclearForces(fw).dim2(); j++)
          (nuclearForces(fw))(i,j).newSample((newProperties->nuclearForces(fw))(i,j).getAverage(), 
              weight);
    }
    
    fwEnergy(fw).newSample(newProperties->fwEnergy(fw).getAverage(), weight);
    fwKineticEnergy(fw).newSample(newProperties->fwKineticEnergy(fw).getAverage(), weight);
    fwPotentialEnergy(fw).newSample(newProperties->fwPotentialEnergy(fw).getAverage(), weight);
    r12(fw).newSample(newProperties->r12(fw).getAverage(), weight);
    r2(fw).newSample(newProperties->r2(fw).getAverage(), weight);
    
  }
}

void QMCFutureWalkingProperties::matchParametersTo( const QMCFutureWalkingProperties &rhs )
{
  if(rhs.calc_forces)
    setCalcForces(rhs.calc_forces,rhs.nuclearForces.get(0).dim1(),rhs.nuclearForces.get(0).dim2());
}

void QMCFutureWalkingProperties::operator = ( const QMCFutureWalkingProperties &rhs )
{
  matchParametersTo(rhs);
  
  //Future walking collectors
  fwEnergy              = rhs.fwEnergy;
  fwKineticEnergy       = rhs.fwKineticEnergy;
  fwPotentialEnergy     = rhs.fwPotentialEnergy;
  r12                   = rhs.r12;
  r2                    = rhs.r2;

  if (calc_forces)
  {
    nuclearForces = rhs.nuclearForces;
  }
}

QMCFutureWalkingProperties QMCFutureWalkingProperties::operator + ( QMCFutureWalkingProperties &rhs )
{
  QMCFutureWalkingProperties result;
  
  matchParametersTo(rhs);
  result.matchParametersTo(rhs);
    
  //Future walking collectors/
  for(int fw=0; fw<numFutureWalking; fw++)
  {
    if (calc_forces)
    {
      for (int i=0; i<nuclearForces(fw).dim1(); i++)
        for (int j=0; j<nuclearForces(fw).dim2(); j++)
          (result.nuclearForces(fw))(i,j) = (nuclearForces(fw))(i,j)
                                          + (rhs.nuclearForces(fw))(i,j);
    }
    
    result.fwEnergy(fw)          = fwEnergy(fw) + rhs.fwEnergy(fw);
    result.fwKineticEnergy(fw)   = fwKineticEnergy(fw) + rhs.fwKineticEnergy(fw);
    result.fwPotentialEnergy(fw) = fwPotentialEnergy(fw) + rhs.fwPotentialEnergy(fw);
    result.r12(fw)               = r12(fw) + rhs.r12(fw);
    result.r2(fw)                = r2(fw) + rhs.r2(fw);
  }

  return result;
}

void QMCFutureWalkingProperties::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCFutureWalkingProperties>" << endl;

  // Nuclear forces
  if (calc_forces)
  {
    for(int fw=0; fw<numFutureWalking; fw++)
    {
      stringstream temp;
      temp << globalInput.flags.future_walking[fw];
      strm << "<NuclearForces" << temp.str() << ">" << endl;
      for (int i=0; i<nuclearForces(fw).dim1(); i++)
        for (int j=0; j<nuclearForces(fw).dim2(); j++)
          (nuclearForces(fw))(i,j).toXML(strm);
      strm << "</NuclearForces" << temp.str() << ">" << endl;
    }
  }
    
  for(int fw=0; fw<numFutureWalking; fw++)
  {
    stringstream temp;
    temp << globalInput.flags.future_walking[fw];

    strm << "<FWEnergy" << temp.str() << ">" << endl;
    fwEnergy(fw).toXML(strm);
    strm << "</FWEnergy" << temp.str() << ">" << endl;
    
    strm << "<FWKineticEnergy" << temp.str() << ">" << endl;
    fwKineticEnergy(fw).toXML(strm);
    strm << "</FWKineticEnergy" << temp.str() << ">" << endl;
    
    strm << "<FWPotentialEnergy" << temp.str() << ">" << endl;
    fwPotentialEnergy(fw).toXML(strm);
    strm << "</FWPotentialEnergy" << temp.str() << ">" << endl;
    
    strm << "<R12" << temp.str() << ">" << endl;
    r12(fw).toXML(strm);
    strm << "</R12" << temp.str() << ">" << endl;
    
    strm << "<R2" << temp.str() << ">" << endl;
    r2(fw).toXML(strm);
    strm << "</R2" << temp.str() << ">" << endl;
  }

  // Close XML
  strm << "</QMCFutureWalkingProperties>" << endl;
}

void QMCFutureWalkingProperties::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;
    
  // Nuclear forces
  if (calc_forces)
  {
    for(int fw=0; fw<numFutureWalking; fw++)
    {
      strm >> temp;
      for (int i=0; i<nuclearForces(fw).dim1(); i++)
        for (int j=0; j<nuclearForces(fw).dim2(); j++)
          (nuclearForces(fw))(i,j).readXML(strm);
      strm >> temp;
    }
  }
    
  for(int fw=0; fw<numFutureWalking; fw++)
  {
    strm >> temp;
    fwEnergy(fw).readXML(strm);
    strm >> temp;

    strm >> temp;
    fwKineticEnergy(fw).readXML(strm);
    strm >> temp;
    
    strm >> temp;
    fwPotentialEnergy(fw).readXML(strm);
    strm >> temp;
    
    strm >> temp;
    r12(fw).readXML(strm);
    strm >> temp;
    
    strm >> temp;
    r2(fw).readXML(strm);
    strm >> temp;
  }
    
  // Close XML
  strm >> temp;
}

ostream& operator <<(ostream& strm, QMCFutureWalkingProperties &rhs)
{ 
  int w1 = 10;
  int w2 = 10;
  int p1 = 2;

  strm << endl << "-------------- FW Total Energy ---------------" << endl;
  for(int fw=0; fw<rhs.numFutureWalking; fw++)
    {
      strm.width(w1);
      strm.precision(p1);
      strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
      strm.width(w2);
      strm << globalInput.flags.future_walking[fw]
	   << ": " << rhs.fwEnergy(fw);
    }

  strm << endl << "------------ FW Kinetic Energy ---------------" << endl;
  for(int fw=0; fw<rhs.numFutureWalking; fw++)
    {
      strm.width(w1);
      strm.precision(p1);
      strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
      strm.width(w2);
      strm << globalInput.flags.future_walking[fw]
	   << ": " << rhs.fwKineticEnergy(fw);
    }
  
  strm << endl << "----------- FW Potential Energy --------------" << endl;
  for(int fw=0; fw<rhs.numFutureWalking; fw++)
    {
      strm.width(w1);
      strm.precision(p1);
      strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
      strm.width(w2);
      strm << globalInput.flags.future_walking[fw]
	   << ": " << rhs.fwPotentialEnergy(fw);
    }
  
  strm << endl << "----------------- FW R12 ---------------------" << endl;
  for(int fw=0; fw<rhs.numFutureWalking; fw++)
    {
      strm.width(w1);
      strm.precision(p1);
      strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
      strm.width(w2);
      strm << globalInput.flags.future_walking[fw]
	   << ": " << rhs.r12(fw);
    }
  
  strm << endl << "------------------ FW R2 ---------------------" << endl;
  for(int fw=0; fw<rhs.numFutureWalking; fw++)
    {
      strm.width(w1);
      strm.precision(p1);
      strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
      strm.width(w2);
      strm << globalInput.flags.future_walking[fw]
	   << ": " << rhs.r2(fw);
    }
  
  // Nuclear forces
  if (rhs.calc_forces)
    {
      strm << endl << "------------ FW Nuclear Forces ---------------" << endl;
      for(int fw=0; fw<rhs.numFutureWalking; fw++)
	{
	  for (int i=0; i<rhs.nuclearForces(fw).dim1(); i++)
	    {
	      for (int j=0; j<rhs.nuclearForces(fw).dim2(); j++)
		{
		  
		  strm.width(5);
		  if(globalInput.flags.nuclear_derivatives == "bin_force_density")
		    {
		      strm << i;
		    } else {
		      strm << globalInput.Molecule.Atom_Labels(i) << i;
		      switch(j)
			{
			case 0: strm << "_x"; break;
			case 1: strm << "_y"; break;
			case 2: strm << "_z"; break;
			}
		    }
		  
		  strm.width(w1);
		  strm.precision(p1);
		  strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
		  strm.width(w2);
		  strm << globalInput.flags.future_walking[fw] 
		       << ": " << (rhs.nuclearForces(fw))(i,j);
		}
	    }
	}
    }
  
  return strm;
}

