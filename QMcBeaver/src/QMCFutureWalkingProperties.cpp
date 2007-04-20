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

  names.allocate(NUM_PROPS);
  //In order to read and write XML tags, spaces
  //are not allowed in names.
  names[FW_It]  = "Normalization_I(t)";
  names[FW_TE]  = "Total_Energy";
  names[FW_KE]  = "Kinetic_Energy";
  names[FW_KEg] = "Kinetic_Energy_g";
  names[FW_PE]  = "Potential_Energy";
  names[FW_R12] = "R12";
  names[FW_R2]  = "R2";
  names[FW_iR]  = "Inverse_R";
  names[FW_iR12]= "Inverse_R12";

  names[FW_TE_2]  = "Total_Energy^2";
  names[FW_KE_2]  = "Kinetic_Energy^2";
  names[FW_KEg_2] = "Kinetic_Energy_g^2";
  names[FW_PE_2]  = "Potential_Energy^2";
  names[FW_R12_2] = "R12^2";
  names[FW_R2_2]  = "R2^2";
  names[FW_iR_2]  = "Inverse_R^2";
  names[FW_iR12_2]= "Inverse_R12^2";

  props.allocate(NUM_PROPS);
  for(int i=0; i<props.size(); i++) 
    props(i).allocate(numFutureWalking);
  
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

  for(int i=0; i<props.size(); i++)
    props(i).deallocate();
  props.deallocate();
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
    for(int i=0; i<props.size(); i++)
      (props(i))(fw).zeroOut();
    
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
    if (calc_forces && (newProperties->nuclearForces(fw))(0,0).getNumberSamples() > 0)
      {
	for (int i=0; i<nuclearForces(fw).dim1(); i++)
	  for (int j=0; j<nuclearForces(fw).dim2(); j++)
	    (nuclearForces(fw))(i,j).newSample((newProperties->nuclearForces(fw))(i,j).getAverage(), 
					       weight);
      }
    
    for(int i=0; i<props.size(); i++)
      if((newProperties->props(i))(fw).getNumberSamples() > 0)
	(props(i))(fw).newSample( (newProperties->props(i))(fw).getAverage(), weight);
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
  props = rhs.props;

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

    for(int i=0; i<props.size(); i++)
      (result.props(i))(fw) = (props(i))(fw) + (rhs.props(i))(fw);
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

    for(int i=0; i<props.size(); i++)
      {
	//No spaces are allowed in a tag
	strm <<  "<FW_" << names[i] << "_" << temp.str() << ">" << endl;
	(props(i))(fw).toXML(strm);
	strm << "</FW_" << names[i] << "_" << temp.str() << ">" << endl;
      }
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
    for(int i=0; i<props.size(); i++)
      {
	strm >> temp;
	(props(i))(fw).readXML(strm);
	strm >> temp;
      }
  }
    
  // Close XML
  strm >> temp;
  if(temp != "</QMCFutureWalkingProperties>")
    {
      clog << "Error: checkpoint read failed in QMCFutureWalkingProperties."
	   << " We expected to read a \"</QMCFutureWalkingProperties>\""
	   << " tag, but found \"" << temp << "\"." << endl;
      exit(0);
    }
}

ostream& operator <<(ostream& strm, QMCFutureWalkingProperties &rhs)
{ 
  int w1 = 10;
  int w2 = 10;
  int p1 = 2;

  if(rhs.numFutureWalking <= 1) return strm;

  for(int i=0; i<rhs.props.size(); i++)
    {
      strm << endl << "-------------------------- FW " << rhs.names[i] << endl;
      for(int fw=0; fw<rhs.numFutureWalking; fw++)
	{
	  if((rhs.props(i))(fw).getNumberSamples() <= 0) break;
	  strm.width(w1);
	  strm.precision(p1);
	  strm << globalInput.flags.dt_effective * globalInput.flags.future_walking[fw] << " <=> ";
	  strm.width(w2);
	  strm << globalInput.flags.future_walking[fw]
	       << ": " << (rhs.props(i))(fw);
	}
    }
  
  // Nuclear forces
  if (rhs.calc_forces)
    {
      strm << endl << "------------ FW Nuclear Forces ---------------" << endl;
      for(int fw=0; fw<rhs.numFutureWalking; fw++)
	{
	  if((rhs.nuclearForces(fw))(0,0).getNumberSamples() <= 0) break;
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

