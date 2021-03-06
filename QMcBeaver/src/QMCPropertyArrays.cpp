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

#include "QMCPropertyArrays.h"
#include "QMCDerivativeProperties.h"

int QMCPropertyArrays::numFutureWalking = -1;

QMCPropertyArrays::QMCPropertyArrays()
{
  /*
    we need to be sure that globalInput has been initialized
    because we need to know how many future walking properties
    we are going to be collecting.
    
    globalInput is a global variable declared in QMCBeaver. i
    did this because it would be extremely tedious and prone
    to error to modify all instances of QMCPropertyArrays to
    require a QMCInput to it's constructor.
    
    lastly, the globalInput needs to be initiated before
    any QMCPropertyArrays are initiated because if we are running
    a parallel calculation, then we want the QMCPropertyArrays on
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

  nBasisFunc   = 0;  
  calc_density = false;
  calc_forces = false;
  zeroOut();
}

QMCPropertyArrays::~QMCPropertyArrays()
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

  der.deallocate();
  for(int i=0; i<hess.dim1(); i++)
    hess(i).deallocate();
  hess.deallocate();

  chiDensity.deallocate();
}

void QMCPropertyArrays::setCalcDensity(bool calcDensity, int nbasisfunctions)
{
  calc_density = calcDensity;

  if(calc_density == true)
    {
      nBasisFunc   = nbasisfunctions;
      chiDensity.allocate(nBasisFunc);
    }
  else
    {
      nBasisFunc = 0;
      chiDensity.deallocate();
    }
}

void QMCPropertyArrays::setCalcForces(bool calcForces, int dim1, int dim2)
{
  calc_forces = calcForces;
  
  if(calc_forces == true)
    {
      nuclearForces.allocate(numFutureWalking);
      for (int fw=0; fw<nuclearForces.dim1(); fw++)
	nuclearForces(fw).allocate(dim1,dim2);
    }
  else
    {
      for (int i=0; i<nuclearForces.dim1(); i++)
	nuclearForces(i).deallocate();
      nuclearForces.deallocate();
    }
}

void QMCPropertyArrays::zeroOut()
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

  if(globalInput.cs_Parameters.dim1() > 1)
    cs_Energies.allocate(globalInput.cs_Parameters.dim1());
  else
    cs_Energies.deallocate();

  int numAI = globalInput.getNumberAIParameters();

  if(globalInput.flags.calculate_Derivatives != 0)
    {
      //second parameter refers to the number of terms necessary
      //to calculate derivatives of the optimization objective function
      der.allocate(numAI,5);
      der = 0.0;
    } else {
      der.deallocate();
    }

  if(globalInput.flags.calculate_Derivatives != 0)
    {
      if( globalInput.flags.optimize_Psi_criteria == "analytical_energy_variance")
	{
	  hessIsSymmetric = true;
	  hess.allocate(1);
	  hess(0).allocate(numAI,numAI);
	} else if(globalInput.flags.optimize_Psi_criteria == "generalized_eigenvector")
	  {
	    hessIsSymmetric = false;
	    hess.allocate(3);
	    for(int i=0; i<hess.dim1(); i++)
	      hess(i).allocate(numAI,numAI);	
	  } else {
	    clog << "Error: unknown optimize_Psi_criteria\n";
	  }

      for(int i=0; i<hess.dim1(); i++)
	hess(i) = 0.0;
    } else {
      for(int i=0; i<hess.dim1(); i++)
	hess(i).deallocate();
      hess.deallocate();
    }

  for(int i=0; i<cs_Energies.dim1(); i++)
    cs_Energies(i).zeroOut();

  numDerHessSamples = 0;

  if (calc_density)
    for (int i=0; i<nBasisFunc; i++)
      chiDensity(i).zeroOut();
}

void QMCPropertyArrays::newSample(QMCPropertyArrays* newProperties, double weight, 
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

  for(int i=0; i<cs_Energies.dim1(); i++)
    if(newProperties->cs_Energies(i).getNumberSamples() > 0)
      cs_Energies(i).newSample(newProperties->cs_Energies(i).getAverage(), weight);

  if (calc_density)
    for (int i=0; i<nBasisFunc; i++)
      if(newProperties->chiDensity(i).getNumberSamples() > 0)
	chiDensity(i).newSample(newProperties->chiDensity(i).getAverage(),weight);
}

void QMCPropertyArrays::matchParametersTo( const QMCPropertyArrays &rhs )
{
  if(rhs.calc_forces)
    setCalcForces(rhs.calc_forces,rhs.nuclearForces.get(0).dim1(),rhs.nuclearForces.get(0).dim2());
  if(rhs.calc_density)
    setCalcDensity(rhs.calc_density,rhs.nBasisFunc);

  zeroOut();
}

void QMCPropertyArrays::operator = ( const QMCPropertyArrays &rhs )
{
  matchParametersTo(rhs);
  
  //Future walking collectors
  props = rhs.props;

  if (calc_forces)
    {
      nuclearForces = rhs.nuclearForces;
    }
  if (calc_density)
    {
      chiDensity = rhs.chiDensity;
    }

  der = rhs.der;
  hess = rhs.hess;
  cs_Energies = rhs.cs_Energies;
  numDerHessSamples = rhs.numDerHessSamples;
}

QMCPropertyArrays QMCPropertyArrays::operator + ( QMCPropertyArrays &rhs )
{
  QMCPropertyArrays result;

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

void QMCPropertyArrays::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCPropertyArrays>" << endl;

  // Chi Density
  if (calc_density)
    {
      strm << "<ChiDensity>" << endl;
      for (int i=0; i<nBasisFunc; i++)
        chiDensity(i).toXML(strm);
      strm << "</ChiDensity>" << endl;
    }

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
  strm << "</QMCPropertyArrays>" << endl;
}

bool QMCPropertyArrays::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;

  // Chi Density
  if (calc_density)
    {
      strm >> temp;
      for (int i=0; i<nBasisFunc; i++)
        chiDensity(i).readXML(strm);
      strm >> temp;
    }
    
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
  if(temp != "</QMCPropertyArrays>")
    {
      clog << "Error: checkpoint read failed in QMCPropertyArrays."
	   << " We expected to read a \"</QMCPropertyArrays>\""
	   << " tag, but found \"" << temp << "\"." << endl;
      return false;
    }
  return true;
}

ostream& operator <<(ostream& strm, QMCPropertyArrays &rhs)
{ 
  int w1 = 10;
  int w2 = 10;
  int p1 = 2;

  if(rhs.cs_Energies.dim1() > 1)
    {
      strm << endl << "------ Correlated Sampling Energies --------" << endl;
      for(int i=0; i<rhs.cs_Energies.dim1(); i++)
	{
	  strm << "Set " << setw(3) << i << ": Sample variance = " << rhs.cs_Energies(i).getSeriallyCorrelatedVariance() << endl;
	  strm << rhs.cs_Energies(i);
	  //rhs.cs_Energies(i).printAll(strm);
	  strm << endl;
	}
    }

  //We almost certainly don't need these printed out
  //Just look for the derivative printout during the optimization step.
  if(false)
    {
      if(rhs.der.dim1() > 0)
	{
	  strm << endl << "------ Objective Function Derivatives --------" << endl;
	  for(int i=0; i<rhs.der.dim1(); i++)
	  {
	    for(int j=0; j<rhs.der.dim2(); j++)
	      {
		strm << "Param " << i << " term " << j << ": " << rhs.der(i,j);
	      }
	    strm << endl;
	  }
	}
    }
  
  if(false)
    {
      if(rhs.hess.dim1() > 0)
	{
	  strm << endl << "------ Objective Function Hessian --------" << endl;
	  
	  for (int i=0; i<rhs.hess.dim1(); i++)
	    {
	      for(int j=0; j<rhs.hess(i).dim1(); j++)
		{
		  int max = rhs.hess(i).dim2();
		  if(rhs.hessIsSymmetric)
		    max = j+1;
		  
		  for(int k=0; k<max; k++)
		    strm << "(" << j << "," << k << "): " << (rhs.hess(i))(j,k);
		}
	      strm << endl;
	    }
	  strm << endl << endl;
	}
    }

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

