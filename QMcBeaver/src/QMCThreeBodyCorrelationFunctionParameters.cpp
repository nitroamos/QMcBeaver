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

#include "QMCThreeBodyCorrelationFunctionParameters.h"
#include "QMCInput.h"

Array1D<double> QMCThreeBodyCorrelationFunctionParameters::getFreeParameters()
{
  return freeParameters;
}

void QMCThreeBodyCorrelationFunctionParameters::operator = (const 
			       QMCThreeBodyCorrelationFunctionParameters & rhs)
{
  NumberOfParameterTypes = rhs.NumberOfParameterTypes;
  NumberOfParameters = rhs.NumberOfParameters;
  NeN = rhs.NeN;
  Nee = rhs.Nee;

  NumberOfFreeParameters = rhs.NumberOfFreeParameters;
  freeParameters = rhs.freeParameters;

  paramDepMatrix = rhs.paramDepMatrix;

  TotalNumberOfParameters = rhs.TotalNumberOfParameters;
  Parameters = rhs.Parameters;

  C = rhs.C;
  cutoff = rhs.cutoff;

  ParticleTypes = rhs.ParticleTypes;
  threeBodyCorrelationFunctionType = rhs.threeBodyCorrelationFunctionType;

  setThreeBodyCorrelationFunction();
  initializeThreeBodyCorrelationFunctionParameters();
}

/**
 *Load QMCThreeBodyCorrelationFunctionParameters of the form: <br>
 * ParticleTypes: He Electron_up Electron_down <br>
 * CorrelationFunctionType: Cambridge
 * NumberOfParameterTypes: 2 <br>
 * NumberOfParametersOfEachType: 2 2 <br>
 * Parameters: 2.0 3.0 4.0 1.0 2.0 1.0 3.0 5.0 <br>
 * C: 3 <br>
 * Cutoff: 2.5 <br>
 */

bool QMCThreeBodyCorrelationFunctionParameters::read(istream & strm)
{
  bool ok = true;
  string temp;
  
  // Read the particle types
  
  string pt1, pt2, pt3;
  strm >> temp >> pt1 >> pt2 >> pt3;
  pt1 = StringManipulation::toFirstUpperRestLower(pt1);
  pt2 = StringManipulation::toFirstUpperRestLower(pt2);
  pt3 = StringManipulation::toFirstUpperRestLower(pt3);

  // Order the particle types

  ParticleTypes.allocate(3);
  
  // The fist particle should be a nucleus, the the second two particles should
  // be electrons.  I am not going to go into all the reordering possibilities.
  // We should just make sure that the input files have the particles in the
  // right order.

  ParticleTypes(0) = pt1;
  ParticleTypes(1) = pt2;
  ParticleTypes(2) = pt3;

  if ( pt2 == "Electron_up")
    {
      if (pt3 != "Electron_up" && pt3 != "Electron_down")
	{
	  cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
	  cerr << "particles " << pt1 << ", " << pt2 << ", " << pt3 << " are ";
	  cerr << "being loaded as a three body correlation function!" << endl;
	  ok = false;
	}
    }
  else if ( pt2 == "Electron_down")
    {
      if (pt3 != "Electron_down")
	{ 
	  // If we have an opposite spin pair, we assume the up electron is 
	  // listed first.
	  cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
	  cerr << "particles " << pt1 << ", " << pt2 << ", " << pt3 << " are ";
	  cerr << "being loaded as a three body correlation function!" << endl;
	  ok = false;
	}
    }
  
  // Load the correlation function type
  
  strm >> temp >> threeBodyCorrelationFunctionType;
  if (threeBodyCorrelationFunctionType != "Cambridge" && 
      threeBodyCorrelationFunctionType != "None")
    {
      cerr << "ERROR: unknown type of ThreeBodyCorrelationFunction is being ";
      cerr << "read.  threeBodyCorrelationFunctionType = ";
      cerr << threeBodyCorrelationFunctionType << endl;
      ok = false;
    }

  // Load the number of parameter types
  
  strm >> temp >> temp;
  NumberOfParameterTypes = atoi(temp.c_str());

  if (NumberOfParameterTypes != 2)
    {
      // There should be two parameter types: NeN and Nee.
      cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
      cerr << "NumberOfParametersTypes = " << NumberOfParameterTypes;
      cerr << ", 2 was expected." << endl;
      ok = false;
    }

  // Load the number of parameters of each type
  // The first number is NeN, then Nee

  NumberOfParameters.allocate(NumberOfParameterTypes);
  strm >> temp;
  
  for (int i=0; i<NumberOfParameterTypes; i++)
    {
      strm >> temp;
      NumberOfParameters(i) = atoi(temp.c_str());
    }
  
  NeN = NumberOfParameters(0);
  Nee = NumberOfParameters(1);
  
  // Right now, we need Nee to be at least 3 in order to have enough parameters
  // to satisfy the EN cusp condition without screwing up the EE cusp condition
  if (NeN < 3 || Nee < 3)
    {
      cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
      cerr << "NeN = " << NeN << ", Nee = " << Nee << ", values of at least ";
      cerr << "3 are needed" << endl;
      ok = false;
    }

  // Set the total number of parameters

  TotalNumberOfParameters = NeN*NeN*Nee;

  // Load the parameters
  
  Parameters.allocate(TotalNumberOfParameters);
  strm >> temp;
  
  for (int i=0; i < TotalNumberOfParameters; i++)
    {
      strm >> temp;
      Parameters(i) = atof(temp.c_str());
    }
  
  // Load the value for C, the exponent of the prefactor.

  strm >> temp >> temp;
  C = atoi(temp.c_str());

  // Load the value for the electron-nucleus cutoff distance for this function
  strm >> temp >> temp;
  cutoff = atof(temp.c_str());

  // Make the matrix of parameter dependencies.  This matrix generates the set
  // of all parameters from the set of free parameters.  The parameters 
  // generated will satisfy symmetry and the cusp conditions
  makeParamDepMatrix();

  // Now we calculate how many of the parameters are free to be optimized.  
  // Some of the parameters are constrained by symmetry or cusp conditions.

  // If everything was free to be optmized, we would have all of the parameters
  // and the cutoff distance

  NumberOfFreeParameters = TotalNumberOfParameters+1; 

  if (globalInput.flags.optimize_NEE_cutoffs == 0)
    NumberOfFreeParameters -= 1;

  // We subtract the number of elements with m>l
  NumberOfFreeParameters -= Nee*NeN*(NeN-1)/2; 

  // For every k from 0 to 2*NeN-2, there is one dependent parameter due to the
  // electron-electron no-cusp condition
  NumberOfFreeParameters -= 2*NeN-1;

  // For every k from 0 to NeN+Nee-2, there is one dependent parameter due to
  // the electron-nucleus no-cusp condition
  NumberOfFreeParameters -= NeN+Nee-1;

  // If we don't allow terms that are proportional to the two particle
  // functions, the number of free parameters is reduced
  if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 0)
    {
      NumberOfFreeParameters -= NeN-1;
      for (int l=1; l<NeN; l++)
	Parameters(l*NeN*Nee) = 0.0;
    }

  if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 0)
    {
      NumberOfFreeParameters -= Nee-2;
      for (int n=2; n<Nee; n++)
	Parameters(n) = 0.0;
    }

  // if the correlation function type is "None"
  // then remove any parameters or constants
  // that were listed
  
  if (threeBodyCorrelationFunctionType == "None")
    {
      NumberOfParameterTypes = 0;
      TotalNumberOfParameters = 0;
      NumberOfParameters.deallocate();
      Parameters.deallocate();
      C = 0;
      cutoff = 0.0;
      NumberOfFreeParameters = 0;
      freeParameters.deallocate();
      paramDepMatrix.deallocate();
    }
  else
    {
      // We multiply the paramDepMatrix by the parameters we loaded in to make
      // sure the parameters satisfy the symmetry and cusp conditions

      Array1D<double> temp;
      temp.allocate(TotalNumberOfParameters);
      temp = 0.0;

      for (int i=0; i<TotalNumberOfParameters; i++)
	for (int j=0; j<TotalNumberOfParameters; j++)
	  temp(i) += paramDepMatrix(i,j)*Parameters(j);

      Parameters = temp;
      checkCuspAndSymmetry();
      makeFreeParameterArray();
    }

  // set the correlation function
  setThreeBodyCorrelationFunction();
  initializeThreeBodyCorrelationFunctionParameters();

  return ok;
}

QMCThreeBodyCorrelationFunctionParameters::QMCThreeBodyCorrelationFunctionParameters()
{
  threeBodyCorrelationFunctionType = "None";
  ThreeBodyCorrelationFunction = 0;

  TotalNumberOfParameters = 0;
  NumberOfFreeParameters = 0;

  NeN = 0;
  Nee = 0;
  C = 0;
  cutoff = 0.0;
}

QMCThreeBodyCorrelationFunctionParameters::QMCThreeBodyCorrelationFunctionParameters(const QMCThreeBodyCorrelationFunctionParameters & rhs)
{
  ThreeBodyCorrelationFunction = 0;
  *this = rhs;
}

string QMCThreeBodyCorrelationFunctionParameters::getParticle1Type()
{
  return ParticleTypes(0);
}

string QMCThreeBodyCorrelationFunctionParameters::getParticle2Type()
{
  return ParticleTypes(1);
}

string QMCThreeBodyCorrelationFunctionParameters::getParticle3Type()
{
  return ParticleTypes(2);
}

int QMCThreeBodyCorrelationFunctionParameters::getNumberOfFreeParameters()
{
  return NumberOfFreeParameters;
}

void QMCThreeBodyCorrelationFunctionParameters::setFreeParameters(Array1D<double> & params)
{
  if( params.dim1() != NumberOfFreeParameters )
    {
      cerr << "ERROR: Parameters of the incorrect size are trying to be set "
	   << "in QMCCorrelationFunctionParameters" << endl;
    }

  freeParameters = params;
  setParameters();

  Array1D<double> temp;
  temp.allocate(TotalNumberOfParameters);
  temp = 0.0;
  
  for (int i=0; i<TotalNumberOfParameters; i++)
    for (int j=0; j<TotalNumberOfParameters; j++)
      temp(i) += paramDepMatrix(i,j)*Parameters(j);
  
  Parameters = temp;

  checkCuspAndSymmetry();
  initializeThreeBodyCorrelationFunctionParameters();
}

void QMCThreeBodyCorrelationFunctionParameters::setParticle1Type(string val)
{
  ParticleTypes(0) = val;
}

void QMCThreeBodyCorrelationFunctionParameters::setParticle2Type(string val)
{
  ParticleTypes(1) = val;
}

void QMCThreeBodyCorrelationFunctionParameters::setParticle3Type(string val)
{
  ParticleTypes(2) = val;
}

double QMCThreeBodyCorrelationFunctionParameters::getCutoffDist()
{
  return cutoff;
}

void QMCThreeBodyCorrelationFunctionParameters::setThreeBodyCorrelationFunction()
{
  if( ThreeBodyCorrelationFunction != 0 )
    {
      delete ThreeBodyCorrelationFunction;
      ThreeBodyCorrelationFunction = 0;
    }
  ThreeBodyCorrelationFunction = QMCThreeBodyCorrelationFunctionFactory::
    threeBodyCorrelationFunctionFactory(threeBodyCorrelationFunctionType);
}

QMCThreeBodyCorrelationFunction * QMCThreeBodyCorrelationFunctionParameters::
                                             getThreeBodyCorrelationFunction()
{
  return ThreeBodyCorrelationFunction;
}

QMCThreeBodyCorrelationFunctionParameters::~QMCThreeBodyCorrelationFunctionParameters()
{
  NumberOfParameters.deallocate();
  Parameters.deallocate();
  freeParameters.deallocate();
  ParticleTypes.deallocate();
  
  if( ThreeBodyCorrelationFunction != 0 )
    {
      delete ThreeBodyCorrelationFunction;
      ThreeBodyCorrelationFunction = 0;
    }
}

void QMCThreeBodyCorrelationFunctionParameters::
                             initializeThreeBodyCorrelationFunctionParameters()
{
  ThreeBodyCorrelationFunction->initializeParameters(NeN,Nee,Parameters,C,
						     cutoff);
}

ostream & operator << (ostream & strm, 
		       QMCThreeBodyCorrelationFunctionParameters & rhs)
{
  strm << "ParticleTypes:\t" << rhs.ParticleTypes(0) << "\t" 
       << rhs.ParticleTypes(1) << "\t" << rhs.ParticleTypes(2) << endl;

  strm << "threeBodyCorrelationFunctionType:\t" 
       << rhs.threeBodyCorrelationFunctionType << endl;

  strm << "NumberOfParameterTypes:\t" << rhs.NumberOfParameterTypes << endl;

  strm << "NumberOfParametersOfEachType:";

  for (int i = 0; i < rhs.NumberOfParameterTypes; i++)
    strm << "  " << rhs.NumberOfParameters(i);
  strm << endl;

  strm << "Parameters:";

  for (int i = 0; i < rhs.TotalNumberOfParameters; i++)
    strm << "  " << rhs.Parameters(i);
  strm << endl;

  strm << "C:  " << rhs.C << endl;
  
  strm << "Cutoff:  " << rhs.cutoff << endl;

  strm << endl;

  return strm;
}

void QMCThreeBodyCorrelationFunctionParameters::makeFreeParameterArray()
{
  // This function takes the full set of parameters and makes the set of 
  // parameters that are free to be optimized.  The decision is based on the 
  // flags of what we are optimizing and the symmetry of the function.

  freeParameters.allocate(NumberOfFreeParameters);
  freeParameters = 0.0;

  int counter = 0;
  
  // If the cutoff distance is being optimized it is the 0th element of the
  // free parameter array
  if (globalInput.flags.optimize_NEE_cutoffs == 1)
    {
      freeParameters(0) = cutoff;
      counter++;
    }

  double element = 0.0;
  // Because of the symmetry due to interchange of the electrons, we only 
  // write the parameters that have l>=m.
  for (int l=0; l<NeN; l++)
    for (int m=0; m<=l; m++)
      for (int n=0; n<Nee; n++)
	{
	  element = Parameters(l*NeN*Nee+m*Nee+n);
	  if (l>0 && m==0 && n==0)
	    {
	      if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 1)
		{
		  freeParameters(counter) = element;
		  counter++;
		}
	    }
	  else if (l==0 && n>1)
	    {
	      if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 1)
		{
		  freeParameters(counter) = element;
		  counter++;
		}
	    }
	  // These parameters are dependent due to the EE no-cusp condition
	  else if (m==0 && n==1)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  else if (l==NeN-1 && n==1)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  // These parameters are dependent due to the EN no-cusp condition
	  else if (l==0 && n==0)
	    {
	      // Do nothing- this parameter is dependent
	    }
	  else if (l<NeN && m==1 && n==0)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  else if (l==NeN-2 && m==1 && n==2)
	    {
	      // Do nothing- this parameter is dependent
	    }
	  else if (l==NeN-1 && m==1 && n>1)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  else
	    {
	      freeParameters(counter) = element;
	      counter++;
	    }
	}

  if (counter != NumberOfFreeParameters)
    {
      // We have miscalculated the number of free parameters.
      cerr << "ERROR in "
	   << "QMCThreeBodyCorrelationFunctionParameters::setFreeParameters()"
	   << ": there should be " << NumberOfFreeParameters << " free "
	   << "parameters, but we find " << counter+1 << endl;
      exit(0);
    }
}

void QMCThreeBodyCorrelationFunctionParameters::setParameters()
{
  // This function takes the set of free parameters and generates the set of 
  // all parameters to send to the correlation function.

  int counter = 0;
  Parameters = 0.0;

  // If the cutoff distance is being optimized it is the 0th element of the
  // free parameter array
  if (globalInput.flags.optimize_NEE_cutoffs == 1)
    {
      cutoff = freeParameters(0);
      counter++;
    }

  int index = -1;
  // We only consider the elements such that l>=m.  The others will be filled
  // in by symmetry
  for (int l=0; l<NeN; l++)
    for (int m=0; m<=l; m++)
      for (int n=0; n<Nee; n++)
	{
	  index = l*NeN*Nee+m*Nee+n;
	  if (l>0 && m==0 && n==0)
	    {
	      if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 1)
		{
		  Parameters(index) = freeParameters(counter);
		  counter++;
		}
	    }
	  else if (l==0 && n>1)
	    {
	      if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 1)
		{
		  Parameters(index) = freeParameters(counter);
		  counter++;
		}
	    }
	  // These parameters are dependent due to the EE no-cusp condition
	  else if (m==0 && n==1)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  else if (l==NeN-1 && n==1)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  // These parameters are dependent due to the EN no-cusp condition
	  else if (l==0 && n==0)
	    {
	      // Do nothing- this parameter is dependent
	    }
	  else if (l<NeN && m==1 && n==0)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  else if (l==NeN-2 && m==1 && n==2)
	    {
	      // Do nothing- this parameter is dependent
	    }
	  else if (l==NeN-1 && m==1 && n>1)
	    {
	      // Do nothing- these parameters are dependent
	    }
	  else 
	    {
	      Parameters(index) = freeParameters(counter);
	      counter++;
	    }
	}
}

void QMCThreeBodyCorrelationFunctionParameters::checkCuspAndSymmetry()
{
  // This function makes sure that the parameters satisfy the cusp and 
  // symmetry conditions.

  // First we will make sure the right elements equal each other
  // then we will run through all the sums to make sure they add to zero

  int index1 = -1;
  int index2 = -1;

  for (int l=0; l<NeN-1; l++)
    for (int m=l+1; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	{
	  index1 = l*NeN*Nee+m*Nee+n;
	  index2 = m*NeN*Nee+l*Nee+n;

	  if ( fabs(Parameters(index1)-Parameters(index2)) > 1e-10 )
	    {
	      cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
		   << "checkCuspAndSymmetry(): symmetry is not satisfied"
		   << endl;
	      exit(0);
	    }
	}

  double sum;
  
  // First we check the electron-electron cusp condition
  for (int k=0; k<2*NeN-1; k++)
    {
      sum = 0.0;
      for (int l=0; l<NeN; l++)
	for (int m=0; m<NeN; m++)
	  if (l+m == k)
	    {
	      index1 = l*NeN*Nee+m*Nee+1;
	      sum += Parameters(index1);
	    }
      if (fabs(sum) > 1e-10)
	{
	  cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
	       << "checkCuspAndSymmetry(): electron-electron cusp condition "
	       << "is not satisfied" << endl;
	  exit(0);
	}
    }

  // Now we check the electron-nucleus cusp condition
  for (int k=0; k<NeN+Nee-1; k++)
    {
      sum = 0.0;
      for (int m=0; m<NeN; m++)
	for (int n=0; n<Nee; n++)
	  if (m+n == k)
	    {
	      index1 = m*Nee+n;
	      index2 = NeN*Nee+m*Nee+n;
	      sum += C*Parameters(index1) - cutoff*Parameters(index2);
	    }
      if (fabs(sum) > 1e-10)
	{
	  cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
	       << "checkCuspAndSymmetry(): electron-nucleus cusp condition "
	       << "is not satisfied" << endl;
	  exit(0);
	}
    }
}

void QMCThreeBodyCorrelationFunctionParameters::makeParamDepMatrix()
{
  // This function will make the matrix of parameter dependencies.  The size of
  // the matrix will be TotalNumberOfParameters x TotalNumberOfParameters
  // The matrix depends on NeN, Nee, and the cusp conditions.  It is made when
  // the function is read and not changed.

  paramDepMatrix.allocate(TotalNumberOfParameters,TotalNumberOfParameters);

  // We start off by assuming all the parameters are free
  for (int i=0; i<TotalNumberOfParameters; i++)
    for (int j=0; j<TotalNumberOfParameters; j++)
      {
	if (i != j)
	  paramDepMatrix(i,j) = 0.0;
	else
	  paramDepMatrix(i,j) = 1.0;
      }

  // We start with the electron-electron cusp

  // The parameter 001 is zero due to the electron-electron cusp (k=0)
  paramDepMatrix(1,1) = 0.0;

  // The parameter 101 is zero due to the electron-electron cusp (k=1)
  paramDepMatrix(NeN*Nee+1,NeN*Nee+1) = 0.0;

  int l_index   = -1;
  int m_index   = -1;
  int n_index   = -1;
  int index1    = -1;
  int index2    = -1;
  int index_dep = -1;

  for (int k=2; k<2*NeN-1; k++)
    {
      // We find the index of the dependent parameter
      if (k<NeN)
	{
	  // let l=k, m=0
	  l_index = k;
	  m_index = 0;
	}
      else if (k>=NeN)
	{
	  // let l=NeN-1, m=k-(NeN-1)
	  l_index = NeN-1;
	  m_index = k-l_index;
	}

      index_dep = l_index*NeN*Nee + m_index*Nee +1;

      paramDepMatrix(index_dep,index_dep) = 0.0;

      // Now we find the indexes of the parameters that determine the dependent
      // one
      if (k%2 == 0)
	{
	  index1 = k*NeN*Nee/2 + k*Nee/2 +1;

	  if (index_dep != index1) 
	    {      
	      if (l_index == m_index)
		paramDepMatrix(index_dep,index1) = -1.0;
	      else if (l_index != m_index)
		paramDepMatrix(index_dep,index1) = -0.5;
	    }
	}

      for (int l=0; l<NeN; l++)
	for (int m=0; m<l; m++)
	  if (l+m == k)
	    {
	      index1 = l*NeN*Nee + m*Nee + 1;

	      if (index_dep != index1)
		{
		  if (l_index == m_index)
		    paramDepMatrix(index_dep,index1) = -2.0;
		  else if (l_index != m_index)
		    paramDepMatrix(index_dep,index1) = -1.0;
		}
	    }
    }

  // Now we do the electron-nucleus cusp
  // The parameter 000 depends on 100 due to the e-n cusp (k=0)
  paramDepMatrix(0,0) = 0.0;
  paramDepMatrix(0,NeN*Nee) = cutoff/C;

  // The parameter 110 depends on 100 due to the e-n cusp (k=1)
  paramDepMatrix(NeN*Nee+Nee,NeN*Nee+Nee) = 0.0;
  paramDepMatrix(NeN*Nee+Nee,NeN*Nee) = C/cutoff;

  for (int k=2; k<NeN+Nee-1; k++)
    {
      // We find the index of the dependent parameter
      if (k<NeN)
	{
	  // let l=k, n=0
	  l_index = k;
	  n_index = 0;
	}
      else if (k==NeN)
	{
	  // let l=NeN-2, k=2
	  l_index = k-2;
	  n_index = 2;
	}
      else if (k>NeN)
	{
	  // let l=NeN-1, n=k-(NeN-1)
	  l_index = NeN-1;
	  n_index = k-l_index;
	}
      index_dep = l_index*NeN*Nee + Nee + n_index;

      paramDepMatrix(index_dep,index_dep) = 0.0;

      if (k<Nee)
	{
	  index1 = k;
	  index2 = NeN*Nee+k;
	  paramDepMatrix(index_dep,index1) = C/cutoff;
	  if (index_dep != index2)
	    paramDepMatrix(index_dep,index2) = -1.0;
	}

      for (int l=1; l<NeN; l++)
	for (int n=0; n<Nee; n++)
	  if (l+n == k)
	    {
	      index1 = l*NeN*Nee + n;
	      index2 = l*NeN*Nee + Nee + n;
	      paramDepMatrix(index_dep,index1) = C/cutoff;
	      if (index_dep != index2)
		paramDepMatrix(index_dep,index2) = -1.0;
	    }
    }

  // If the two body Jastrows are not to be reproduced, the terms proportional
  // to them will be set to zero
  if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 0)
    for (int l=1; l<NeN; l++)
      {
	index1 = l*NeN*Nee;
	paramDepMatrix(index1,index1) = 0.0;
      }

  if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 0)
    for (int n=2; n<Nee; n++)
      paramDepMatrix(n,n) = 0.0;

  // For every parameter that is not free, we substitute its definition in
  // terms of free parameters in the definitions of the parameters that depend
  // on it
  for (int i=0; i<TotalNumberOfParameters; i++)
    {
      if ( fabs(paramDepMatrix(i,i)-1.0 ) > 1e-10 )
	for (int j=0; j<TotalNumberOfParameters; j++)
	  {
	    for (int k=0; k<TotalNumberOfParameters; k++)
	      paramDepMatrix(j,k) += paramDepMatrix(j,i)*paramDepMatrix(i,k);
	    paramDepMatrix(j,i) = 0.0;
	  }
    }
	      
  // We can enforce symmetry by setting the elements with m>l equal to the ones
  // with l>m

  for (int l=0; l<NeN-1; l++)
    for (int m=l+1; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	{
	  index1 = l*NeN*Nee + m*Nee + n;
	  index2 = m*NeN*Nee + l*Nee + n;
	  for (int i=0; i<TotalNumberOfParameters; i++)
	    paramDepMatrix(index1,i) = paramDepMatrix(index2,i);
	}
}
