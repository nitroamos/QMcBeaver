#include "QMCThreeBodyCorrelationFunctionParameters.h"
#include "QMCInput.h"
#include <iomanip>

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
	  clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
	  clog << "particles " << pt1 << ", " << pt2 << ", " << pt3 << " are ";
	  clog << "being loaded as a three body correlation function!" << endl;
	  ok = false;
	}
    }
  else if ( pt2 == "Electron_down")
    {
      if (pt3 != "Electron_down")
	{ 
	  // If we have an opposite spin pair, we assume the up electron is 
	  // listed first.
	  clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
	  clog << "particles " << pt1 << ", " << pt2 << ", " << pt3 << " are ";
	  clog << "being loaded as a three body correlation function!" << endl;
	  ok = false;
	}
    }
  
  // Load the correlation function type
  
  strm >> temp >> threeBodyCorrelationFunctionType;
  if (threeBodyCorrelationFunctionType != "Cambridge" && 
      threeBodyCorrelationFunctionType != "None")
    {
      clog << "ERROR: unknown type of ThreeBodyCorrelationFunction is being ";
      clog << "read.  threeBodyCorrelationFunctionType = ";
      clog << threeBodyCorrelationFunctionType << endl;
      ok = false;
    }

  // Load the number of parameter types
  strm >> temp;
  strm >> temp;
  NumberOfParameterTypes = atoi(temp.c_str());

  if (NumberOfParameterTypes != 2)
    {
      // There should be two parameter types: NeN and Nee.
      clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
      clog << "NumberOfParametersTypes = " << NumberOfParameterTypes;
      clog << ", 2 was expected." << endl;
      ok = false;
      return ok;
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
      clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters: ";
      clog << "NeN = " << NeN << ", Nee = " << Nee << ", values of at least ";
      clog << "3 are needed" << endl;
      ok = false;
    }

  // Set the total number of parameters

  TotalNumberOfParameters = NeN*NeN*Nee;

  // Load the parameters
  
  Parameters.allocate(TotalNumberOfParameters);
  strm >> temp;

  int pi = Parameters.read(strm,0.0,"3 particle Jastrow");
  
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

  //cout << "param dep matrix:" << endl;  
  //paramDepMatrix.printArray2D(cout,0,0,-1,' ',false);
  //paramDepMatrix.printArray2D(cout,2,2,-1,' ',false);

  // Now we calculate how many of the parameters are free to be optimized.  
  // Some of the parameters are constrained by symmetry or cusp conditions.

  // If everything was free to be optmized, we would have all of the parameters
  // and the cutoff distance

  NumberOfFreeParameters = TotalNumberOfParameters+1; 

  if (globalInput.flags.optimize_L == 0)
    {
      NumberOfFreeParameters -= 1;
    }
  else
    {
      /*
	until i can figure out a way to reliably optimize 2 particle cutoffs, i don't
	see the point here.
      */
      cerr << "Error: sorry... I didn't program NEE cutoff derivatives! :-)" << endl;
      exit(0);
    }

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
      setTotalParameters(Parameters);
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

int QMCThreeBodyCorrelationFunctionParameters::getNumberOfTotalParameters()
{
  return TotalNumberOfParameters;
}

void QMCThreeBodyCorrelationFunctionParameters::zeroOutDerivatives()
{
  pt_a = 0.0;
  pt3_xxa = 0.0;
  for(int i=0; i<pt2_xa.dim1(); i++)
    pt2_xa(i) = 0.0;
}

void QMCThreeBodyCorrelationFunctionParameters::getFree(const Array1D<double> & total,
							Array1D<double> & free) 
{
  tempArray.allocate(TotalNumberOfParameters);
  tempArray = 0.0;
  
  //make apply the transformation to the array in the total basis set
  for (int i=0; i<TotalNumberOfParameters; i++)
    for (int j=0; j<TotalNumberOfParameters; j++)
      tempArray(i) += paramDepMatrix(i,j)*total.get(j);
  
  totalToFree(tempArray, free);  
}

void QMCThreeBodyCorrelationFunctionParameters::setTotalParameters(const Array1D<double> & total)
{
  tempArray.allocate(TotalNumberOfParameters);
  tempArray = 0.0;
  
  // We multiply the paramDepMatrix by the parameters we loaded in to make
  // sure the parameters satisfy the symmetry and cusp conditions
  for (int i=0; i<TotalNumberOfParameters; i++)
    for (int j=0; j<TotalNumberOfParameters; j++)
      tempArray(i) += paramDepMatrix(i,j)*total.get(j);
  
  Parameters = tempArray;
  
  if(!checkCuspAndSymmetry())
    {
      clog << "Error: checkCuspAndSymmetry failed in setTotalParameters" << endl;
      exit(0);
    }

  getFree(Parameters, freeParameters);
}

void QMCThreeBodyCorrelationFunctionParameters::setFreeParameters(const Array1D<double> & free)
{
  if( free.dim1() != NumberOfFreeParameters )
    {
      clog << "ERROR: Parameters of the incorrect size are trying to be set "
	   << "in setFreeParameters" << endl;
    }

  freeParameters = free;
  freeToTotal(freeParameters, Parameters);

  tempArray.allocate(TotalNumberOfParameters);
  tempArray = 0.0;
  
  for (int i=0; i<TotalNumberOfParameters; i++)
    for (int j=0; j<TotalNumberOfParameters; j++)
      tempArray(i) += paramDepMatrix(i,j)*Parameters(j);
  
  Parameters = tempArray;

  if(!checkCuspAndSymmetry())
    {
      clog << "Error: checkCuspAndSymmetry failed in setFreeParameters" << endl;
      exit(0);
    }
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

  pt_a.allocate(TotalNumberOfParameters);
  pt2_xa.allocate(TotalNumberOfParameters);
  pt3_xxa.allocate(TotalNumberOfParameters);

  for(int i=0; i<pt2_xa.dim1(); i++)
    pt2_xa(i).allocate(globalInput.WF.getNumberElectrons(),3); 
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

  strm << "Parameters:" << endl;

  for (int i = 0; i < rhs.TotalNumberOfParameters; i++)
    {
      if(i % (rhs.NeN) == 0 && i != 0)
	strm << endl;
      double val = rhs.Parameters(i);
      if(val < 0.0) strm << " -";
      else          strm << "  ";
      strm << setw(20) << left << fabs(val);
    }
  strm << endl;

  strm << "C:  " << rhs.C << endl;
  
  strm << "Cutoff:  " << rhs.cutoff << endl;

  strm << endl;

  return strm;
}


void QMCThreeBodyCorrelationFunctionParameters::totalDerivativesToFree(int shift,
								       Array1D<double> & p_a,
								       Array1D< Array2D<double> > & p2_xa,
								       Array1D<double> & p3_xxa) const
{
  int counter = 0;

  // If the cutoff distance is being optimized it is the 0th element of the
  // free parameter array
  if (globalInput.flags.optimize_L == 1)
    {
      p_a(0) = pt_a.get(0);
      counter++;
    }

  /*
    The derivatives were calculated wrt the total parameters.
    However, we will only optimize wrt the free parameters, so
    we need to translate them.
  */
  for (int l=0; l<NeN; l++)
    for (int m=0; m<=l; m++)
      for (int n=0; n<Nee; n++)
	if(isFree(l,m,n))
	  {
	    int index = l*NeN*Nee+m*Nee+n;
	    for (int i=0; i<TotalNumberOfParameters; i++)
	      {
		/*
		  pd is the coefficient in the linear combination
		  of the dependence of this free parameter on
		  this total parameter.
		*/
		double pd = paramDepMatrix.get(i,index);
		if( fabs(pd) < 1.0e-50) continue;

		p_a(shift+counter)    += pd * pt_a.get(i);
		p3_xxa(shift+counter) += pd * pt3_xxa.get(i);

		for(int e=0; e<globalInput.WF.getNumberElectrons(); e++)
		  for(int xyz=0; xyz<3; xyz++)
		    (p2_xa(shift+counter))(e,xyz) += pd * (pt2_xa.get(i))(e,xyz);
	      }
	    counter++;
	  }

  if (counter != NumberOfFreeParameters)
    {
      // We have miscalculated the number of free parameters.
      clog << "ERROR in "
	   << "QMCThreeBodyCorrelationFunctionParameters::totalDerivativesToFree()"
	   << ": there should be " << NumberOfFreeParameters << " free "
	   << "parameters, but we found " << counter << endl;
      exit(0);
    }
}

void QMCThreeBodyCorrelationFunctionParameters::totalToFree(const Array1D<double> & total,
							    Array1D<double> & free) const
{
  // This function takes the full set of parameters and makes the set of 
  // parameters that are free to be optimized.  The decision is based on the 
  // flags of what we are optimizing and the symmetry of the function.

  free.allocate(NumberOfFreeParameters);
  free = 0.0;
  int counter = 0;

  // If the cutoff distance is being optimized it is the 0th element of the
  // free parameter array
  if (globalInput.flags.optimize_L == 1)
    {
      free(0) = cutoff;
      counter++;
    }

  // Because of the symmetry due to interchange of the electrons, we only 
  // write the parameters that have l>=m.
  for (int l=0; l<NeN; l++)
    for (int m=0; m<=l; m++)
      for (int n=0; n<Nee; n++)
	if(isFree(l,m,n))
	  {
	    free(counter) = total.get(l*NeN*Nee+m*Nee+n);
	    counter++;
	  }

  if (counter != NumberOfFreeParameters)
    {
      // We have miscalculated the number of free parameters.
      clog << "ERROR in "
	   << "QMCThreeBodyCorrelationFunctionParameters::totalToFree()"
	   << ": there should be " << NumberOfFreeParameters << " free "
	   << "parameters, but we found " << counter << endl;
      exit(0);
    }
}

void QMCThreeBodyCorrelationFunctionParameters::freeToTotal(const Array1D<double> & free,
							    Array1D<double> & total)
{
  // This function takes the set of free parameters and generates the set of 
  // all parameters to send to the correlation function.

  int counter = 0;
  total = 0.0;

  // If the cutoff distance is being optimized it is the 0th element of the
  // free parameter array
  if (globalInput.flags.optimize_L == 1)
    {
      cutoff = free.get(0);
      counter++;
    }

  // We only consider the elements such that l>=m.  The others will be filled
  // in by symmetry
  for (int l=0; l<NeN; l++)
    for (int m=0; m<=l; m++)
      for (int n=0; n<Nee; n++)
	if(isFree(l,m,n))
	  {
	    total(l*NeN*Nee+m*Nee+n) = free.get(counter);
	    counter++;
	  }
}

inline bool QMCThreeBodyCorrelationFunctionParameters::isFree(int l, int m, int n) const
{
  if (m>l)
    // These parameters are dependent by symmetry
    return false;

  if (l>0 && m==0 && n==0)
    {
      if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 1)
	return true;
      return false;
    }
  else if (l==0 && n>1)
    {
      if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 1)
	return true;
      return false;
    }
  // These parameters are dependent due to the EE no-cusp condition
  else if (m==0 && n==1)
    {
      // Do nothing- these parameters are dependent
      return false;
    }
  else if (l==NeN-1 && n==1)
    {
      // Do nothing- these parameters are dependent
      return false;
    }
  // These parameters are dependent due to the EN no-cusp condition
  else if (l==0 && n==0)
    {
      // Do nothing- this parameter is dependent
      return false;
    }
  else if (l<NeN && m==1 && n==0)
    {
      // Do nothing- these parameters are dependent
      return false;
    }
  else if (l==NeN-2 && m==1 && n==2)
    {
      // Do nothing- this parameter is dependent
      return false;
    }
  else if (l==NeN-1 && m==1 && n>1)
    {
      // Do nothing- these parameters are dependent
      return false;
    }

  return true;
}

bool QMCThreeBodyCorrelationFunctionParameters::checkCuspAndSymmetry()
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
	      clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
		   << "checkCuspAndSymmetry(): symmetry is not satisfied"
		   << endl;
	      return false;
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
	  clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
	       << "checkCuspAndSymmetry(): electron-electron cusp condition "
	       << "is not satisfied" << endl;
	  return false;
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
	  clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
	       << "checkCuspAndSymmetry(): electron-nucleus cusp condition "
	       << "is not satisfied" << endl;
	  return false;
	}
    }
  return true;
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

  int l_index   = -1;
  int m_index   = -1;
  int n_index   = -1;
  int index1    = -1;
  int index2    = -1;
  int index_dep = -1;

  // We know that all parameters with m>l are dependent by symmetry
  // These rows will get filled in at the end of this function
  for (int l=0; l<NeN-1; l++)
    for (int m=l+1; m<NeN; m++)
      for (int n=0; n<Nee; n++)
        {
          index1 = l*NeN*Nee + m*Nee +n;
          paramDepMatrix(index1,index1) = 0.0;
        }

  // We start with the electron-electron cusp

  // The parameter 001 is zero due to the electron-electron cusp (k=0)
  paramDepMatrix(1,1) = 0.0;

  // The parameter 101 is zero due to the electron-electron cusp (k=1)
  paramDepMatrix(NeN*Nee+1,NeN*Nee+1) = 0.0;

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

  // Here we use a while loop to make sure that all dependent parameters are 
  // defined in terms of independent parameters.  Dependent params should not 
  // be in the definition of any parameter.

  bool goodDefs = false;
  int counter = 0;
  
  while (goodDefs == false)
    {
      counter++;
      if (counter > 100)
        {
          cerr << "ERROR in QMCThreeBodyCorrelationFunctionParameters::make"
               << "ParamDepMatrix(): dependent parameters are still used in "
               << "definitions after 100 iterations- definitions may be "
               << "circular" << endl;
          exit(0);
        }

      goodDefs = true;
      for (int l=0; l<NeN; l++)
        for (int m=0; m<=l; m++)
          for (int n=0; n<Nee; n++)
            {
              index1 = l*NeN*Nee + m*Nee + n;
              
              if (isFree(l,m,n) == true)
                {
                  // If this parameter is free, its diagonal element in the 
                  // matrix should be 1.0 and it should have nothing else in
                  // its row.

                  if ( fabs(paramDepMatrix(index1,index1)-1.0) > 1e-10 )
                    {
                      cerr << "ERROR in QMCThreeBodyCorrelationFunction"
                           << "Parameters::makeParamDepMatrix(): parameter "
                           << l << m << n << " should be free, but its "
                           << "diagonal matrix element is "
                           << paramDepMatrix(index1,index1) << endl;
                      exit(0);
                    }

                  for (int i=0; i<TotalNumberOfParameters; i++)
                    if (index1 != i && fabs(paramDepMatrix(index1,i)) > 1e-10)
                      {
                        cerr << "ERROR in QMCThreeBodyCorrelationFunction"
                             << "Parameters::makeParamDepMatrix(): parameter "
                             << l << m << n << " should be free, but element "
                             << "(" << index1 << "," << i << ") is "
                             << paramDepMatrix(index1,i) << endl;
                        exit(0);
                      }
                }
              else if (isFree(l,m,n) == false)
                {
                  // If this parameter is dependent, its entire column should
                  // be zeros.  Any nonzero elements in its row should
                  // correspond to free parameters.
                  if ( fabs(paramDepMatrix(index1,index1)) > 1e-10 )
                    {
                      cerr << "ERROR in QMCThreeBodyCorrelationFunction"
                           << "Parameters::makeParamDepMatrix(): parameter "
                           << l << m << n << "should be dependent, but "
                           << "its diagonal element is "
                           << paramDepMatrix(index1,index1) << endl;
                      exit(0);
                    }
                  // Now we check the column
                  for (int i=0; i<TotalNumberOfParameters; i++)             
                    if (index1 != i && fabs(paramDepMatrix(i,index1)) > 1e-10)
                      {
                        goodDefs = false;
                        for (int j=0; j<TotalNumberOfParameters; j++)
                          paramDepMatrix(i,j) += paramDepMatrix(i,index1) *
                            paramDepMatrix(index1,j);
                        paramDepMatrix(i,index1) = 0.0;                     
                      }
                  // Finally we check the row
                  for (int a=0; a<NeN; a++)
                    for (int b=0; b<NeN; b++)
                      for (int c=0; c<Nee; c++)
                        {
                          index2 = a*NeN*Nee + b*Nee + c;
                          if (isFree(a,b,c) == false && 
                              fabs(paramDepMatrix(index1,index2)) > 1e-10)
                            {
                              goodDefs = false;
                              for (int i=0; i<TotalNumberOfParameters; i++)
                                paramDepMatrix(index1,i) += 
                                  paramDepMatrix(index1,index2) * 
                                  paramDepMatrix(index2,i);
                              paramDepMatrix(index1,index2) = 0.0;
                            }
                        }
                }
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
