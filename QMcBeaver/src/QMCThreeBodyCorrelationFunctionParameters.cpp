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
  isFree  = rhs.isFree;

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
  int start = strm.tellg();   
  // Read the particle types
  
  string pt1, pt2, pt3;
  strm >> temp;

  if(temp != "ParticleTypes:")
    {
      strm.seekg(start);
      return false;
    }

  strm >> pt1 >> pt2 >> pt3;
  pt1 = StringManipulation::toFirstUpperRestLower(pt1);
  pt2 = StringManipulation::toFirstUpperRestLower(pt2);
  pt3 = StringManipulation::toFirstUpperRestLower(pt3);

  // Order the particle types
  
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
  
  strm >> temp;

  if(temp != "threeBodyCorrelationFunctionType:")
    {
      //it's not a 3 particle jastrow
      strm.seekg(start);
      return false;
    }

  strm >> threeBodyCorrelationFunctionType;
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

  Array1D<double> fileParameters(TotalNumberOfParameters);
  Parameters.allocate(TotalNumberOfParameters);
  isFree.allocate(TotalNumberOfParameters);
  strm >> temp;

  int pi = fileParameters.read(strm,0.0,"3 particle Jastrow");
  
  /**
     It is convenient to be able to sort parameters in the input
     file in such a way that increasing Nee or NeN will preserve optimized
     parameters. To do this, the input file orders parameters according to
     ascending powers.

     In the input file, the no coefficient for a term with a power of 2
     will appear ahead a term where the highest power is a 1. That is, the
     order will look like:

     terms with max power 0
     terms with max power 1
     terms with max power 2
     ...
     terms with max power max(NeN,Nee)

     This way, if you increase Nee or NeN, the coefficients for all the new terms
     will appear at the end of the list.

     If you degrease Nee or NeN, then the terms will be messed up.

     In summary, in order to keep the same coefficient for an independent (l,m,n) term
     when you increase Nee or NeN, you will not have to change the input file.
     The dependent parameters may change.
  */
  int counter = 0;
  for(int p=0; p<max(NeN,Nee); p++)
    for (int n=0; n<Nee; n++)
      for (int l=0; l<NeN; l++)
	for (int m=0; m<=l; m++)
	  if( (l == p || m == p || n == p) &&
	      (l <= p && m <= p && n <= p))
	    {
	      int index = map(l,m,n);
	      Parameters(index) = fileParameters(counter);
	      counter++;

	      //If there are symmetric elements, make sure they show up right next to each other.
	      if(l != m)
		{
		  index = map(m,l,n);
		  Parameters(index) = fileParameters(counter);
		  counter++;
		}
	    }
  fileParameters.deallocate();
  if(counter != TotalNumberOfParameters)
    {
      clog << "Error: TotalNumberOfParameters = " << TotalNumberOfParameters;
      clog << " but counter = " << counter << endl;
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

  //Count the number of free parameters
  NumberOfFreeParameters = 0;
  for(int p=0; p<isFree.dim1(); p++)
    if(isFree(p)) NumberOfFreeParameters++;
  
  //and add one for the cutoff parameter
  NumberOfFreeParameters++;
  
  // if the correlation function type is "None"
  // then remove any parameters or constants
  // that were listed
  
  if (threeBodyCorrelationFunctionType == "None")
    {
      NumberOfParameterTypes = 0;
      TotalNumberOfParameters = 0;
      NumberOfParameters.deallocate();
      Parameters.deallocate();
      isFree.deallocate();
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
  ParticleTypes.allocate(3);
  TotalNumberOfParameters = 0;
  NumberOfFreeParameters = 0;

  NeN = 0;
  Nee = 0;
  C = 0;
  cutoff = 0.0;

  setThreeBodyCorrelationFunction();
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

  /*
    When running correlated sampling, we don't need to run this every time.
    So how about we only check it when we are doing something else.
  */
  if(globalInput.cs_Parameters.dim1() == 0)
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

bool QMCThreeBodyCorrelationFunctionParameters::isUsed()
{
  if(threeBodyCorrelationFunctionType == "None")
    return false;
  return true;
}

QMCThreeBodyCorrelationFunctionParameters::~QMCThreeBodyCorrelationFunctionParameters()
{
  NumberOfParameters.deallocate();
  Parameters.deallocate();
  isFree.deallocate();
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

  pt_a.allocate(TotalNumberOfParameters+1);
  pt2_xa.allocate(TotalNumberOfParameters+1);
  pt3_xxa.allocate(TotalNumberOfParameters+1);

  for(int i=0; i<pt2_xa.dim1(); i++)
    pt2_xa(i).allocate(globalInput.WF.getNumberElectrons(),3); 
}

ostream & operator << (ostream & strm, 
		       QMCThreeBodyCorrelationFunctionParameters & rhs)
{
  int Nee = rhs.Nee;
  int NeN = rhs.NeN;

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

  int counter = 0;
  int endlInd = 0;
  for(int p=0; p<max(NeN,Nee); p++)
    {
      for (int n=0; n<Nee; n++)
	for (int l=0; l<NeN; l++)
	  for (int m=0; m<=l; m++)
	    if( (l == p || m == p || n == p) &&
		(l <= p && m <= p && n <= p))
	      {
		int index = rhs.map(l,m,n);
		if(endlInd % (NeN) == 0 && endlInd != 0)
		  {
		    strm << endl;
		    endlInd = 0;
		  }
		double val = rhs.Parameters(index);
		if(val < 0.0) strm << " -";
		else          strm << "  ";
		strm << setw(20) << left << fabs(val);	      
		
		counter++;
		endlInd++;

		if(l != m)
		  {
		    index = rhs.map(l,m,n);
		    if(endlInd % (NeN) == 0 && endlInd != 0)
		      {
			strm << endl;
			endlInd = 0;
		      }
		    val = rhs.Parameters(index);
		    if(val < 0.0) strm << " -";
		    else          strm << "  ";
		    strm << setw(20) << left << fabs(val);	      
		    
		    counter++;
		    endlInd++;
		  }
	      }
      strm << endl << endl;
      endlInd = 0;
    }

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
  p_a(shift) += pt_a.get(0);
  for(int e=0; e<globalInput.WF.getNumberElectrons(); e++)
    for(int xyz=0; xyz<3; xyz++)
      (p2_xa(shift))(e,xyz) += (pt2_xa.get(0))(e,xyz);  
  p3_xxa(shift) += pt3_xxa.get(0);
  int counter = 1;

  /*
    The derivatives were calculated wrt the total parameters.
    However, we will only optimize wrt the free parameters, so
    we need to translate them.
  */
  for (int l=0; l<NeN; l++)
    for (int m=0; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	if(isFree.get(map(l,m,n)))
	  {
	    for (int i=0; i<TotalNumberOfParameters; i++)
	      {
		double pd = paramDepMatrix.get(i,map(l,m,n));
		if( fabs(pd) < 1.0e-50) continue;
		
		p_a(shift+counter)    += pd * pt_a.get(i+1);
		p3_xxa(shift+counter) += pd * pt3_xxa.get(i+1);
		
		for(int e=0; e<globalInput.WF.getNumberElectrons(); e++)
		  for(int xyz=0; xyz<3; xyz++)
		    (p2_xa(shift+counter))(e,xyz) += pd * (pt2_xa.get(i+1))(e,xyz);
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
  free(0) = cutoff;
  counter++;

  // Because of the symmetry due to interchange of the electrons, we only 
  // write the parameters that have l>=m.
  for (int l=0; l<NeN; l++)
    for (int m=0; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	if(isFree.get(map(l,m,n)))
	  {
	    free(counter) = total.get(map(l,m,n));
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
  cutoff = free.get(0);
  counter++;

  // We only consider the elements such that l>=m.  The others will be filled
  // in by symmetry
  for (int l=0; l<NeN; l++)
    for (int m=0; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	if(isFree(map(l,m,n)))
	  {
	    total(map(l,m,n)) = free.get(counter);
	    counter++;
	  }
}

bool QMCThreeBodyCorrelationFunctionParameters::checkCuspAndSymmetry()
{
  double tolerance = 1e-10;
  int numFree = 0;
  for (int l=0; l<NeN; l++)
    for (int m=0; m<NeN; m++)
      for (int n=0; n<Nee; n++)       
	{
	  int c = map(l,m,n);
	  int label = 100*l + 10*m + n;
	  //the variable is independent if the row is identity
	  Array2D<double> tranPDM = paramDepMatrix;
	  tranPDM.transpose();
	  bool isInd = c == tranPDM.isDependent(c);
	  if(isInd) numFree++;
	  if(isInd != isFree(c))
	    {
	      cout << "Error for variable " << label << ": isFree = "
		   << isFree(c) << " while paramDepMatrix says independent." << endl;
	      cout << "paramDepMatrix is flawed..." << endl;
	      return false;
	    }
	  //the variable is dependent if the column is all zero
	  double sum = 0;
	  int dependsOnDependent = -1;
	  double dependVal = 0;
	  for(int r=0; r<paramDepMatrix.dim1(); r++)
	    for (int ll=0; ll<NeN; ll++)
	      for (int mm=0; mm<NeN; mm++)
		for (int nn=0; nn<Nee; nn++)       		  
		  {
		    int r = map(ll,mm,nn);
		    int labelr = 100*ll + 10*mm + nn;
		    if(!isFree(r) && fabs(paramDepMatrix(r,c)) > tolerance)
		      {
			dependsOnDependent = labelr;
			dependVal = paramDepMatrix(r,c);
		      }
		    sum += fabs(paramDepMatrix(r,c));
		  }
	  
	  if(!isFree(c) && dependsOnDependent != -1)
	    {
	      cout << "Error: dependent parameter " << label
		   << " depends on dependent parameter " << dependsOnDependent
		   << " with coefficient " << dependVal << endl;
	    }
	  
	  if(sum < tolerance && isFree(c))
	    {
	      cout << "Error for variable " << label << ": isFree = "
		   << isFree(c) << " while paramDepMatrix says dependent." << endl;
	      cout << "paramDepMatrix is flawed..." << endl;
	      return false;
	    }
	}

  if(NumberOfFreeParameters - 1 != numFree)
    {
      cout << "Error: incorrect number of free parameters. We counted " << numFree;
      //we're not including cutoff in the counting here
      cout << " but there should have been " << (NumberOfFreeParameters-1) << endl;
    }

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
	  index1 = map(l,m,n);
	  index2 = map(m,l,n);
	  double diff = Parameters(index1)-Parameters(index2);

	  //make it a relative measurement
	  if(fabs(Parameters(index1)) > tolerance) diff /= Parameters(index1);

	  if ( fabs(diff) > tolerance )
	    {
	      clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
		   << "checkCuspAndSymmetry(): symmetry is not satisfied, difference = " << diff << endl;
	      clog << "   Parameters(map(" << l << "," << m << "," << n << ")) = " << setw(20) << Parameters(index1) << endl;
	      clog << "   Parameters(map(" << m << "," << l << "," << n << ")) = " << setw(20) << Parameters(index2) << endl;
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
	      index1 = map(l,m,1);
	      sum += Parameters(index1);
	    }
      if (fabs(sum) > tolerance)
	{
	  clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
	       << "checkCuspAndSymmetry():" << endl
	       << "(l+m=" << k << ") electron-electron cusp condition "
	       << "is not satisfied, sum = " << sum << endl;
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
	      index1 = map(0,m,n);
	      index2 = map(1,m,n);
	      sum += C*Parameters(index1) - Parameters(index2);
	    }
      if (fabs(sum) > tolerance)
	{
	  clog << "ERROR in QMCThreeBodyCorrelationFunctionParameters::"
	       << "checkCuspAndSymmetry():" << endl
	       << "(m+n=" << k << ") electron-nucleus cusp condition "
	       << "is not satisfied, sum = " << sum << endl;
	  return false;
	}
    }
  return true;
}

void QMCThreeBodyCorrelationFunctionParameters::gaussianParamDepMatrix()
{
  /*
    There is a difference between the definition of Nee and NeN as used
    in this code, and as used in the PRB 2004 paper. In order to make the
    connection clear, I'm going to use the definitions from the paper.

    That is, NeN is the number of terms, nn is the highest power used.
  */
  int nn = NeN-1;
  int ne = Nee-1;

  double tolerance = 1e-50;
  int numConstraints = 0;
  int numVariables = Nee * NeN * Nee;
  numConstraints += Nee*NeN*(NeN-1)/2; 
  numConstraints += 2*nn + 1;
  numConstraints += nn + ne + 1;

  if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 0)
    numConstraints += nn;

  if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 0)
    numConstraints += ne;

  Array2D<double> constraints(numConstraints,numVariables);
  Array1D<int> labels(numVariables);
  Array1D<int> numbers(numVariables);

  constraints = 0.0;
  for (int l=0; l<NeN; l++)
    for (int m=0; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	{
	  int index = map(l,m,n);
	  labels(index) = 100*l + 10*m + n;
	  numbers(index) = index;

	  if(n == 1)
	    {
	      //this is equation 24
	      int k = l + m;
	      constraints(k,index) = 1.0;
	    }
	  int k0 = 2*nn+1;
	  if(l == 0)
	    {
	      //this is equation 27, for ri cusp
	      int k = m + n + k0;
	      constraints(k,index) = C;
	    }
	  if(l == 1)
	    {
	      //this is equation 27, for ri cusp
	      int k = m + n + k0;
	      constraints(k,index) = -1.0;
	    }
	  k0 += nn+ne+1;
	  int num = NeN * (NeN - 1) / 2;
	  if(l > m)
	    {
	      //impose symmetry between l and m
	      int k = l*(l-1)/2 + m + n*num + k0;
	      constraints(k,index) = -1.0;
	    }
	  if(l < m)
	    {
	      //impose symmetry between l and m
	      int k = m*(m-1)/2 + l + n*num + k0;
	      constraints(k,index) = 1.0;
	    }
	  k0 += Nee * num;
	  if(globalInput.flags.reproduce_EE_with_NEE_jastrow == 0)
	    {
	      int k = n - 1 + k0;
	      if(l == 0 && m == 0 && n > 0)
		{
		  constraints(k,index) = 1.0;
		}
	      k0 += ne;
	    }
	  if(globalInput.flags.reproduce_NE_with_NEE_jastrow == 0)
	    {
	      int k = l - 1 + k0;
	      if(l > 0 && m == 0 && n == 0)
		{
		  constraints(k,index) = 1.0;
		}
	      k0 += nn;
	    }
	}

  /*
    These are all the variables that we know will be zero, so solve for them
    right away. There is some overlap already.
  */
  if (globalInput.flags.reproduce_NE_with_NEE_jastrow == 0)
    for(int l=1; l<NeN; l++)
      constraints.rref(map(l,0,0));

  if (globalInput.flags.reproduce_EE_with_NEE_jastrow == 0)
    for(int n=1; n<Nee; n++)
      constraints.rref(map(0,0,n));

  constraints.rref(map(0,0,1));
  constraints.rref(map(1,0,1));
  constraints.rref(map(nn,nn,1));
  constraints.rref(map(nn,nn-1,1));

  /*
    Now solve for all the parameters that we know are
    symmetric.
  */
  for (int n=0; n<Nee; n++)
    for (int l=0; l<NeN; l++)
      for (int m=0; m<NeN; m++)
	{
	  int index = map(l,m,n);
	  if(l < m) constraints.rref(index);
	}

  /*
    At this point, there are 3 nn + ne - 2 more
    parameters to make dependent.

    We can never choose terms lln (l > 1, n != 1) to be dependent.
  */
  if(NeN == 3 && Nee == 3)
    {
      //my goal here is to eliminate as many terms with zeros as i can
      constraints.rref(map(1,0,0));
      constraints.rref(map(2,0,0));
      constraints.rref(map(2,0,2));
      constraints.rref(map(2,0,1));
      constraints.rref(map(1,1,0));
      constraints.rref(map(1,0,2));
      constraints.rref(map(0,1,2));
    } else if(NeN == 4 && Nee == 4)
      {
	//these three are connected to 0,0,0
	//so if we solve for these, 0,0,0 is guaranteed
	//to be independent
	constraints.rref(map(1,1,0));
	constraints.rref(map(1,0,0));
	constraints.rref(map(0,1,0));

	constraints.rref(map(2,0,0));
	constraints.rref(map(3,0,0));	
	constraints.rref(map(1,1,0));
	constraints.rref(map(3,0,1));
	constraints.rref(map(2,0,1));
	constraints.rref(map(3,0,2));
	constraints.rref(map(3,0,3));
	constraints.rref(map(2,0,2));
      }

  //this helps identify constraints that are available for solving for variables
  for(int con=0; con<constraints.dim1(); con++)
    {
      if(constraints.constraintUsed(con) == -1)
	printf("constraint %3i is available\n",con);
    }

  //negative zero is annoying in printouts
  for(int i=0; i<constraints.dim1(); i++)
    for(int j=0; j<constraints.dim2(); j++)
      if(fabs(constraints(i,j)) < 1e-100) constraints(i,j) = 0.0;

  int cs = 3;
  if(false)
    {
      for(int i=0; i<constraints.dim2(); i++)
	if(constraints.isDependent(i) == -1)
	  printf("var L %3i N %3i is independent, has value = %20.10g\n",labels(i),i,Parameters(i));

      cout << " ";
      numbers.printArray1D(cout,2,cs,-1,' ',false);
      cout << endl << " ";
      labels.printArray1D(cout,2,cs,-1,' ',false);
      cout << endl;
      constraints.printArray2D(cout,2,cs,-1,' ',false);
      cout << endl;
    }

  /*
    Now convert these constraints into our parameter dependency
    matrix. This matrix is defined by:
    Look across the row, you'll see what this parameter depends on
    or
    Look down the column to see who uses this parameter
    
    If the parameter is independent, it will have a 1 on the diagonal, and
    0 for the rest of the row.
    If the parameter is dependent, it will have a 0 on the entire column.

    Notice: for terms which have a symmetric counter part, we expect the 
    l < m term to be dependent. Loops were programmed to scan only l > m for
    independent parameters.
  */
  paramDepMatrix.allocate(TotalNumberOfParameters,TotalNumberOfParameters);
  paramDepMatrix.setToIdentity();

  for(int j=0; j<constraints.dim2(); j++)
    {
      //If it is a dependent variable, then it will tell us the row
      //with the dependencies
      int theOne = constraints.isDependent(j);
      isFree(j) = theOne == -1 ? true : false;

      if(!isFree(j))
	{
	  /*
	    The dependency multiplied by -1 because the all the variables
	    are on one side of the equation,
	    depdendent variable + independent variables = 0
	    so by negating it, we are solving
	    dependent variable = - indepdendent variables
	  */
	  for(int d=0; d<constraints.dim2(); d++)
	    if(fabs(constraints(theOne,d)) > tolerance)
	      paramDepMatrix(j,d) = -1.0 * constraints(theOne,d);
	  paramDepMatrix(j,j) = 0.0;
	}
    }

  if(constraints.numDependent() != constraints.dim1())
    {
      cout << "Error: you didn't use up all the constraints! There are too few dependent variables." << endl;
      cout << " number   dependent = " << constraints.numDependent() << endl;
      cout << " number constraints = " << constraints.dim1() << endl;
    }

  if(false)
    {
      cs = 3;
      cout << " ";
      numbers.printArray1D(cout,2,cs,-1,' ',false);
      cout << endl << " ";
      isFree.printArray1D(cout,2,cs,-1,' ',false);
      cout << endl << " ";
      labels.printArray1D(cout,2,cs,-1,' ',false);
      cout << endl;
      paramDepMatrix.printArray2D(cout,2,cs,-1,' ',false);
      cout << endl;
    }
}

void QMCThreeBodyCorrelationFunctionParameters::makeParamDepMatrix()
{
  gaussianParamDepMatrix();
  return;

  for (int l=0; l<NeN; l++)
    for (int m=0; m<NeN; m++)
      for (int n=0; n<Nee; n++)
	{
	  bool free = true;

	  if (m>l)
	    {
	      // These parameters are dependent by symmetry
	      free = false;
	    }
	  else if (l>0 && m==0 && n==0)
	    {
	      if (globalInput.flags.reproduce_NE_with_NEE_jastrow != 1)
		free = false;
	    }
	  else if (l==0 && n>1)
	    {
	      if (globalInput.flags.reproduce_EE_with_NEE_jastrow != 1)
		free = false;
	    }
	  // These parameters are dependent due to the EE no-cusp condition
	  else if (m==0 && n==1)
	    {
	      // Do nothing- these parameters are dependent
	      free = false;
	    }
	  else if (l==NeN-1 && n==1)
	    {
	      // Do nothing- these parameters are dependent
	      free = false;
	    }
	  // These parameters are dependent due to the EN no-cusp condition
	  else if (l==0 && n==0)
	    {
	      // Do nothing- this parameter is dependent
	      free = false;
	    }
	  else if (l<NeN && m==1 && n==0)
	    {
	      // Do nothing- these parameters are dependent
	      free = false;
	    }
	  else if (l==NeN-2 && m==1 && n==2)
	    {
	      // Do nothing- this parameter is dependent
	      free = false;
	    }
	  else if (l==NeN-1 && m==1 && n>1)
	    {
	      // Do nothing- these parameters are dependent
	      free = false;
	    }

	  isFree(map(l,m,n)) = free;
	}

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

  //if the cutoff term looks like (r L - 1.0)
  double shift = 1.0;

  //if the cutoff term looks like (r - L)
  //double shift = cutoff;

  // Now we do the electron-nucleus cusp
  // The parameter 000 depends on 100 due to the e-n cusp (k=0)
  paramDepMatrix(0,0) = 0.0;
  paramDepMatrix(0,NeN*Nee) = shift/C;

  // The parameter 110 depends on 100 due to the e-n cusp (k=1)
  paramDepMatrix(NeN*Nee+Nee,NeN*Nee+Nee) = 0.0;
  paramDepMatrix(NeN*Nee+Nee,NeN*Nee) = C/shift;

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
	  paramDepMatrix(index_dep,index1) = C/shift;
	  if (index_dep != index2)
	    paramDepMatrix(index_dep,index2) = -1.0;
	}

      for (int l=1; l<NeN; l++)
	for (int n=0; n<Nee; n++)
	  if (l+n == k)
	    {
	      index1 = l*NeN*Nee + n;
	      index2 = l*NeN*Nee + Nee + n;
	      paramDepMatrix(index_dep,index1) = C/shift;
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
              
              if(isFree(map(l,m,n)))
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
              else if (isFree(map(l,m,n)) == false)
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
                          if (isFree(map(a,b,c)) == false && 
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
