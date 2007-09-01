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

#include "QMCCorrelationFunctionParameters.h"

Array1D < double > QMCCorrelationFunctionParameters::getParameters()
{
  return Parameters;
}

Array1D < Complex > QMCCorrelationFunctionParameters::getPoles()
{
  return CorrelationFunction->getPoles();
}

void QMCCorrelationFunctionParameters::operator = (const 
				 QMCCorrelationFunctionParameters & rhs)
{
  NumberOfParameterTypes = rhs.NumberOfParameterTypes;
  NumberOfParameters = rhs.NumberOfParameters;
  Parameters = rhs.Parameters;
  BeginningIndexOfParameterType = rhs.BeginningIndexOfParameterType;
  TotalNumberOfParameters = rhs.TotalNumberOfParameters;

  NumberOfConstantTypes = rhs.NumberOfConstantTypes;
  NumberOfConstants = rhs.NumberOfConstants;
  Constants = rhs.Constants;
  BeginningIndexOfConstantType = rhs.BeginningIndexOfConstantType;
  TotalNumberOfConstants = rhs.TotalNumberOfConstants;

  ParticleTypes = rhs.ParticleTypes;
  CorrelationFunctionType = rhs.CorrelationFunctionType;

  setCorrelationFunction();
  initializeCorrelationFunctionParameters();
}

/**
 *Load QMCCorrelationFunctionParameters of the form: <br>
 * ParticleTypes: He Electron_up <br>
 * CorrelationFunctionType: Pade
 * NumberOfParameterTypes: 3 <br>
 * NumberOfParametersOfEachType: 4 2 1 <br>
 * Parameters: 1.0 2.0 3.0 4.0 1.0 2.0 1.0 <br>
 * NumberOfConstantTypes: 3 <br>
 * NumberOfConstantsOfEachType: 4 2 1 <br>
 * Constants: 1.0 2.0 3.0 4.0 1.0 2.0 1.0
 */

bool QMCCorrelationFunctionParameters::read(istream & strm, bool nucCuspReplacement)
{
  bool ok = true;
  string temp;
  
  // Read the particle types
  
  string pt1, pt2;
  strm >> temp >> pt1 >> pt2;
  pt1 = StringManipulation::toFirstUpperRestLower(pt1);
  pt2 = StringManipulation::toFirstUpperRestLower(pt2);

  // Order the particle types

  ParticleTypes.allocate(2);
  
  if( pt1 == "Electron_up" )
    {
      ParticleTypes(0) = pt1;
      ParticleTypes(1) = pt2;
    }
  else if( pt2 == "Electron_up" )
    {
      ParticleTypes(0) = pt2;
      ParticleTypes(1) = pt1;
    }
  else if( pt1 == "Electron_down" )
    {
      ParticleTypes(0) = pt1;
      ParticleTypes(1) = pt2;
    }
  else if( pt2 == "Electron_down" )
    {
      ParticleTypes(0) = pt2;
      ParticleTypes(1) = pt1;
    }
  else
    {
      cerr << "ERROR: Two non-electron particles (" << pt1 << ", " << pt2 << ") in one correlation function!"
	   << endl;
      ok = false;
    }
  
  // Load the correlation function type
  
  strm >> temp >> CorrelationFunctionType;

  // Load the number of parameter types
  
  strm >> temp >> temp;
  NumberOfParameterTypes = atoi(temp.c_str());
  
  // Load the number of paramters of each type
  
  NumberOfParameters.allocate(NumberOfParameterTypes);
  strm >> temp;
  
  for (int i = 0; i < NumberOfParameterTypes; i++)
    {
      strm >> temp;
      NumberOfParameters(i) = atoi(temp.c_str());
    }
  
  // Find the beginning index of each parameter type
  
  BeginningIndexOfParameterType.allocate(NumberOfParameterTypes);
  
  for (int i = 0; i < NumberOfParameterTypes; i++)
    {
      if (i == 0)
 	{
 	  BeginningIndexOfParameterType(0) = 0;
 	}
      else
 	{
 	  BeginningIndexOfParameterType(i) =
 	    BeginningIndexOfParameterType(i - 1) + 
 	    NumberOfParameters(i-1);
 	}
    }
  
  // Set the total number of parameters
  
  if( NumberOfParameterTypes > 0 )
    {
      TotalNumberOfParameters = BeginningIndexOfParameterType(
        NumberOfParameterTypes - 1) +
	NumberOfParameters(NumberOfParameterTypes - 1);
    }
  else
    {
      TotalNumberOfParameters = 0;
    }

  // Load the parameters
  
  Parameters.allocate(TotalNumberOfParameters);
  strm >> temp;
  
  int pi = Parameters.read(strm,0.0,"2 particle Jastrow");

  // Load the number of Constant types  
  strm >> temp >> temp;
  NumberOfConstantTypes = atoi(temp.c_str());

  // Load the number of Constants of each type
  
  NumberOfConstants.allocate(NumberOfConstantTypes);
  strm >> temp;
  
  for (int i = 0; i < NumberOfConstantTypes; i++)
    {
      strm >> temp;
      NumberOfConstants(i) = atoi(temp.c_str());
    }
  
  // Find the beginning index of each constant type
  
  BeginningIndexOfConstantType.allocate(NumberOfConstantTypes);
  
  for (int i = 0; i < NumberOfConstantTypes; i++)
    {
      if (i == 0)
 	{
 	  BeginningIndexOfConstantType(0) = 0;
 	}
      else
 	{
 	  BeginningIndexOfConstantType(i) =
 	    BeginningIndexOfConstantType(i - 1) + 
 	    NumberOfConstants(i-1);
 	}
    }
  
  // Set the total number of constants
  
  if( NumberOfConstantTypes > 0 )
    {
      TotalNumberOfConstants = BeginningIndexOfConstantType(
        NumberOfConstantTypes - 1) +
	NumberOfConstants(NumberOfConstantTypes - 1);
    }
  else
    {
      TotalNumberOfConstants = 0;
    }

  // Load the constants
  
  Constants.allocate(TotalNumberOfConstants);
  strm >> temp;
  
  for (int i = 0; i < TotalNumberOfConstants; i++)
    {
      strm >> temp;
      Constants(i) = atof(temp.c_str());
    }

  if( ParticleTypes(1) != "Electron_up" && ParticleTypes(1) != "Electron_down"
      && nucCuspReplacement && CorrelationFunctionType != "None")
    {
      /*
	If the electron nuclear cusp condition is already fulfilled
	(i.e. replace_electron_nucleus_cusps = 1), then the question
	is whether we want to automatically turn off any Electron-Nuclear
	Jastrow that might conflict with this.

	Previously, it would be automatically turned off. However, with the
	Cambridge Jastrow, it's probably advantageous to use an Electron Nucleus
	Jastrow with a cusp condition of 0.
      */
      //clog << "WARNING: switching CorrelationFunctionType to \"None\" for (" << pt1 << ", " << pt2 << ")" << 
      //" since we're replacing the cusps." << endl;
      //CorrelationFunctionType = "None";
    }

  // if the correlation function type is "None"
  // then remove any parameters or constants
  // that were listed
  
  if( CorrelationFunctionType == "None" )
    {
      NumberOfParameterTypes = 0;
      TotalNumberOfParameters = 0;
      NumberOfParameters.deallocate();
      Parameters.deallocate();
      BeginningIndexOfParameterType.deallocate();

      NumberOfConstantTypes = 0;
      TotalNumberOfConstants = 0;
      NumberOfConstants.deallocate();
      Constants.deallocate();
      BeginningIndexOfConstantType.deallocate();

    }

  // set the correlation function
  setCorrelationFunction();
  initializeCorrelationFunctionParameters();
  return ok;
}

QMCCorrelationFunctionParameters::QMCCorrelationFunctionParameters()
{
  CorrelationFunctionType = "None";
  CorrelationFunction = 0;
}

QMCCorrelationFunctionParameters::QMCCorrelationFunctionParameters(const 
			       QMCCorrelationFunctionParameters & rhs)
{
  CorrelationFunction = 0;
  * this = rhs;
}

string QMCCorrelationFunctionParameters::getParticle1Type()
{
  return ParticleTypes(0);
}

string QMCCorrelationFunctionParameters::getParticle2Type()
{
  return ParticleTypes(1);
}

int QMCCorrelationFunctionParameters::getTotalNumberOfParameters()
{
  return TotalNumberOfParameters;
}

void QMCCorrelationFunctionParameters::setParameters(Array1D<double> & params)
{
  if( params.dim1() != Parameters.dim1() )
    {
      cerr << "ERROR: Parameters of the incorrect size are trying to be set "
	   << "in QMCCorrelationFunctionParameters" << endl;
    }

  Parameters = params;
  initializeCorrelationFunctionParameters();
}

void QMCCorrelationFunctionParameters::setParticle1Type(string val)
{
  ParticleTypes(0) = val;
}

void QMCCorrelationFunctionParameters::setParticle2Type(string val)
{
  ParticleTypes(1) = val;
}

void QMCCorrelationFunctionParameters::setCorrelationFunction()
{
  if( CorrelationFunction != 0 )
    {
      delete CorrelationFunction;
      CorrelationFunction = 0;
    }
  CorrelationFunction = QMCCorrelationFunctionFactory::
    correlationFunctionFactory(CorrelationFunctionType);
}

QMCCorrelationFunction * QMCCorrelationFunctionParameters::
                                             getCorrelationFunction()
{
  return CorrelationFunction;
}

QMCCorrelationFunctionParameters::~QMCCorrelationFunctionParameters()
{
  NumberOfParameters.deallocate();
  Parameters.deallocate();
  ParticleTypes.deallocate();
  BeginningIndexOfParameterType.deallocate();
  
  if( CorrelationFunction != 0 )
    {
      delete CorrelationFunction;
      CorrelationFunction = 0;
    }
}

bool QMCCorrelationFunctionParameters::isSingular()
{
  return CorrelationFunction->isSingular();
}

void QMCCorrelationFunctionParameters::
                                     initializeCorrelationFunctionParameters()
{
  CorrelationFunction->initializeParameters(BeginningIndexOfParameterType,
					    Parameters,
					    BeginningIndexOfConstantType,
					    Constants);
}

ostream & operator << (ostream & strm, QMCCorrelationFunctionParameters & rhs)
{
  strm << "ParticleTypes:\t" << rhs.ParticleTypes(0) << "\t" 
       << rhs.ParticleTypes(1) << endl;

  strm << "CorrelationFunctionType:\t" << rhs.CorrelationFunctionType << endl;

  strm << "NumberOfParameterTypes:\t" << rhs.NumberOfParameterTypes << endl;

  strm << "NumberOfParametersOfEachType:";

  for (int i = 0; i < rhs.NumberOfParameterTypes; i++)
    {
      strm << "\t" << rhs.NumberOfParameters(i);
    }
  strm << endl;

  strm << "Parameters:";

  for (int i = 0; i < rhs.TotalNumberOfParameters; i++)
    {
      strm << "\t" << rhs.Parameters(i);
    }
  
  strm << endl;

  strm << "NumberOfConstantTypes:\t" << rhs.NumberOfConstantTypes << endl;

  strm << "NumberOfConstantsOfEachType:";

  for (int i = 0; i < rhs.NumberOfConstantTypes; i++)
    {
      strm << "\t" << rhs.NumberOfConstants(i);
    }
  strm << endl;

  strm << "Constants:";

  for (int i = 0; i < rhs.TotalNumberOfConstants; i++)
    {
      strm << "\t" << rhs.Constants(i);
    }

  strm << endl;
  
  return strm;
}
