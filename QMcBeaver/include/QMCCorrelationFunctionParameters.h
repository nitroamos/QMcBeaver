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

#ifndef QMCCorrelationFunctionParameters_H
#define QMCCorrelationFunctionParameters_H

#include <iostream>
#include <string>

#include "Array1D.h"
#include "QMCCorrelationFunctionFactory.h"
#include "StringManipulation.h"

using namespace std;


/**
  This is a collection of parameters and related functions which describe
  the interaction of two particles of specific types.  For example, an instance
  of this class could hold the information describing the interaction of
  an up spin electron and a hydrogen nucleus or two down spin electrons.

  The interactions are parameterized in terms of "parameters" and "constants."
  "parameters" are modified during optimizations, and "constants" are not.
  */

class QMCCorrelationFunctionParameters
{
public:
  /**
    Creates an instance of the class.
    */

  QMCCorrelationFunctionParameters();


  /** 
    Creates an instance of the class that is identical to another instance
    of the class.

    @param rhs object to copy
    */

  QMCCorrelationFunctionParameters(const QMCCorrelationFunctionParameters & 
				   rhs);


  /**
    Deallocates all of the memory used by the object and prepares it to be 
    destroyed.
    */

  ~QMCCorrelationFunctionParameters();


  /**
    Gets the parameters describing the particle-particle interactions.

    @return parameters describing particle-particle interactions.
    */

  Array1D < double > getParameters();
  
  /**
     Gets the poles of the correlation function.

     @return poles of the correlation function.
  */
  Array1D < Complex > getPoles();


  /**
    Gets the first particle type in a particle1-particle2 interaction described
    by this object.

    @return particle type
    */

  string getParticle1Type();


  /** 
    Gets the second particle type in a particle1-particle2 interaction 
    described  by this object.

    @return particle type
    */

  string getParticle2Type();


  /**
    Gets the total number of parameters used to describe the 
    particle-particle  interaction.

    @return total number of parameters
    */

  int getTotalNumberOfParameters();


  /**
    Gets the parameterized QMCCorrelationFunction used in QMCJastrow to
    describe the particular particle-particle interaction when
    calculating the Jastrow function.

    @return function describing getParticle1Type()-getParticle2Type() 
    interactions
    */

  QMCCorrelationFunction * getCorrelationFunction();    
  
  /**
     Is this correlation function "None"?
  */
  bool isUsed();

  /**
    Sets the parameters describing the particle-particle interaction.

    @param params new set of parameters
    */

  void setParameters(Array1D<double> & params);


  /**
    Sets the type of particle1 for the particular particle-particle
    interaction described by this object.
    */
  
  void setParticle1Type(string val);
    
  
  /**
    Sets the type of particle2 for the particular particle-particle
    interaction described by this object.
    */

  void setParticle2Type(string val);    
  

  /**
    Returns true if the parameterized correlation function described
    by this object is singular on the positive real axis and false
    otherwise.

    @return true if the current parameterization of the correlation 
    function is singular on the positive real axis and false otherwise
    */

  bool isSingular();


  /**
    Sets two QMCCorrelationFunctionParameters objects equal.

    @param rhs object to set this object eqal to
    */

  void operator = (const QMCCorrelationFunctionParameters & rhs);


  /**
    Loads the state of the object from an input stream.

    @param strm input stream
    @param nucCuspReplacement indicates whether we need to switch
    off any electron-nucleus jastrows from the input file
    */
  
  bool read(istream & strm, bool nucCuspReplacement);


  /**
    Writes the state of the object to an output stream.
    */

  friend ostream & operator << (ostream & strm, 
				QMCCorrelationFunctionParameters & rhs);

private:
  void initializeCorrelationFunctionParameters();    
  void setCorrelationFunction();
  
  Array1D < string > ParticleTypes;

  int NumberOfParameterTypes;
  Array1D < int > NumberOfParameters;
  Array1D < double > Parameters;
  Array1D < int > BeginningIndexOfParameterType;
  int TotalNumberOfParameters;

  int NumberOfConstantTypes;
  Array1D < int > NumberOfConstants;
  Array1D < double > Constants;
  Array1D < int > BeginningIndexOfConstantType;
  int TotalNumberOfConstants;

  string CorrelationFunctionType;
  QMCCorrelationFunction* CorrelationFunction;
};

#endif
