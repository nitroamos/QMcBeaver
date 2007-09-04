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

#ifndef QMCThreeBodyCorrelationFunctionParameters_H
#define QMCThreeBodyCorrelationFunctionParameters_H

#include <iostream>
#include <string>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCThreeBodyCorrelationFunctionFactory.h"
#include "StringManipulation.h"

using namespace std;

/**
  This is a collection of parameters and related functions which describe
  the interaction of a nucleus and two electrons.  For example, an instance
  of this class could hold the information describing the interaction of
  a hydrogen nucleus, an up spin electron, and a down spin electron.

  The interactions are parameterized in terms of "parameters" and "constants."
  "parameters" are modified during optimizations, and "constants" are not.
*/

class QMCThreeBodyCorrelationFunctionParameters
{
 public:
  /**
    Creates an instance of the class.
  */
  QMCThreeBodyCorrelationFunctionParameters();

  /** 
    Creates an instance of the class that is identical to another instance
    of the class.

    @param rhs object to copy
  */
  QMCThreeBodyCorrelationFunctionParameters(
		        const QMCThreeBodyCorrelationFunctionParameters & rhs);

  /**
    Deallocates all of the memory used by the object and prepares it to be 
    destroyed.
  */
  ~QMCThreeBodyCorrelationFunctionParameters();

  /**
    Gets the free parameters describing the three body interactions.
    Not all of the parameters are free due to cusp and symmetry constraints.

    @return free parameters describing the three body interactions
  */
  Array1D<double> getFreeParameters();

  /**
    Gets the first particle type in the three body interaction described
    by this object.  The first particle should be a nucleus.

    @return particle1 type
  */
  string getParticle1Type();

  /** 
    Gets the second particle type in the three body interaction 
    described  by this object.  This particle should be an up or down electron.

    @return particle2 type
  */
  string getParticle2Type();

  /** 
    Gets the third particle type in the three body interaction 
    described  by this object.  This particle should be an up or down electron.

    @return particle3 type
  */
  string getParticle3Type();

  /**
    Gets the total number of free parameters used to describe the three body
    interaction.  Not all of the parameters are free because of cusp and
    symmetry constraints

    @return total number of free parameters
  */
  int getNumberOfFreeParameters();

  /**
     Gets the total number of parameters.
     @return total number of parameters
  */
  int getNumberOfTotalParameters();

  /**
    Gets the parameterized QMCThreeBodyCorrelationFunction used in QMCJastrow 
    to describe the particular three body interaction when calculating the 
    Jastrow function.

    @return function describing the three body interactions
  */

  QMCThreeBodyCorrelationFunction * getThreeBodyCorrelationFunction();    

  /**
    Sets the type of particle1 for the particular three body interaction 
    described by this object.
  */  
  void setParticle1Type(string val);
    
  /**
    Sets the type of particle2 for the particular particle-particle
    interaction described by this object.
  */
  void setParticle2Type(string val);    

  /**
    Sets the type of particle3 for the particular particle-particle
    interaction described by this object.
  */
  void setParticle3Type(string val);    

  /**
    Sets two QMCThreeBodyCorrelationFunctionParameters objects equal.

    @param rhs object to set this object equal to
  */
  void operator = (const QMCThreeBodyCorrelationFunctionParameters & rhs);

  /**
    Loads the state of the object from an input stream.

    @param strm input stream
  */
  bool read(istream & strm);

  /**
    Writes the state of the object to an output stream.
  */
  friend ostream & operator << (ostream & strm, 
			      QMCThreeBodyCorrelationFunctionParameters & rhs);

  /**
    Gets the cutoff distance for this three body correlation function.
  */
  double getCutoffDist();

  void zeroOutDerivatives();
  Array1D<double> pt_a;
  Array1D< Array2D<double> > pt2_xa;
  Array1D<double> pt3_xxa;

  //private:
  void initializeThreeBodyCorrelationFunctionParameters();    
  void setThreeBodyCorrelationFunction();
  
  /**
     convert an array in the total basis into an array in the free basis
     free -> total (transformed)
  */
  void getFree(const Array1D<double> & total,
	       Array1D<double> & free);

  /**
     takes a new set of free parameters and recreates the total parameter array
     free -> total (transformed)
  */
  void setFreeParameters(const Array1D<double> & free);

  /**
     takes a new set of free parameters and recreates the total parameter array
     total -> total (transformed) -> free
  */
  void setTotalParameters(const Array1D<double> & total);
  
  /**
     maps the (transformed) total parameter array to the free array
     (transformed) total -> free
  */
  void totalToFree(const Array1D<double> & total,
		   Array1D<double> & free) const;

  void totalDerivativesToFree(int shift,
			      Array1D<double> & p_a,
			      Array1D< Array2D<double> > & p2_xa,
			      Array1D<double> & p3_xxa) const;

  /**
     returns true if the coefficient associated with
     this particular l, m, and n is free.
  */
  bool isFree(int l, int m, int n) const;

  /**
     maps the free parameter array to the total parameter array
     free -> total (transformed)
  */
  void freeToTotal(const Array1D<double> & free,
		   Array1D<double> & total);

  /**
     creates the transformation matrix which can apply the cusp and symmetry conditions
     to the total parameter array. this creates certain dependencies between the parameters.
     total -> (transformed) total
  */
  void makeParamDepMatrix();

  /**
     verifies that total parameter array statisfies the cusp and symmetry conditions
     (transformed) total -> true
     total -> false
  */
  bool checkCuspAndSymmetry();

  Array1D<string> ParticleTypes;

  int NumberOfParameterTypes;
  Array1D<int> NumberOfParameters;
  Array1D<double> Parameters;
  Array1D<double> tempArray;
  int TotalNumberOfParameters;

  Array1D<double> freeParameters;
  int NumberOfFreeParameters;

  int NeN;
  int Nee;

  int C;
  double cutoff;

  Array2D<double> paramDepMatrix;

  string threeBodyCorrelationFunctionType;
  QMCThreeBodyCorrelationFunction* ThreeBodyCorrelationFunction;
};

#endif