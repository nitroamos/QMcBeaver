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

#ifndef QMCJastrowParameters_H
#define QMCJastrowParameters_H

#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "Array1D.h"
#include "QMCCorrelationFunctionParameters.h"

using namespace std;


/**
  This class contains all of the parameters and correlation functons from 
  which the Jastrow function is composed.
 
  The wavefunction is assumed to be of the form
  \f[
  \Psi_{QMC} = \Psi_{Trial}J
  \f]
  where \f$\Psi_{Trial}\f$ is a wavefunction calculated using a standard QM
  method and
  \f[
  J=exp(\sum{u_{i,j}(r_{i,j})})
  \f]
  is a Jastrow type correlation function.  \f$u_{ij}(r_{ij})\f$ are 
  QMCCorrelationFunction describing the interactions of particles 
  \f$i\f$ and \f$j\f$.  The correlation functions are parameterized to 
  allow optimization.  This class contains the functions and their specific
  parameterizations.   The interactions are parameterized in terms of 
  "parameters" and "constants."  "parameters" are modified during 
  optimizations, and "constants" are not.
 */

class QMCJastrowParameters
{
public:
  /**
    Creates an instance of the class.
  */
  QMCJastrowParameters();

  /** 
    Creates an instance of the class that is identical to another instance
    of the class.

    @param rhs object to copy
  */
  QMCJastrowParameters(const QMCJastrowParameters & rhs);

  /**
    Sets the parameters describing the particle-particle interactions.

    @param params new set of parameters
  */
  void setParameterVector(Array1D<double> & params);

  /**
    Gets the parameters describing the particle-particle interactions.

    @return parameters describing particle-particle interactions.
  */
  Array1D<double> getParameters();

  int getNumberEEParameters();
  int getNumberNEParameters();
  int getNumberNEupParameters();

  /**
     Gets the poles of the correlation functions.

     @return poles of the correlation functions.
  */
  Array1D<Complex> getPoles();

  // calculates a penalty function for getting singular parameters
  double calculate_penalty_function();

  // calculates a penalty function for getting singular parameters
  static double calculate_penalty_function(Array1D<Complex> & poles);

  /**
    Gets the QMCCorrelationFunctionParameters describing up-down electron
    interactions.

    @return up-down electron interaction parameters
  */
  QMCCorrelationFunctionParameters * getElectronUpElectronDownParameters();

  /**
    Gets the QMCCorrelationFunctionParameters describing up-up electron
    interactions.

    @return up-up electron interaction parameters
  */
  QMCCorrelationFunctionParameters * getElectronUpElectronUpParameters();


  /**
    Gets the QMCCorrelationFunctionParameters describing down-down electron
    interactions.

    @return down-down electron interaction parameters
  */
  QMCCorrelationFunctionParameters * getElectronDownElectronDownParameters();


  /**
    Gets an array of QMCCorrelationFunctionParameters describing 
    up electron-nuclear interactions.

    @return up electron-nuclear interaction parameters
  */
  Array1D<QMCCorrelationFunctionParameters> * getElectronUpNuclearParameters();


  /**
    Gets an array of QMCCorrelationFunctionParameters describing 
    down electron-nuclear interactions.

    @return down electron-nuclear interaction parameters
  */
  Array1D<QMCCorrelationFunctionParameters> * 
                                            getElectronDownNuclearParameters();
 
  /**
    Gets an array which is a list of all the different types of nuclei
    in the molecule being calculated.
  */
  Array1D<string> * getNucleiTypes();

  /**
    Sets two QMCJastrowParameters objects equal.

    @param rhs object to set this object eqal to
  */
  void operator=( const QMCJastrowParameters & rhs );

  /**
    Loads the state of the object from a file.

    @param nucleitypes list of the different kinds of nuclei
    @param linkparams true if nuclear-electron interactions are strictly the 
    same and false otherwise
    @param nucCuspReplacement indicates whether we need to shut off any
    electron-nucleus Jastrow functions in the input file
    @param nelup number of up spin electrons
    @param neldn numer of down spin electrons
    @param runfile name of the file to be loaded
  */
  void read(Array1D<string> & nucleitypes, bool linkparams, bool nucCuspReplacement,
	    int nelup, int neldn, string runfile);

  /**
    Writes the state of the object to an output stream.
  */
  friend ostream & operator<<(ostream &strm, QMCJastrowParameters & rhs);
    
private:
  int NumberOfEEParameters;
  int NumberOfNEParameters;
  int NumberOfNEupParameters;
  Array1D<QMCCorrelationFunctionParameters> EupNuclear;
  Array1D<QMCCorrelationFunctionParameters> EdnNuclear;
  QMCCorrelationFunctionParameters EupEdn;
  QMCCorrelationFunctionParameters EupEup;
  QMCCorrelationFunctionParameters EdnEdn;
  bool EquivalentElectronUpDownParams;
  Array1D<string> NucleiTypes;
  int NumberOfElectronsUp;
  int NumberOfElectronsDown;
};

#endif
