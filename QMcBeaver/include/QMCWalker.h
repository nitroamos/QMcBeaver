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

#ifndef QMCWALKER_H
#define QMCWALKER_H

#include "IeeeMath.h"
#include "QMCFunctions.h"
#include "QMCInitializeWalkerFactory.h"
#include "QMCProperties.h"
#include "MathFunctions.h"
#include "QMCGreensRatioComponent.h"

using namespace std;

/**
  An instantaneous snapshot of all 3N electronic corrdinates for a system.
  This is the same as the "walker" or "psip" discussed in QMC literature.
*/

class QMCWalker
{
public:
  /**
    Creates a new uninitialized instance of this class.
  */
  QMCWalker();

  /**
    Creates a new instance of this class and makes it equivalent to another
    instance of this class.
    @param rhs object to set this equal to.
  */
  QMCWalker( const QMCWalker & rhs );

  /**
    Deallocates the memory allocated by this object.
  */
  ~QMCWalker();

  /**
    Initializes and allocates memory for the walker.  The electronic 
    configuration for the walker is not set.  To do this 
    <code>initializeWalkerPosition</code> must be used to generate a new
    walker, or <code>read</code> must be used to read this walkers state
    from a stream.
    @param input data input to control the calculation.
  */
  void initialize(QMCInput *input);

  /**
    Initializes the electronic configuration for this walker using an
    algorithm from QMCInitializeWalkerFactory.  If a singular walker is
    generated, upto 100 configurations are generated until one is not 
    singular.
  */
  void initializeWalkerPosition(); 

  /**
    Proposes a trial walker move and accepts or rejects it.
  */
  void propagateWalker();

  /**
    Calculates the observables for this walker and adds them to the input
    QMCProperties.
    @param props properties to which this walkers current observable values are
     added.
  */
  void calculateObservables( QMCProperties & props );

  /**
    Sets two QMCWalker objects equal.
    @param rhs object to set this object equal to.
  */
  void operator=(const QMCWalker & rhs );

  /**
    Gets the weight for this walker.
    @return weight for this walker.
  */
  double getWeight();

  /**
    Sets the weight for this walker.
    @param val value to set the weight equal to.
  */
  void setWeight(double val);

  /**
    Determines if the trial wavefunction is singular for this walker.
    @return <code>true</code> if the trial wavefunction is singular for 
    this walker, and <code>false</code> otherwise.
  */
  bool isSingular();

  /**
    Writes the state of this object to an XML stream.
    @param strm XML stream
  */
  void toXML(ostream& strm);

  /**
    Loads the state of this object from an XML stream.  The input stream must
    be formatted exactly like the output from <code>toXML</code> because
    it is not intelligent.
    @param strm XML stream
  */
  void readXML(istream& strm);

  /**
    Writes the state of this walker to a stream in a format that is suitable
    for correlated sampling.  This writes out more information than 
    <code>toXML</code> so that parts of the wavefunction do not have to
    be reevaluated every time properties are calculated using correlated
    sampling.
    @param strm stream to write correlated sampling information to.
  */ 
  void writeCorrelatedSamplingConfiguration(ostream& strm);  

  /**
    Gets the value of the local energy estimator for this walker.
  */
  double getLocalEnergyEstimator();

private:
  double weight;
  int age;

  double localEnergy;
  double kineticEnergy;
  double potentialEnergy;
  double distanceMovedAccepted;
  double AcceptanceProbability;

  QMCFunctions QMF;

  Array2D <double> R;

  // Magnitude of the proposed move
  double dR2;

  QMCWalker *TrialWalker;
  QMCWalker *OriginalWalker;

  bool move_accepted;

  QMCInput *Input;

  void setAcceptanceProbability(double p);

  void evaluate();
  void createChildWalkers();
  void calculateMoveAcceptanceProbability(double GreensRatio);
  void acceptOrRejectMove();

  /**
    Randomly moves the electrons to their new locations and returns the 
    Green's function for the forward move.
    @return Greens's function for the forward move of the electrons.
  */
  QMCGreensRatioComponent moveElectrons();

  /**
    Randomly moves the electrons to their new locations without using
    importance sampling.
    @return Greens's function for the forward move of the electrons.
  */
  QMCGreensRatioComponent moveElectronsNoImportanceSampling();

  /**
    Randomly moves the electrons to their new locations using
    importance sampling.
    @return Greens's function for the forward move of the electrons.
  */
  QMCGreensRatioComponent moveElectronsImportanceSampling();

  /**
    Randomly moves the electrons to their new locations using
    importance sampling and the accelerated Metropolis algorithm.
    @return Greens's function for the forward move of the electrons.
  */
  QMCGreensRatioComponent moveElectronsUmrigar93ImportanceSampling();

  /**
    Calculates the reverse Green's function for the proposed move of the
    electrons.
    @return Green's function for the reverse move of the electrons.
  */
  QMCGreensRatioComponent calculateReverseGreensFunction();

  /**
    Calculates the reverse Green's function for the proposed move of the
    electrons when no importance sampling is used.
    @return Green's function for the reverse move of the electrons.
  */
  QMCGreensRatioComponent calculateReverseGreensFunctionNoImportanceSampling();

  /**
    Calculates the reverse Green's function for the proposed move of the
    electrons when importance sampling is used.
    @return Green's function for the reverse move of the electrons.
  */
  QMCGreensRatioComponent calculateReverseGreensFunctionImportanceSampling();

  /**
    Calculates the reverse Green's function for the proposed move of the
    electrons when umrigar93 importance sampling is used.
    @return Green's function for the reverse move of the electrons.
  */
  QMCGreensRatioComponent \
                   calculateReverseGreensFunctionUmrigar93ImportanceSampling();

  // reweight the walker after a move
  void reweight_walker();

  Array2D<double> * getR();
  int getAge();
  double getAcceptanceProbability();

  /**
    Calculates the observables for this walker.
  */
  void calculateObservables();
};

#endif






