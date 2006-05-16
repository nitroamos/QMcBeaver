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
#include "QMCFutureWalkingProperties.h"

using namespace std;

/**
  An instantaneous snapshot of all 3N electronic coordinates for a system.
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
    generated, up to 100 configurations are generated until one is not 
    singular.
  */
  void initializeWalkerPosition(QMCFunctions & QMF); 

  /**
    Proposes a trial walker move and accepts or rejects it. This method has
    been broken into 2 parts. The first part (this function) moves the
    electrons and then returns pointers to the new positions along
    with pointers to the QMCWalkerData structs so they can be filled.

    The second part is the processPropagation function in this class.
    See QMCRun::propagateWalkers for how the two functions work together.

    The whole point of the * & (reference to pointer) is that the two
    parameters are pointers, and are both outputs of this function.

    Another way to do this might be to pass in the actual array of
    pointers, and tell initializePropagation which index to assign,
    but this was a good exercise in understanding pointers...

    @param (output) data to put a pointer to this walkerData
    @param (output) R to put a pointer to this' new configuration
  */
  void initializePropagation(QMCWalkerData * &data, Array2D<double> * &R);

  /**
     This function completes the processing. The forwardGreensFunction was
     stored and because pointers were given with a call to 
     initializePropagation, no parameters need to be passed. This function 
     should not be called without first calling initializePropagation.
  */
  void processPropagation(QMCFunctions & QMF);

  /**
    Calculates the observables for this walker and adds them to the input
    QMCProperties.
    @param props properties to which this walkers current observable values are
     added.
  */
  void calculateObservables( QMCProperties & props );

  /**
    Calculates the observables for this walker and adds them to the input
    QMCFutureWalkingProperties.
    @param props properties to which this walkers current observable values are
     added.
  */
  void calculateObservables( QMCFutureWalkingProperties & props );

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
    Calculates the distance between each pair of electrons and records it in 
    the appropriate histogram.  This will be used to evaluate pair density 
    functions for DFT development.
  */
  void calculateElectronDensities(double max_pair_distance, double dr, 
			  Array1D<double> &pll_spin, Array1D<double> &opp_spin,
                                     Array1D< Array1D<double> > &alpha_density,
                                     Array1D< Array1D<double> > &beta_density);

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
    Gets the value of the local energy estimator for this walker.
  */
  double getLocalEnergyEstimator();

  /**
    Gets the positions of the electrons.
  */
  Array2D<double> * getR();
  
  /**
    Sets the positions of the electrons.
    @param temp_R the positions of the electrons
  */
  void setR(Array2D<double>& temp_R);

  /**
    Gets the walkerData for this walker
    @return walkerData for this walker
  */
  QMCWalkerData* getWalkerData();
  
  void resetFutureWalking(int whichStage, int whichBlock);
  void resetFutureWalking();
    
private:
  double weight;
  int age;
  
  /**
    These energies are calculated by QMCWalker using the acceptance probability
    as opposed to those calculated by QMCFunction
  */
  double localEnergy;
  double kineticEnergy;
  double potentialEnergy;
  double r12;
  double r2;

  double neEnergy;
  double eeEnergy;

  static const double maxFWAsymp;
  
  Array1D<int> numFWSteps;
  
  enum fwStage { ACCUM, ASYMP, DONE };
  
  Array2D<fwStage> isCollectingFWResults;  

  /*
    These have dimensions of (numFW, 2)
  */
#ifdef USE_QMCPROPERTY
  Array2D<QMCProperty> fwNormalization;
  Array2D<QMCProperty> fwR12;
  Array2D<QMCProperty> fwR2;
  Array2D<QMCProperty> fwKineticEnergy;
  Array2D<QMCProperty> fwPotentialEnergy;
  
  // (nuc dim1, nuc dim2) (numFW, 2)
  Array2D< Array2D<QMCProperty> > fwNuclearForces;
#else
  Array2D<double> fwNormalization;
  Array2D<double> fwEnergy;
  Array2D<double> fwKineticEnergy;
  Array2D<double> fwPotentialEnergy;
  Array2D<double> fwR12;
  Array2D<double> fwR2;
  
  // (nuc dim1, nuc dim2) (numFW, 2)
  Array2D< Array2D<double> > fwNuclearForces;
#endif
  
  double distanceMovedAccepted;
  double AcceptanceProbability;

  /**
     walkerData is meant to hold all the necessary data given by a QMCFunction
     to complete the iteration.
  */
  QMCWalkerData walkerData;
  
  /**
     This data is simply meant to carry the forwardGreensFunction between the 2
     propagate walker functions.
  */
  QMCGreensRatioComponent forwardGreensFunction;

  Array2D<double> R;

  // Magnitude of the proposed move
  double dR2;

  QMCWalker *TrialWalker;
  QMCWalker *OriginalWalker;

  bool move_accepted;

  QMCInput* Input;

  void setAcceptanceProbability(double p);
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

  int getAge();
  double getAcceptanceProbability();
  
  /**
    Reweight the walker after a move
  */
  void reweight_walker();

  /**
    Calculates the observables for this walker.
  */
  void calculateObservables();
};

#endif
