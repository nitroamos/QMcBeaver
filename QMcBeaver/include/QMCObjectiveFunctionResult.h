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

#ifndef QMCOBJECTIVEFUNCTIONRESULTS_H
#define QMCOBJECTIVEFUNCTIONRESULTS_H

#include "IeeeMath.h"
#include "QMCInput.h"

//This number will be returned as a result if any number gets
//too big or small or nan.
#define MAX_RESULT_VALUE 1e8

/**
   Results from the evaluation of an objective function during a QMC 
   calculation.  These results can then be used for numerical optimization
   or other functions.
*/

class QMCObjectiveFunctionResult 
{
 public:
  /**
     Creates a new uninitialized instance of this class.
  */
  QMCObjectiveFunctionResult();
  ~QMCObjectiveFunctionResult();

  /**
     Creates and initializes a new instance of this class.

     @param input data input to control the calculation.
     @param energyAve calculated energy value
     @param energyVar calculated energy variance
     @param logWeightAve average value of the natural log of the 
     statistical weights of the configurations.
     @param logWeightVar variance in the above quantity.
     @param poles of the correlation functions
  */ 
  QMCObjectiveFunctionResult(QMCInput *input,
			     double energyAve, double energyVar,
			     double logWeightAve, double logWeightVar, 
			     int numSamples,
			     Array1D<Complex> & poles);

  /**
     Creates a new instance of this class and makes it equivalent to another
     instance of this class.

     @param rhs object to set this equal to.
  */
  QMCObjectiveFunctionResult(QMCObjectiveFunctionResult & rhs);

  /**
     Gets the average value of the natural log of the statistical weights for 
     the configurations used in this function evaluation.

     @return average value of the natural log of the statistical weights.
  */
  double getLogWeightsAve();

  /**
     Gets the variance of the natural log of the statistical weights for 
     the configurations used in this function evaluation.

     @return variance of the natural log of the statistical weights.
  */
  double getLogWeightsVar();

  /**
     Gets the calculated average energy value.

     @return calculated average energy value.
  */
  double getEnergyAve() const;

  /**
     Gets the calculated energy variance.

     @return calculated energy variance.
  */
  double getEnergyVar() const;

  int getNumberSamples() const;

  /**
     Gets a score for this function evaluation.  Better scores have lower
     values.  The algorithm used for arriving at the scoris is determined
     by the input data.  The convergence of a numerical optimization can
     be modified by changing the score functions.

     @return score for the function evaluation.
  */
  double getScore() const;
  
  /**
     Gets a score for this function evaluation that is to be used in
     calculating the derivative in a numerical optimization.  
     The algorithm used for arriving at this score is determined
     by the input data.  The convergence of a numerical optimization can
     be modified by changing the score functions.

     @return score for the derivative evaluation.
  */
  double getDerivativeScore();
  
  /**
     Sets two QMCObjectiveFunctionResult objects equal.

     @param rhs object to set this object equal to.
  */
  void operator=(const QMCObjectiveFunctionResult &rhs);

  /**
     Prints the contents of this object in a human readable format.
  */
  friend ostream& operator<<(ostream & strm, 
			     const QMCObjectiveFunctionResult & rhs);

 private:
  int numSamples;
  double log_weights_ave;
  double log_weights_var;
  double energy_ave;
  double energy_var;

  double score;
  double score_for_derivative;

  QMCInput *Input;

  void set_log_weights_ave(double wa);
  void set_log_weights_var(double wv);
  void set_energy_ave(double ea);
  void set_energy_var(double ev);

  void set_score(Array1D<Complex> & poles);
  void set_score_for_derivative(Array1D<Complex> & poles);

  //This will calculate the penalty function and return a score for a 
  //parameter set
  double mikes_score_function();
  double mikes_penalty(double x);
  double mikes_penalty_scaler(double x);

  double calculate_monkey_spank();
  double calculate_umrigar88();
};

#endif
