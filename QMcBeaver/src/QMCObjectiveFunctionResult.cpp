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

#include "QMCObjectiveFunctionResult.h"

QMCObjectiveFunctionResult::QMCObjectiveFunctionResult()
{
}

QMCObjectiveFunctionResult::QMCObjectiveFunctionResult(QMCInput *input, 
						       double energyAve, 
						       double energyVar, 
						       double logWeightAve,
						       double logWeightVar,
						      Array1D<Complex> & poles)
{
  Input = input;

  set_energy_ave( energyAve );
  set_energy_var( energyVar );
  set_log_weights_ave( logWeightAve );
  set_log_weights_var( logWeightVar );

  set_score(poles);
  set_score_for_derivative(poles);
}

QMCObjectiveFunctionResult::QMCObjectiveFunctionResult(
					 QMCObjectiveFunctionResult & rhs)
{
  *this = rhs;
}

double QMCObjectiveFunctionResult::getLogWeightsAve()
{
  return log_weights_ave;
}

double QMCObjectiveFunctionResult::getLogWeightsVar()
{
  return log_weights_var;
}

double QMCObjectiveFunctionResult::getEnergyAve()
{
  return energy_ave;
}

double QMCObjectiveFunctionResult::getEnergyVar()
{
  return energy_var;
}

double QMCObjectiveFunctionResult::getScore()
{
  return score;
}

double QMCObjectiveFunctionResult::getDerivativeScore()
{
  return score_for_derivative;
}

void QMCObjectiveFunctionResult::set_log_weights_ave(double wa)
{
  if(IeeeMath::isnan(wa) != 0) 
    {
      log_weights_ave = MAX_RESULT_VALUE;
    }
  else if(fabs(wa) > MAX_RESULT_VALUE)
    {
      log_weights_ave = MAX_RESULT_VALUE;
    }
  else
    {
      log_weights_ave = wa;
    }
}

void QMCObjectiveFunctionResult::set_log_weights_var(double wv)
{
  if(IeeeMath::isnan(wv) != 0) 
    {
      log_weights_var = MAX_RESULT_VALUE;
    }
  else if(fabs(wv) > MAX_RESULT_VALUE)
    {
      log_weights_var = MAX_RESULT_VALUE;
    }
  else if(wv < 0.0)
    {
      log_weights_var=MAX_RESULT_VALUE;
    }
  else
    {
      log_weights_var = wv;
    }
}

void QMCObjectiveFunctionResult::set_energy_ave(double ea)
{
  if(IeeeMath::isnan(ea) != 0) 
    {
      energy_ave = MAX_RESULT_VALUE;
    }
  else if(fabs(ea) > MAX_RESULT_VALUE)
    {
      energy_ave = MAX_RESULT_VALUE;
    }
  else
    {
      energy_ave=ea;
    }
}

void QMCObjectiveFunctionResult::set_energy_var(double ev)
{
  if(IeeeMath::isnan(ev) != 0) 
    {
      energy_var=MAX_RESULT_VALUE;
    }
  else if(fabs(ev) > MAX_RESULT_VALUE)
    {
      energy_var=MAX_RESULT_VALUE;
    }
  else if(ev < 0.0)
    {
      energy_var=MAX_RESULT_VALUE;
    }
  else
    {
      energy_var=ev;
    }
}

void QMCObjectiveFunctionResult::operator=(QMCObjectiveFunctionResult &rhs)
{
  this->score                = rhs.score;
  this->score_for_derivative = rhs.score_for_derivative;
  this->log_weights_ave      = rhs.log_weights_ave;
  this->log_weights_var      = rhs.log_weights_var;
  this->energy_ave           = rhs.energy_ave;
  this->energy_var           = rhs.energy_var;
  this->Input                = rhs.Input;
}

void QMCObjectiveFunctionResult::set_score(Array1D<Complex> & poles)
{
  //take these properties and make a score out of them
  //this is where we can add our penalty function stuff if needed
 
  if(Input->flags.optimize_Psi_criteria=="energy_variance")
    {
      score = getEnergyVar();
    }
  else if(Input->flags.optimize_Psi_criteria=="energy_average")
    {
      score = getEnergyAve();
    }
  else if(Input->flags.optimize_Psi_criteria=="umrigar88")
    {
      score = calculate_umrigar88();
    }
  else if(Input->flags.optimize_Psi_criteria=="monkey_spank")
    {
      score = calculate_monkey_spank();
    }
  else
    {
      cerr<<"ERROR: in Input.flags.optimize_Psi_criteria"<<endl;
      cerr<<Input->flags.optimize_Psi_criteria<<endl;
      exit(1);
    }

  // add the penalty function
  score += calculate_penalty_function(poles);
}

void QMCObjectiveFunctionResult::set_score_for_derivative(
						     Array1D<Complex> & poles)
{
  //take these properties and make a score out of them
  //this is where we can add our penalty function stuff if needed
 
  if(Input->flags.numerical_derivative_surface=="energy_variance")
    {
      score_for_derivative = getEnergyVar();
    }
  else if(Input->flags.numerical_derivative_surface=="energy_average")
    {
      score_for_derivative = getEnergyAve();
    }
  else if(Input->flags.numerical_derivative_surface=="umrigar88")
    {
      score_for_derivative = calculate_umrigar88();
    }
  else if(Input->flags.numerical_derivative_surface=="monkey_spank")
    {
      score_for_derivative = calculate_monkey_spank();
    }
  else
    {
      cerr<<"ERROR: in Input.flags.numerical_derivative_surface"<<endl;
      cerr<<Input->flags.numerical_derivative_surface<<endl;
      exit(1);
    }

  // add the penalty function
  score_for_derivative += calculate_penalty_function(poles);
}

//This will calculate the penalty function and return a score for a 
//parameter set
double QMCObjectiveFunctionResult::mikes_score_function()
{
  //It is assumed that: 
  //energy_ave has no inherent boundaries but we would like it to get smaller
  //energy_var>=0.0 and we would like it to go to zero
  //log_weights_ave has no inherent boundaries but we would like to see it 
  //go to 0.0 log_weights_var>=0.0 and we would like to keep it near 0.0

  //log_weights_ave=dev_log_weights_ave is where the penalty function will 
  //start to dominate
  //  double dev_log_weights_ave 
  //      = 1.0*Input->flags.optimize_Psi_barrier_parameter;

  //log_weights_var=dev_log_weights_var is where the penalty function will 
  //start to dominate
  double dev_log_weights_var = 0.1*Input->flags.optimize_Psi_barrier_parameter;

  double LV;       //a measure of the lack of validity for a particular input

  /////////////////////////////////////////////////
  //ACTUALLY START CALCULATING PENALTIES AND SCORES
  /////////////////////////////////////////////////

  //goal in constructing f1,f2,f3,f4 is to keep them all fi(x)>=1.0
  double f1,f2,f3,f4;
  //calculate the contribution of energy_ave
  f1=1.0; //no contribution to the score

  //calculate the contribution of energy_var
  f2=1.0+energy_var;  //f2>=1.0 everywhere 

  //calculate the contribution of log_weights_ave
  //LV = mikes_penalty_scaler(log_weights_ave/dev_log_weights_ave);
  //f3 = mikes_penalty(LV);
  f3 = 1.0;

  //calculate the contribution of log_weights_var
  LV = mikes_penalty_scaler(log_weights_var/dev_log_weights_var);
  f4 = mikes_penalty(LV);

  //since all the fi>=1.0 the scores will be above 1.0
#ifdef DEBUG_MIKES_FUNCTION
  cerr << "mike's function: " << f1 << "\t" << f2 << "\t" << f3
       << "\t" << f4 << endl;
#endif

  return f1*f2*f3*f4;
}

//this will take a number and scale it in the following way
//0->1
//0.5->very near 1.0
//1->2
//higher->blows up
double QMCObjectiveFunctionResult::mikes_penalty(double x)
{
  double result;
  int steepness=4;
  result=exp(log(2.0)*pow(x,steepness*2.0));
  return result;
}

//this will take a number and scale it downward in the following way:
//near linear under one
//1->1
//bounded by MAX for x->infinity
double QMCObjectiveFunctionResult::mikes_penalty_scaler(double x)
{
  double MAX=1.5;
  //this is so low because it will work with penalty()
  //to give us a bound on how large the penalties will get
  //once the LV goes beyond 1.3 with a steepness of 4 is rather
  //nutty and pushes machine precision much too far

  if( IeeeMath::isnan(x) != 0 ) return MAX;
  if( fabs(x) >= MAX) return MAX;

  if( IeeeMath::isnan(x) != 0 )
    {
      cerr << x << endl;
      cerr <<"ERROR QMCObjectiveFunctionResult::mikes_penalty_scaler(double x)"
	   << endl;
      exit(1);
    }
  double result = MAX*(1.0-exp(log(1.0-1.0/MAX)*fabs(x)));
  return result;
}

// From our ass...
double QMCObjectiveFunctionResult::calculate_monkey_spank()
{
  return mikes_score_function();
}

// From umrigar's 1998 phys rev let paper on QMC optimization
double QMCObjectiveFunctionResult::calculate_umrigar88()
{
  double temp = Input->flags.energy_estimated_original - getEnergyAve();
  double val = getEnergyVar() + temp*temp;
  return val;
}

// calculates a penalty function for getting singular parameters
double QMCObjectiveFunctionResult::calculate_penalty_function(
						     Array1D<Complex> & poles)
{
  double penalty = 0.0;

  for(int i=0; i<poles.dim1(); i++)
    {
      // calculate the distance of the pole from the positive real axis
      double distance = 0.0;

      if( poles(i).real() > 0 )
	{
	  distance = fabs(poles(i).imaginary());
	}
      else
	{
	  distance = poles(i).abs();
	}

      penalty -= log( distance );
    }

  return penalty * Input->flags.singularity_penalty_function_parameter;
}

ostream& operator<<(ostream & strm, QMCObjectiveFunctionResult & rhs)
{
  strm << "log(weights):     " << rhs.getLogWeightsAve() << "+/-" 
       << sqrt(rhs.getLogWeightsVar()) << endl;

  strm << "energy:           " << rhs.getEnergyAve() << "+/-"
       << sqrt(rhs.getEnergyVar()) << endl;

  strm << "score:            " << rhs.getScore() << endl;
  strm << "derivative score: " << rhs.getDerivativeScore() << endl;

  return strm;
}
