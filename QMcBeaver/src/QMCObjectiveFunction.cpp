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

#include "QMCObjectiveFunction.h"

void QMCObjectiveFunction::initialize(QMCInput *Ip, int configsToSkip)
{
  Input   = Ip;
  RAEC.initialize(Input, configsToSkip);
}

Array1D<QMCObjectiveFunctionResult> QMCObjectiveFunction::evaluate(
		                          Array1D< Array1D<double> > & Params) 
{ 
  // calculate properties from the configuration file
  Array1D<QMCProperties> Properties;
  RAEC.rootCalculateProperties(Params, Properties);

  // process the properties and put the results into RESULTS
  Array1D<QMCObjectiveFunctionResult> RESULTS(Params.dim1());
  
  for(int i=0;i<Params.dim1();i++)
    {
      Input->JP.setParameterVector(Params(i));
      Array1D<Complex> poles = Input->JP.getPoles();

      QMCObjectiveFunctionResult newResult(Input,
		   Properties(i).energy.getAverage(),
		   Properties(i).energy.getSeriallyCorrelatedVariance(),
		   Properties(i).logWeights.getAverage(),
		   Properties(i).logWeights.getSeriallyCorrelatedVariance(),
		   poles);

      RESULTS(i) = newResult;
    }

  return RESULTS;
}

// Return the result of the objective function evaluation
QMCObjectiveFunctionResult QMCObjectiveFunction::evaluate(Array1D<double> &
						                        Params)
{
  Array1D< Array1D<double> > P(1);
  P(0) = Params;

  Array1D< QMCObjectiveFunctionResult > Result; 
  Result = evaluate(P);
  QMCObjectiveFunctionResult Final_Result = Result(0);

  return Final_Result;
}

// The gradient of the objective function
Array1D< Array1D<double> > QMCObjectiveFunction::grad(Array1D< 
						   Array1D<double> > &Params)
{
  cerr << "ERROR: Call to QMCObjectiveFunction::grad!" << endl;
  cerr << "This function is not currently implemented." << endl;
  exit(0);

  return 0;
}
  
// The gradient of the objective function
Array1D<double> QMCObjectiveFunction::grad(Array1D<double> &Params)
{
  // By default use a numerical gradient
  Array1D<double> GRAD;
  numerical_grad(Params, GRAD);
  return GRAD;
}


void QMCObjectiveFunction::numerical_grad(Array1D<double> &Params,
					  Array1D<double> &GRAD)
{
  //FIX this epsilon!!!!!!!! back to 0.001 or something 
  double epsilon = 0.001;
  int dimension  = Params.dim1();

  // Create a set of parameters displaced in every direction.
  // The last vector is the original set of Params.
  Array1D< Array1D<double> > Params_Vector(dimension+1);

  for(int i=0; i<dimension+1; i++)
    {
      Params_Vector(i) = Params;
    }

  for(int i=0; i<dimension; i++)
    {
      Params_Vector(i)(i) += epsilon;
    }

  Array1D< QMCObjectiveFunctionResult > Results = evaluate(Params_Vector);

  GRAD.allocate(dimension);
  for(int i=0; i<dimension; i++)
    {
      GRAD(i) = (Results(i).getDerivativeScore()-
		 Results(dimension).getDerivativeScore())/epsilon;
    }

}










