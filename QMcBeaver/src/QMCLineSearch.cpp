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

#include "QMCLineSearch.h"
#include "QMCInput.h"
#include <iomanip>

int QMCLineSearch::optStep = -1;

QMCLineSearch::QMCLineSearch(QMCObjectiveFunction *function, 
			     QMCLineSearchStepLengthSelectionAlgorithm * stepAlg,
			     int MaxSteps, 
			     double tol)
{
  dim           = 0;
  OF            = function;
  stepLengthAlg = stepAlg;
  maximumSteps  = MaxSteps;
  epsilon       = tol;
}

double QMCLineSearch::stepLength(Array1D<double> & x,Array1D<double> & p,
				 Array1D<double> & g, double f)
{
  if(stepLengthAlg != 0)
    return stepLengthAlg->stepLength(OF,x,p,g,f); 
  return 1.0;
}

Array1D<double> QMCLineSearch::searchDirection()
{
  calculateHessian();

  // calculate the search direction, p_k
  
  Array1D<double> p_k(dim);
  
  Array2D<double> & invHess = inverseHessian.back();
  Array1D<double> & grad    = gradient.back();

  for(int i=0; i<dim; i++)
    {
      p_k(i) = 0.0;
      for(int j=0; j<dim; j++)
	{
	  p_k(i) -= invHess(i,j) * grad(j);
	}
    }
  return p_k;
}

Array1D<double> QMCLineSearch::optimize(Array1D<double> & InitialGuess,
					double InitialValue,
					Array1D<double> & InitialGradient,
					Array2D<double> & InitialHessian,
					double a_diag_factor)
					
{
  if(fabs(a_diag_factor - 1.0) < 1.0e-10)
    optStep++;

  cout.setf(ios::scientific);
  cout << endl;

  dim = InitialGuess.dim1();
  bool useInitialHess = InitialHessian.dim1() == dim && InitialHessian.dim2() == dim;
  bool useInitialGrad = InitialGradient.dim1() == dim;

  if( (!useInitialGrad && InitialGradient.dim1() > 0) ||
      (!useInitialHess && InitialHessian.dim1() > 0))
    {
      clog << "Warning: were you intending to use InitialGradient/Hessian?" << endl
	   << "    InitialGuess.dim1() == " << dim << endl
	   << " InitialGradient.dim1() == " << InitialGradient.dim1() << endl
	   << "  InitialHessian.dim1() == " << InitialHessian.dim1() << endl << endl;
    }

  double a_diag = 0.0;
  if(globalInput.flags.optimize_Psi_method == "automatic")
    if(optStep < 6)
      {
	/*
	  Rotate the hessian in the direction of steepest
	  descent, described in PRL 94, 150201 (2005)
	  supposed to make a_diag smaller as optimization progresses...
	  
	  If the hessian is the identity matrix, then we have
	  the steepest descent method. Linear convergence?
	  
	  If the hessian is the second derivatives, then
	  conjugate gradient. Quadratic convergence?
	  
	  I think the point is compensating for deficiencies in
	  the approximate hessian, as well as for nonlinearities
	  in the objective function.


	  Maybe we want to use a different factor for the different
	  parameters?
	*/
	a_diag = 0.01;
	
	//use Steepest_Descent for a couple iterations first
	//useInitialHess = false;
      } else if(optStep < 12)
	{
	  a_diag = 0.01;
	} else if(optStep < 23)
	  {
	    a_diag = 2.0e-7;
	  }

  if( fabs(a_diag) > 0.0){
    a_diag *= a_diag_factor;
    for(int d=0; d<dim; d++)
      InitialHessian(d,d) = InitialHessian(d,d) + a_diag;
    
    cout << "Notice: Adding a_diag = " << a_diag
	 << " to the diagonal elements of the hessian" << endl;

    /*
    int numJP = globalInput.JP.getNumberEEParameters() +
      globalInput.JP.getNumberNEParameters();
    int numCI = globalInput.WF.getNumberDeterminants() - 1;
    if(globalInput.flags.optimize_Jastrows == 0)
      {
	cout << "Notice: we are not optimizing Jastrows this optimization step." << endl;
	for(int ai=0; ai<numJP; ai++)
	  {
	    InitialHessian.setToIdentity(ai);
	    InitialGradient(ai) = 0.0;
	  }
      }
    if(globalInput.flags.optimize_CI == 0)
      {
	cout << "Notice: we are not optimizing CI coefficients this optimization step." << endl;
	for(int ai=0; ai<numCI; ai++)
	  {
	    InitialHessian.setToIdentity(ai+numJP);
	    InitialGradient(ai+numJP) = 0.0;
	  }
      }
    */
  }
  
  // Allocate the point used in the search and initialize its value.
  x.push_back(InitialGuess);

  cout << "Beginning Line Search Optimization step " << optStep << " for " << dim
       << " parameters, max steps = " << maximumSteps << ", with: " << endl;
  cout << "        InitialValue = " << InitialValue << endl;
  globalInput.printAIParameters(cout,"InitialGuess",20,InitialGuess,false);
  if(useInitialGrad)
    globalInput.printAIParameters(cout,"InitialGradient",20,InitialGradient,true);

  if(useInitialHess)
    {
      double nsym = InitialHessian.nonSymmetry();
      cout << "      InitialHessian   (non symmetry = " << nsym << ") =" << endl;
      if(InitialHessian.dim1() < 20)
	cout << InitialHessian;
      else
	cout << "<too large, not printed>" << endl;
    }

  bool ok = false;
  Array2D<double> InitialInverseHessian(dim,dim);
  if(useInitialHess)
    {
      double det = 0.0;
      InitialHessian.determinant_and_inverse(InitialInverseHessian,
					     det, &ok);

      //we're transposing because the determinant_and_inverse returns
      //the transpose of what we want, so we need to untranspose it.
      InitialInverseHessian.transpose();
      if(!ok)
	{
	  cerr << "Error: InitialHessian can't be inverted!";
	  cerr << "   det = " << det << endl;
	  cerr << "   Using identity matrix instead." << endl;
	}
    }

  if(!ok || !useInitialHess)
    InitialInverseHessian.setToIdentity();
  inverseHessian.push_back(InitialInverseHessian);

  if(useInitialGrad)
    gradient.push_back(InitialGradient);

  cout << endl << endl;
  for(int i=0; i<maximumSteps; i++)
    {
      cout << endl << "\tIteration: " << i << endl;

      // Calculate the function value, step length, and search direction.
      //we may or may not have been given an inital value
      if(i == 0 && fabs(InitialValue - 99.0) > 1e-10)
	f.push_back(InitialValue);
      else
	f.push_back(OF->evaluate(x[i]).getScore());
      
      // if converged quit
      if( i>0 && fabs(1.0-f[i]/f[i-1]) < epsilon ) 
	{
	  cout << "Line Search Optimization Has Converged in " 
	       << i << " Iterations... " << endl; 
	  break;
	}

      if(i > 0 || !useInitialGrad)
	gradient.push_back(getObjectiveFunction()->grad(x[i]));

      if(i > 0)
	{
	  Array2D<double> new_inverseHessian(dim,dim);
	  new_inverseHessian.setToIdentity();
	  inverseHessian.push_back(new_inverseHessian);
	}

      /*
	this is where the hessian is determined
	and then A^-1 * b is calculated
      */
      Array1D<double> p_k   = searchDirection();

      /*
	Some methods have fancy ways of modifying the distance
	we move in the search direction vector.
      */
      double alpha_k = stepLength(x[i],
				  p_k,
				  gradient[i],
				  f[i]);

      cout << setw(20) << "Objective Value:";
      cout.precision(12);
      cout.width(20);
      cout << f[i] << endl;
      globalInput.printAIParameters(cout,"Parameters:",20,x[i],false);
      globalInput.printAIParameters(cout,"Gradient:",20,gradient[i],false);
      globalInput.printAIParameters(cout,"Search Direction:",20,p_k,false);
      cout << setw(20) << "stepLength:";
      cout.precision(12);
      cout.width(20);
      cout << alpha_k << endl;

      if(inverseHessian.back().isIdentity())
	{
	  cout << "\t\tInverseHessian:   (identity matrix)\n";
	} else {
	  double nsym = inverseHessian.back().nonSymmetry();
	  cout << "\t\tInverseHessian    (non symmetry = " << nsym << "):" << endl;
	  if(inverseHessian.back().dim1() < 20)
	    cout << inverseHessian.back();
	  else
	    cout << "<too large, not printed>" << endl;
	}

      // Calculate the next step
      Array1D<double> x_new = x.back();
      for(int j=0; j<dim; j++)
	{
	  x_new(j) += alpha_k * p_k(j);
	}
      x.push_back(x_new);
    }

  cout << endl << "Ending Line Search Optimization... " << endl;

  cout << endl;
  cout << setw(20) << "Best Objective Value (step = " << optStep << "):";
  cout.precision(12);
  cout.width(20);
  cout << f.back() << endl;
  globalInput.printAIParameters(cout,"Best Parameters:",20,x.back(),false);
  cout << endl << endl;

  /*
  for(unsigned int i=0; i<x.size(); i++)
    cout << "             x[" << i << "] = " << x[i];
  for(unsigned int i=0; i<gradient.size(); i++)
    cout << "      gradient[" << i << "] = " << gradient[i];
  for(unsigned int i=0; i<inverseHessian.size(); i++)
    cout << "inverseHessian[" << i << "] =\n" << inverseHessian[i];
  cout << endl;
  //*/

  return x.back();
}

QMCObjectiveFunction * QMCLineSearch::getObjectiveFunction()
{
  return OF;
}
