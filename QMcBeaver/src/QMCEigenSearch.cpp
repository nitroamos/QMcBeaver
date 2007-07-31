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

#include "QMCEigenSearch.h"
#include "QMCInput.h"
#include <iomanip>

QMCEigenSearch::QMCEigenSearch(QMCObjectiveFunction *function, 
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

Array1D<double> QMCEigenSearch::optimize(Array1D<double> & CurrentParams,
					 QMCDerivativeProperties & dp,
					 double a_diag_factor,
					 int optStep)  
{
  cout.setf(ios::scientific);
  cout << endl;
  
  dim = CurrentParams.dim1();
  
  double InitialValue = dp.getParameterValue();
  Array1D<double> gradient = dp.getParameterGradient();
  Array2D<double> hamiltonian = dp.getParameterHamiltonian();
  Array2D<double> overlap = dp.getParameterOverlap();
    
  x.push_back(CurrentParams);
  f.push_back(InitialValue);  
  
  cout << "Beginning Generalized Eigenvalue Optimization step " << optStep << " for " << dim
       << " parameters, max steps = " << maximumSteps << ", with: " << endl;
  cout << "        InitialValue = " << InitialValue << endl;
  globalInput.printAIParameters(cout,"Parameters",20,CurrentParams,false);
  globalInput.printAIParameters(cout,"Gradient",20,gradient,true);
  
  static double a_diag = globalInput.flags.a_diag;      
  if(globalInput.flags.optimize_Psi_method == "automatic" &&
     fabs(a_diag_factor - 1.0) < 1e-10)
    if(optStep == 6)
      {
	a_diag *= 1e-1;
      } else if(optStep == 12)
	{
	  a_diag *= 1e-2;
	} else if(optStep == 23)
	  {
	    a_diag *= 1e-5;
	  }
  
  if( fabs(a_diag) > 0.0)
    {
      for(int d=1; d<dim+1; d++)
	hamiltonian(d,d) = hamiltonian(d,d) + a_diag * a_diag_factor;
      
      cout << "Notice: Adding a_diag = " << (a_diag * a_diag_factor)
	   << " to the diagonal elements of the hamiltonian" << endl;
    }

  cout << endl << endl;

  if(hamiltonian.dim1() < 20)
    {
      cout << "Hamiltonian:\n" << hamiltonian << endl;
      cout << "Overlap:\n" << overlap << endl;
    }

  Array2D<double> eigvec;
  Array1D<Complex> eigval;
  eigvec.generalizedEigenvectors(hamiltonian,
				 overlap,
				 eigval);

  if(eigvec.dim1() < 20)
    cout << "Eigenvectors:\n" << eigvec << endl;

  int best_idx = 0;
  double best_val = 0.0;
  for(int i=0; i<dim+1; i++)
    {
      cout << "Eigenvalue(" << i << "): " << eigval(i) << endl;
      double val = eigval(i).real();
      if( fabs(eigval(i).imaginary()) < 1e-50 &&
	  val < best_val &&
	  !IeeeMath::isNaN(val) &&
	  fabs(val) < 1.0e5)
	{
	  best_val = eigval(i).real();
	  best_idx = i;
	}
    }
  cout << "Lowest eigenvalue: " << best_val << endl;
  
  Array1D<double> delta_x(dim);
  for(int i=0; i<delta_x.dim1(); i++)
    delta_x(i) = eigvec[best_idx][i+1];

  delta_x /= eigvec[best_idx][0];
  globalInput.printAIParameters(cout,"Eigenvector",20,delta_x,true); 
  
  /*
    Some methods have fancy ways of modifying the distance
    we move in the search direction vector.
  */
  double rescale = 1.0;  
  if(stepLengthAlg != 0)
    {
      Array2D<double> fresh_overlap = dp.getParameterOverlap();
      rescale = stepLengthAlg->stepLength(OF,
					  delta_x,
					  delta_x,
					  delta_x,
					  fresh_overlap,
					  0.5);
    }

  cout << setw(20) << "rescaling factor:";
  cout.precision(12);
  cout.width(20);
  cout << rescale << endl;
    
  // Calculate the next step
  Array1D<double> x_new = x.back();
  for(int j=0; j<dim; j++)
    x_new(j) += rescale * delta_x(j);

  cout << endl << "Ending Generalized Eigenvalue Search Optimization... " << endl;

  globalInput.printAIParameters(cout,"Delta params",20,delta_x,true);   
  return x_new;
}

QMCObjectiveFunction * QMCEigenSearch::getObjectiveFunction()
{
  return OF;
}
