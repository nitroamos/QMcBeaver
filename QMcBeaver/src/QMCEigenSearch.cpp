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
#include <sstream>
#include "CubicSpline.h"
#include "IeeeMath.h"

bool QMCEigenSearch::currentlyHalf = true;
int QMCEigenSearch::orig_steps = 0;
Array1D<double> QMCEigenSearch::orig_params;
vector<double> QMCEigenSearch::adiag_tests;
Array2D<double> QMCEigenSearch::hamiltonian;
Array2D<double> QMCEigenSearch::overlap;

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

double QMCEigenSearch::get_a_diag(QMCDerivativeProperties & dp, double a_diag_factor)
{
  static double a_diag = globalInput.flags.a_diag;      

  if(globalInput.flags.a_diag > 0)
    {
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
      
      double cutoff = 1.0e-10;
      if( fabs(a_diag) < cutoff )
	a_diag = cutoff;
      
      return a_diag * a_diag_factor;
    }

  Array1D<double> samplesE = dp.getCorrelatedSamples(0);
  Array1D<double> samplesV = dp.getCorrelatedSamples(1);

  if(samplesV.dim1() == 0)
    {
      //clog << "Error: we were supposed to have been collecting correlated samples!" << endl;
      return a_diag;
    }

  Array1D<double> xvals(adiag_tests.size());
  Array1D<double> yvals(adiag_tests.size());
  xvals = 0.0;
  yvals = 0.0;

  double eng_frac = 0.95;
  double var_frac = 1.0 - eng_frac;
  printf("Correlated Sampling Results, Energy = %g\% Variance = %g\%:\n",100*eng_frac,100*var_frac);
  printf("       a_diag = %15s Energy %20.10g Sample Variance %20.10e\n",
	 "guide", samplesE(0), samplesV(0));

  double best = -99;
  int best_idx = -1;
  for(int cs=0; cs<xvals.dim1(); cs++)
    {
      int si = cs+1;
      double objv = (samplesV(si) - samplesV(0))/samplesV(0);
      double obje = (samplesE(si) - samplesE(0));
      double obj = eng_frac * obje + var_frac * objv;
      
      xvals(cs) = log(adiag_tests[cs]);
      yvals(cs) = obj;
      
      if(!IeeeMath::isNaN(obj))
	if(obj < best || best_idx < 0)
	  {
	    best_idx = cs;
	    best     = obj;
	    a_diag   = adiag_tests[cs];
	  }
      
      printf("CS %3i a_diag = %15.5g Energy %20.10g Sample Variance %20.10e objE %15.5e objV %15.5e obj %15.5e\n",
	     cs, adiag_tests[cs], samplesE(si), samplesV(si), obje, objv, obj);
    }

  cout << "Best a_diag = " << a_diag << " with objective value " << best << " idx = " << best_idx << endl;
  return a_diag;

  CubicSpline findmin;
  findmin.initializeWithFunctionValues(xvals,yvals,0,0);

  cout << "CubicSpline data:" << endl;
  for(int i=0; i<xvals.dim1(); i++)
    printf("%20.10e  %20.10e\n",xvals(i),yvals(i));

  /*
    I'm having trouble identifying the global minimum. I think
    the problem might be with the cubic spline -- each
    segment has a local minima?
  */
  double min1, min2, min3;
  double fmin1, fmin2, fmin3;

  if(best_idx-1 >= 0)
    min1  = findmin.minimum(xvals(best_idx-1),xvals(best_idx));
  else
    min1  = xvals(best_idx);
  
  if(best_idx+1 < xvals.dim1())
    min2  = findmin.minimum(xvals(best_idx),xvals(best_idx+1));
  else
    min2  = xvals(best_idx);

  min3  = findmin.minimum(xvals(0),xvals(xvals.size()-1));

  fmin1 = findmin.function(min1);  
  fmin2 = findmin.function(min2);
  fmin3 = findmin.function(min3);

  double min, fmin;
  min = fmin1 < fmin2 ? min1 : min2;
  fmin = findmin.function(min);
  min = fmin < fmin3 ? min : min3;
  fmin = findmin.function(min);
  min = fmin < best ? min : xvals(best_idx);
  fmin = findmin.function(min);

  printf("min1 = %20.10e fmin1 = %20.10e\n",min1,fmin1);
  printf("min2 = %20.10e fmin2 = %20.10e\n",min2,fmin2);
  printf("min3 = %20.10e fmin3 = %20.10e\n",min3,fmin3);
  printf("min  = %20.10e fmin  = %20.10e\n",min,fmin);

  a_diag = exp(min);
  cout << "Best a_diag = " << a_diag << " with objective value " << fmin << endl;

  return a_diag;
}

Array1D<double> QMCEigenSearch::optimize(Array1D<double> & CurrentParams,
					 QMCDerivativeProperties & dp,
					 double a_diag_factor,
					 int step)  
{
  optStep = step;
  stepinfo.clear();
  stepinfo << "(Step = " << setw(3) << optStep << "):";
  stepinfo.precision(12);
  
  cout.setf(ios::scientific);
  cout << endl;
  
  dim = CurrentParams.dim1();
  
  double InitialValue = dp.getParameterValue();
    
  x.push_back(CurrentParams);
  f.push_back(InitialValue);  

  if(globalInput.flags.a_diag > 0)
    {
      currentlyHalf = false;
      hamiltonian = dp.getParameterHamiltonian();
      overlap = dp.getParameterOverlap();
      orig_params = CurrentParams;
      orig_steps = globalInput.flags.max_time_steps;
    }

  cout << "Beginning Generalized Eigenvector Optimization ";
  cout << (currentlyHalf ? "half" : "full");
  cout << " step " << optStep << " for " << dim
       << " parameters, max steps = " << maximumSteps << ", with: " << endl;
  cout << "        InitialValue = " << InitialValue << endl;
  globalInput.printAIParameters(cout,"Parameters",20,CurrentParams,false);
  Array1D<double> gradient = dp.getParameterGradient();
  globalInput.printAIParameters(cout,"Gradient",20,gradient,true);

  Array1D<double> params;
  if(currentlyHalf)
    {
      orig_params = CurrentParams;
      orig_steps = globalInput.flags.max_time_steps;
      hamiltonian = dp.getParameterHamiltonian();
      overlap = dp.getParameterOverlap();
      adiag_tests.clear();

      if(false)
	{
	  double fac = sqrt(10.0);
	  adiag_tests.push_back(0.0);
	  adiag_tests.push_back(1e-15);
	  adiag_tests.push_back(1e-12);
	  adiag_tests.push_back(1e-10);
	  adiag_tests.push_back(1.0e-8);
	  while(adiag_tests[adiag_tests.size()-1] < 1000)
	    adiag_tests.push_back(fac * adiag_tests[adiag_tests.size()-1]);
	} else {
	  adiag_tests.push_back(1.0e-9);
	  while(adiag_tests[adiag_tests.size()-1] < 1000)
	    adiag_tests.push_back(100.0 * adiag_tests[adiag_tests.size()-1]);	  
	}

      int numTests = adiag_tests.size();
      globalInput.cs_Parameters.allocate(adiag_tests.size()+1);
      for(int cs=0; cs<numTests; cs++)
	{
	  Array1D<double> aParams = getParameters(dp,adiag_tests[cs],false);
	  globalInput.cs_Parameters(cs+1) = aParams;

	  stepinfo << endl;
	  stringstream temp;
	  temp.setf(ios::scientific);
	  temp << "CS " << adiag_tests[cs];
	  //cout << endl;
	  //globalInput.printAIParameters(cout,temp.str(),20,aParams,true);
	}

      globalInput.cs_Parameters(0) = CurrentParams;

      //whichever wavefunction is going to do the guiding needs
      //to be at the 0th index
      params = globalInput.cs_Parameters(0);
      if(orig_steps < 10000)
	{
	  globalInput.flags.max_time_steps = min(2000,orig_steps);
	  //globalInput.flags.max_time_steps = 2000;
	}
      else
	{
	  globalInput.flags.max_time_steps = min(20000,(int)(0.2*orig_steps));
	}
      globalInput.flags.max_time_steps += globalInput.flags.equilibration_steps;
      globalInput.flags.calculate_Derivatives = 0;
    }
  else
    {      
      globalInput.cs_Parameters.deallocate();
      double a_diag = get_a_diag(dp,a_diag_factor);
      params = getParameters(dp,a_diag,true); 
      globalInput.flags.max_time_steps = orig_steps;
      globalInput.flags.calculate_Derivatives = 1;
    }

  currentlyHalf = !currentlyHalf;
  cout << endl << "Ending Generalized Eigenvector Search Optimization... " << endl << endl;
  cout << stepinfo.str() << endl;
  return params;
}

Array1D<double> QMCEigenSearch::getParameters(QMCDerivativeProperties & dp, double a_diag, bool verbose)
{
  Array2D<double> ham  = hamiltonian;
  Array2D<double> olap = overlap;

  if( fabs(a_diag) > 0.0)
    {
      for(int d=1; d<dim+1; d++)
	ham(d,d) = ham(d,d) + a_diag;
      
      stepinfo << " a_diag = " << setw(20) << setprecision(10) << a_diag;
    }

  int largestPrintableMatrix = 10;
  if(verbose)
    {
      if(ham.dim1() <= largestPrintableMatrix)
	{
	  cout << "Hamiltonian:" << endl;
	  ham.printArray2D(cout,12,7,-1,',',true);
	  cout << "Overlap:\n" << endl;
	  olap.printArray2D(cout,12,7,-1,',',true);
	} else if(ham.dim1() < 20) {
	  cout << "Hamiltonian:" << endl;
	  ham.printArray2D(cout,2,7,-1,',',true);
	  cout << "Overlap:\n" << endl;
	  olap.printArray2D(cout,2,7,-1,',',true);
	}
    }

  Array2D<double> eigvec;
  Array1D<Complex> eigval;
  eigvec.generalizedEigenvectors(ham,
				 olap,
				 eigval);

  if(verbose)
    {
      if(eigvec.dim1() <= largestPrintableMatrix)
	{
	  cout << "Eigenvectors:" << endl;
	  eigvec.printArray2D(cout,12,7,-1,',',true);
	} else if(ham.dim1() < 20) {
	  cout << "Eigenvectors:" << endl;
	  eigvec.printArray2D(cout,2,7,-1,' ',true);
	}
    }

  int best_idx = 0;
  double best_val = 0.0;
  cout.precision(12);
  for(int i=0; i<dim+1; i++)
    {
      //cout << "Eigenvalue(" << setw(3) << i << "): " << eigval(i) << endl;
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
  stepinfo << " lowest eigenvalue = " << setw(20) << best_val;
  
  Array1D<double> delta_x(dim);
  for(int i=0; i<delta_x.dim1(); i++)
    delta_x(i) = eigvec[best_idx][i+1];

  delta_x /= eigvec[best_idx][0];

  if(verbose)
    globalInput.printAIParameters(cout,"Eigenvector",20,delta_x,true); 
  
  /*
    Some methods have fancy ways of modifying the distance
    we move in the search direction vector.
  */
  double rescale = 1.0;  
  if(stepLengthAlg != 0)
    {
      Array2D<double> fresh_overlap = overlap;
      double ksi = globalInput.flags.ksi;
      stepinfo << " ksi = " << setw(5) << ksi;
      rescale = stepLengthAlg->stepLength(OF,
					  delta_x,
					  delta_x,
					  delta_x,
					  fresh_overlap,
					  ksi);
    }
  stepinfo << " rescaling factor = " << setw(20) << rescale;

  Array1D<double> x_new = orig_params;

  if(IeeeMath::isNaN(rescale) || rescale < 0.05)
    {
      cout << "Warning: invalid rescale = " << rescale << endl;
      x_new = 0.0;
      return x_new;
    } else {
      //cout << "Warning: ok rescale = " << rescale << endl;
    }
  // Calculate the next step
  delta_x *= rescale;



  if(verbose)
    globalInput.printAIParameters(cout,"Delta params",20,delta_x,true);
  
  for(int j=0; j<dim; j++)
    x_new(j) += delta_x(j);

  return x_new;
}

QMCObjectiveFunction * QMCEigenSearch::getObjectiveFunction()
{
  return OF;
}
