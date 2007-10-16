#include <iostream>
#include "QMCLinearizeStepLength.h"
#include "QMCInput.h"
#include <iomanip>

using namespace std;

  
bool QMCLinearizeStepLength::isLinear(int ai)
{
  int numJW = globalInput.JP.getNumberJWParameters();
  int numCI = globalInput.WF.getNumberCIParameters();
  int numOR = globalInput.WF.getNumberORParameters();

  if(ai >= numJW && ai < (numJW + numCI))
    return true;
  return false;
}

double QMCLinearizeStepLength::rescalingJCP()
{
  //see JCP 126 084102 (2007) eq 34 for description
  Array1D<double> N(dp.dim1());
  double denom = 1.0;

  for(int j=0; j<N.dim1(); j++)
    for(int k=0; k<N.dim1(); k++)
      if(!isLinear(j) && !isLinear(k))
	denom += dp(j) * dp(k) * S(j+1,k+1);
  if(verbose) printf("den = %20.10e",denom);
  denom = sqrt(denom);
  denom = (1.0 - ksi) + ksi*denom;
  if(verbose) printf(" -> den = %20.10e\n",denom);
  if(verbose) cout << "Rescaling denominator: " << denom << endl;
  for(int ai=0; ai<N.dim1(); ai++)
    {      
      if(isLinear(ai))
	{
	  N(ai) = S(ai+1,0);
	} else {
	  
	  double numerator = 0;
	  for(int j=0; j<N.dim1(); j++)
	    if(!isLinear(j))
	      numerator += dp(j) * S(ai+1,j+1);
	  if(verbose) printf("ai %3i num = %20.10e",ai,numerator);
	  numerator *= (1.0 - ksi);
	  
	  N(ai) = - numerator / denom;
	  if(verbose) printf(" Nai = %20.10e\n",N(ai));
	}
    }

  if(verbose) globalInput.printAIParameters(cout,"JCP terms",20,N,!true);

  double sum = 0;
  for(int i=0; i<N.dim1(); i++)
    sum += dp(i) * N(i);
  sum = 1.0 - sum;

  return 1.0 / sum;
}

double QMCLinearizeStepLength::rescalingPRL()
{
  Array1D<double> N(dp.dim1());

  double sum_S0jPj = 0;
  for(int j=0; j<N.dim1(); j++)
    if(!isLinear(j))
      sum_S0jPj += S(0,j+1) * dp(j);

  double sum_SjkPjPk = 0;
  for(int j=0; j<N.dim1(); j++)
    for(int k=0; k<N.dim1(); k++)
      if(!isLinear(j) && !isLinear(k))
	sum_SjkPjPk += S(j+1,k+1) * dp(j) * dp(k);
  
  double D     = sqrt(1.0 + 2.0 * sum_S0jPj + sum_SjkPjPk);
  double denom = ksi * D + (1.0 - ksi)*(1.0 + sum_S0jPj);

  if(verbose) printf("S0jPj = %20.10e SjkPjPk = %20.10e D = %20.10e denom = %20.10e\n",
		     sum_S0jPj, sum_SjkPjPk, D, denom);

  for(int i=0; i<N.dim1(); i++)
    {      
      double sum_SijPj = 0;
      for(int j=0; j<N.dim1(); j++)
	if(!isLinear(j))
	  sum_SijPj += S(i+1,j+1) * dp(j);

      double numerator = ksi * D * S(0,i+1) + (1.0 - ksi)*(S(0,i+1) + sum_SijPj);

      N(i) = - numerator / denom;
      
      if(verbose) printf("%3i S0i = %20.10e SijPj = %20.10e num = %20.10e N = %20.10e\n",
			 i, S(0,i+1), sum_SijPj, numerator, N(i));
    }

  if(verbose) globalInput.printAIParameters(cout,"PRL terms",20,N,!true);
  double sum = 0;
  for(int i=0; i<N.dim1(); i++)
    sum += dp(i) * N(i);
  sum = 1.0 - sum;

  return 1.0 / sum;
}

double QMCLinearizeStepLength::stepLength(QMCObjectiveFunction *function, 
					  Array1D<double> & delta_x,
					  Array1D<double> & unused1,
					  Array1D<double> & unused2,
					  Array2D<double> & overlap,
					  double ksi)
{
  this->ksi = ksi;
  S  = overlap;
  dp = delta_x;
  /*
    ksi = 0 is supposed to be safer
    ksi = 1 is supposed to be the "SR method"
  */
  verbose = !true;

  if(ksi > 1.0 || ksi < 0.0)
    {
      cerr << "Warning: bad value for ksi in QMCLinearizeStepLength!" << endl;
      cerr << " ksi = " << ksi << endl;
      ksi = 0.5;
      cerr << "changing, ksi is now = " << ksi << endl;
    }

  double jcp_scale = rescalingJCP();
  double prl_scale = rescalingPRL();
  double scale = jcp_scale;
  //double scale = prl_scale;

  return scale;
}
