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

double QMCLinearizeStepLength::stepLength(QMCObjectiveFunction *function, 
					  Array1D<double> & delta_x,
					  Array1D<double> & unused1,
					  Array1D<double> & unused2,
					  Array2D<double> & overlap,
					  double ksi)
{
  //see JCP 126 084102 (2007) for description

  Array1D<double> N(delta_x.dim1());
  double denom;

  for(int ai=0; ai<N.dim1(); ai++)
    {
      
      if(isLinear(ai))
	{
	  N(ai) = overlap(ai+1,0);
	} else {
	  
	  double numerator = 0;
	  for(int j=0; j<N.dim1(); j++)
	    if(!isLinear(j))
	      numerator += delta_x(j) * overlap(ai+1,j+1);
	  numerator *= (1.0 - ksi);

	  denom = 1.0;
	  for(int j=0; j<N.dim1(); j++)
	    for(int k=0; k<N.dim1(); k++)
	      if(!isLinear(j) && !isLinear(k))
		denom += delta_x(j) * delta_x(k) * overlap(j+1,k+1);
	  denom = sqrt(denom);
	  denom = (1.0 - ksi) + ksi*denom;
	  
	  N(ai) = - numerator / denom;
	}
    }

  //it appears to be the same for all the ai
  cout << "Rescaling denominator: " << denom << endl;

  globalInput.printAIParameters(cout,"Rescaling terms",20,N,!true);
  
  double sum = 0;
  for(int i=0; i<N.dim1(); i++)
    sum += delta_x(i) * N(i);
  sum = 1.0 - sum;
  return 1.0 / sum;
}