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

#include "QMCHarmonicOscillator.h"

QMCHarmonicOscillator::QMCHarmonicOscillator(QMCInput *input)
{
  Input = input;

  //This defines the eigenfunction
  w = 1.0;
  a = w/2.0;

  //Now we change the wavefunction a little bit
  a *= 1.01;

  cout << "Using Harmonic Oscilator w = " << w << endl;
  cout << "and                      a = " << a << endl;

  //try to be a little bit efficient
  a2     = a*a;
  w2div2 = 0.5*w*w;
}

void QMCHarmonicOscillator::evaluate(Array1D<QMCWalkerData *> &walkerData, Array1D<Array2D<double> * > &xData, int num, bool writeConfig)
{
  for(int i=0; i<num; i++)  
    {
      double psi  = 0;
      double ke   = 0;
      double pe   = 0;
      double grad = 0;

      for(int e=0; e<xData(i)->dim1(); e++)
	{
	  for(int d=0; d<xData(i)->dim2(); d++)
	    {
	      double x = (*xData(i))(e,d);
	      double x2 = x*x;

	      psi += - a * x2;

	      // = -1/2 (4 a^2 x^2 - 2 a)
	      ke += a - 2.0*a2*x2;
	      pe += w2div2*x2;

	      grad = -2.0 * x * a;

	      walkerData(i)->gradPsiRatio(e,d) = grad;
	    }
	}
      
      walkerData(i)->psi = exp(psi);
      walkerData(i)->neEnergy = 0;
      walkerData(i)->eeEnergy = 0;
      walkerData(i)->kineticEnergy = ke;
      walkerData(i)->potentialEnergy = pe;
      walkerData(i)->localEnergy = ke + pe;
      walkerData(i)->isSingular = false;
      walkerData(i)->modifiedGradPsiRatio = walkerData(i)->gradPsiRatio;
      
      if(false)
	{
	  printf("Psi %20.10e ", psi );
	  printf("KE %20.10e ", (*walkerData(i)).kineticEnergy );
	  printf("PE %20.10e ", (*walkerData(i)).potentialEnergy );
	  //printf("EE %20.10e ", (*walkerData(i)).eeEnergy );
	  printf("TE %20.10e\n", (*walkerData(i)).localEnergy );
	  cout << walkerData(i)->modifiedGradPsiRatio << endl;
	}
    }
}





