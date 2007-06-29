// This will become the array of decorrelation objects.

#include "QMCEquilibrationArray.h"

QMCEquilibrationArray::QMCEquilibrationArray()
{
  calc_density = false;
  nBasisFunc   = 0;
  zeroOut();
}

void QMCEquilibrationArray::zeroOut()
{
  for (int i=0; i<EQ; i++)
    Eq_Array[i].zeroOut();

  decorr_objects = 1;
  Eq_Array[0].setStartingStep(1);
}

void QMCEquilibrationArray::setCalcDensity(bool calcDensity, 
					                   int nbasisfunctions)
{
  calc_density = calcDensity;
  nBasisFunc = nbasisfunctions;

  for (int i=0; i<EQ; i++)
    {
      assert(0);
      //Eq_Array[i].getProperties()->setCalcDensity(calcDensity, nbasisfunctions);
    }
}

void QMCEquilibrationArray::setCalcForces(bool calcForces, int dim1, int dim2)
{
  calc_forces = calcForces;
  for (int i=0; i<EQ; i++)
    {
      assert(0);
      /*
      Eq_Array[i].getProperties()->setCalcForces
	(calc_forces, dim1,dim2);
      */
    }
}

void QMCEquilibrationArray::newSample(QMCProperties * timeStepProps,
                                      double totalWeights, int nWalkers)
{ 
  // Adds the new sample to the active elements of the array.
  for (int i=0; i<decorr_objects; i++)
    Eq_Array[i].getProperties()->newSample(timeStepProps, totalWeights, 
					   nWalkers);

  // Activates the ith element on the 2^ith step.
  unsigned long check = power(2,decorr_objects);

  if (Eq_Array[0].getProperties()->energy.getNumberSamples() == check-1)
    {
      Eq_Array[decorr_objects].setStartingStep(check);
      decorr_objects++;
    }
}

QMCProperties * QMCEquilibrationArray::chooseDecorrObject()
{
  int index = getDecorrObjectIndex();
  return Eq_Array[index].getProperties();
}

// Gets the index of the element of the array with the lowest variance for the 
// energy.
int QMCEquilibrationArray::getDecorrObjectIndex()
{
  int obj = 0;
  while(true)
    {
      if (Eq_Array[obj].getProperties()->energy.getStandardDeviation() < 99)
	break;
      else
        {
	  obj++;
	  if (obj == decorr_objects) return 0;
	}
    }

  int nSamples_tot = Eq_Array[0].getProperties()->energy.getNumberSamples();
  int nSamples_element = -1;

  double value = Eq_Array[0].getProperties()->energy.getAverage() + 50*
    sqrt(Eq_Array[0].getProperties()->energy.getSeriallyCorrelatedVariance()/
    (nSamples_tot-1));

  obj = 0;

  for (int i=1; i<decorr_objects; i++)
    {
      nSamples_element = Eq_Array[i].getProperties()->energy.getNumberSamples();
      double test = Eq_Array[i].getProperties()->energy.getAverage() +
        50*nSamples_tot*sqrt(Eq_Array[i].getProperties()->energy.getSeriallyCorrelatedVariance()/(nSamples_element-1))/nSamples_element;

      if (IeeeMath::isNaN(test)) continue;
      if (test < value)
	{
	  value = test;
          obj = i;
	}
    }
  return obj;
}

Stopwatch * QMCEquilibrationArray::getPropagationStopwatch()
{
  int index = getDecorrObjectIndex();
  return Eq_Array[index].getPropagationStopwatch();
}

Stopwatch * QMCEquilibrationArray::getEquilibrationStopwatch()
{
  int index = getDecorrObjectIndex();
  return Eq_Array[index].getEquilibrationStopwatch();
}

void QMCEquilibrationArray::startTimers()
{
  for (int i=0; i<decorr_objects; i++)
    Eq_Array[i].getPropagationStopwatch()->start();
  for (int j=decorr_objects; j<EQ; j++)
    Eq_Array[j].getEquilibrationStopwatch()->start();
}

void QMCEquilibrationArray::stopTimers()
{
  for (int i=0; i<EQ; i++)
    {
      if (Eq_Array[i].getPropagationStopwatch()->isRunning())
	Eq_Array[i].getPropagationStopwatch()->stop();
      else if (Eq_Array[i].getEquilibrationStopwatch()->isRunning())
	Eq_Array[i].getEquilibrationStopwatch()->stop();
    }
}

void QMCEquilibrationArray::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCEquilibrationArray>" << endl;
  strm << "<decorr_objects>\n" << decorr_objects << "\n</decorr_objects>\n";
  for (int i=0; i<decorr_objects; i++) 
    Eq_Array[i].getProperties()->toXML(strm);
  strm << "</QMCEquilibrationArray>" << endl;
}

void QMCEquilibrationArray::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;
  strm >> temp;
  strm >> temp;

  decorr_objects = atoi(temp.c_str());

  strm >> temp;

  for (int i=0; i<decorr_objects; i++) 
    Eq_Array[i].getProperties()->readXML(strm);
  
  strm >> temp;

  for (int i=0; i<decorr_objects; i++)
    {
      long step = power(2,i);
      Eq_Array[i].setStartingStep(step);
    }

  for (int i=decorr_objects; i<EQ; i++)
    Eq_Array[i].getProperties()->zeroOut();
}

long QMCEquilibrationArray::power(int a, int b)
{
  long result = 1;
  for (int i=0; i<b; i++)
    result *= a;
  return result;
}


