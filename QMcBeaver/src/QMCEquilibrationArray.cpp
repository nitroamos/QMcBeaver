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
    {
      Eq_Array[i].zeroOut();
    }

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
      Eq_Array[i].getProperties()->setCalcDensity
                                                (calcDensity, nbasisfunctions);
    }
}

void QMCEquilibrationArray::newSample(QMCProperties * timeStepProps,
                                      double totalWeights, int nWalkers)
{ 
  //  cout << "Adding sample " 
  //       << Eq_Array[0].getProperties()->energy.getNumberSamples() + 1
  //       << ".  There are " << decorr_objects << " Properties objects alive."
  //       << endl;

  // Adds the new sample to the active elements of the array.
  for (int i=0; i<decorr_objects; i++)
    {
      Eq_Array[i].getProperties()->energy.newSample
        (timeStepProps->energy.getAverage(), totalWeights);

      Eq_Array[i].getProperties()->kineticEnergy.newSample
	(timeStepProps->kineticEnergy.getAverage(), totalWeights);

      Eq_Array[i].getProperties()->potentialEnergy.newSample
        (timeStepProps->potentialEnergy.getAverage(), totalWeights);

      Eq_Array[i].getProperties()->acceptanceProbability.newSample
        (timeStepProps->acceptanceProbability.getAverage(), totalWeights);

      Eq_Array[i].getProperties()->distanceMovedAccepted.newSample
        (timeStepProps->distanceMovedAccepted.getAverage(), totalWeights);

      Eq_Array[i].getProperties()->distanceMovedTrial.newSample
	(timeStepProps->distanceMovedTrial.getAverage(), totalWeights);

      Eq_Array[i].getProperties()->logWeights.newSample
	(timeStepProps->logWeights.getAverage(), nWalkers);
      if (calc_density)
	for (int j=0; j<nBasisFunc; j++)
	  {
	    Eq_Array[i].getProperties()->chiDensity(j).newSample
	      (timeStepProps->chiDensity(j).getAverage(), totalWeights);
	  }
    }

  // Activates the ith element on the 2^ith step.
  long check = power(2,decorr_objects);

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
  //  cout << "Choosing Properties object." << endl;

  int obj = 0;
  while(true)
    {
      if (Eq_Array[obj].getProperties()->energy.getStandardDeviation() < 1000)
        {
	  break;
	}
      else
        {
	  obj++;
	  if (obj == decorr_objects) return 0;
	}
    }

  double value = Eq_Array[0].getProperties()->energy.getAverage() + 
    3*Eq_Array[0].getProperties()->energy.getSeriallyCorrelatedVariance()/
    sqrt(Eq_Array[0].getProperties()->energy.getNumberSamples()-1.0);
  obj = 0;

  for (int i=1; i<decorr_objects; i++)
    {
      double test = Eq_Array[i].getProperties()->energy.getAverage() +
	3*Eq_Array[i].getProperties()->energy.getSeriallyCorrelatedVariance()/
	sqrt(Eq_Array[i].getProperties()->energy.getNumberSamples()-1.0);
      if (IeeeMath::isNaN(test)) continue;
      if (test < value)
	{
	  value = test;
          obj = i;
	}
    }
  // cout << "Object " << obj << " has the lowest standard deviation." << endl;
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
    {
      Eq_Array[i].getPropagationStopwatch()->start();
    }
  for (int j=decorr_objects; j<EQ; j++)
    {
      Eq_Array[j].getEquilibrationStopwatch()->start();
    }
}

void QMCEquilibrationArray::stopTimers()
{
  for (int i=0; i<EQ; i++)
    {
      if (Eq_Array[i].getPropagationStopwatch()->isRunning())
        {
	  Eq_Array[i].getPropagationStopwatch()->stop();
	}
      else if (Eq_Array[i].getEquilibrationStopwatch()->isRunning())
	{
	  Eq_Array[i].getEquilibrationStopwatch()->stop();
	}
    }
}

void QMCEquilibrationArray::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCEquilibrationArray>" << endl;
  strm << "<decorr_objects>\n" << decorr_objects << "\n</decorr_objects>\n";
  for (int i=0; i<EQ; i++) Eq_Array[i].getProperties()->toXML(strm);
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

  for (int i=0; i<EQ; i++) 
    {
      Eq_Array[i].getProperties()->readXML(strm);
    }
  
  strm >> temp;

  for (int i=0; i<decorr_objects; i++)
    {
      long step = power(2,i);
      Eq_Array[i].setStartingStep(step);
    }
}

long QMCEquilibrationArray::power(int a, int b)
{
  long result = 1;
  for (int i=0; i<b; i++)
    {
      result *= a;
    }
  return result;
}


