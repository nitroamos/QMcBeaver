// This will become the array of decorrelation objects.

#include "QMCEquilibrationArray.h"

QMCEquilibrationArray::QMCEquilibrationArray()
{
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

void QMCEquilibrationArray::newSample(QMCProperties * timeStepProps,
                                      double totalWeights, int nWalkers)
{ 
  //  cout << "Adding sample " 
  //       << Eq_Array[0].getProperties()->energy.getNumberSamples() + 1
  //       << ".  There are " << decorr_objects << " Properties objects alive."
  //       << endl;


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
    }

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

  double value = Eq_Array[0].getProperties()->energy.getStandardDeviation() + 
   Eq_Array[0].getProperties()->energy.getStandardDeviationStandardDeviation();
  obj = 0;

  for (int i=1; i<decorr_objects; i++)
    {
      double test = Eq_Array[i].getProperties()->energy.getStandardDeviation()
 + Eq_Array[i].getProperties()->energy.getStandardDeviationStandardDeviation();

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


