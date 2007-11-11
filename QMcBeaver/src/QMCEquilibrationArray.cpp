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

  // I want to eliminate elements 1, 2, 3, 4, 5, and 6 from the array.  The
  // reason is that element 7 has 128 equilibration steps, so the previous
  // contain mostly redundant data. 

  // Activates the ith element on the 2^ith step.
  unsigned long check = power(2,decorr_objects+6);

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

  // If there are less than 256 samples, return object 0.
  if (decorr_objects < 3) return 0;

  // Here we look for the plateau in the objective function.
  // To decide which objective function to use, we compare the expectation 
  // value with 128 equilibration steps to the expectation value with none.
  double sgn = Eq_Array[1].getProperties()->energy.getAverage() - 
    Eq_Array[0].getProperties()->energy.getAverage();

  if (sgn < 0.0)
    {
      // If the expectation value decreases as the number of equilibration 
      // steps increases, we add the standard deviation to the expectation 
      // value for the objective function.

      double value = Eq_Array[0].getProperties()->energy.getAverage() +
        Eq_Array[0].getProperties()->energy.getBlockStandardDeviation(0);

      // The objective function will generally be decreasing, so the plateau 
      // occurs when the value of object i is greater than the value of object 
      // i-1.

      for (int i=1; i<decorr_objects; i++)
        {
          double test = Eq_Array[i].getProperties()->energy.getAverage() + 
            Eq_Array[i].getProperties()->energy.getBlockStandardDeviation(0);

          if (IeeeMath::isNaN(test)) continue;
          if (test > value) 
            return i-1;
          else
            value = test;
        }
    }

  else if (sgn > 0.0)
    {
      // If the expectation value increases as the number of equilibration
      // steps increases, we subtract the standard deviation from the 
      // expectation value for the objective function.
  
      double value = Eq_Array[0].getProperties()->energy.getAverage() -
        Eq_Array[0].getProperties()->energy.getBlockStandardDeviation(0);

      for (int i=1; i<decorr_objects; i++)
        {
          double test = Eq_Array[i].getProperties()->energy.getAverage() -
            Eq_Array[i].getProperties()->energy.getBlockStandardDeviation(0);

          if (IeeeMath::isNaN(test)) continue;
          if (test < value)
            return i-1;
          else
            value = test;
        }
    }

  // If we get to this point, the energy is changing monotonically through our
  // array and has not equilibrated according to our algorithm.  We return
  // the last element that has at least half as many samples as the one before
  // it.  If the last element has a small number of samples and a crazy value, 
  // this will prevent it from being chosen.

  unsigned long nSamples_last = Eq_Array[decorr_objects-1].getProperties()->energy.getNumberSamples();
  unsigned long nSamples_prev = Eq_Array[decorr_objects-2].getProperties()->energy.getNumberSamples();

  if (nSamples_last > nSamples_prev/2.0)
    return decorr_objects-1;
  else
    return decorr_objects-2;
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

bool QMCEquilibrationArray::readXML(istream& strm)
{
  string temp;

  // Open XML
  strm >> temp;
  if (temp != "<QMCEquilibrationArray>")
    return false;

  strm >> temp;
  if (temp != "<decorr_objects>")
    return false;
  strm >> temp;
  decorr_objects = atoi(temp.c_str());
  strm >> temp;
  if (temp != "</decorr_objects>")
    return false;

  for (int i=0; i<decorr_objects; i++) 
    if (!Eq_Array[i].getProperties()->readXML(strm))
      return false;
  
  strm >> temp;
  if (temp != "</QMCEquilibrationArray>")
    return false;

  Eq_Array[0].setStartingStep(1);

  for (int i=1; i<decorr_objects; i++)
    {
      long step = power(2,i+6);
      Eq_Array[i].setStartingStep(step);
    }

  for (int i=decorr_objects; i<EQ; i++)
    Eq_Array[i].getProperties()->zeroOut();
  
  return true;
}

long QMCEquilibrationArray::power(int a, int b)
{
  long result = 1;
  for (int i=0; i<b; i++)
    result *= a;
  return result;
}


