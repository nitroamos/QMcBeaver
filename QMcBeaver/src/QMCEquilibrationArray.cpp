// This will become the array of decorrelation objects.

#include "QMCEquilibrationArray.h"

QMCEquilibrationArray::QMCEquilibrationArray()
{
  over_ten_thousand = 0;
  calc_density = false;
  nBasisFunc   = 0;

  ten_thousand_samples.allocate(EQ,2);
  zeroOut();
}

void QMCEquilibrationArray::zeroOut()
{
  over_ten_thousand = 0;
  ten_thousand_samples = 0.0;
  for (int i=0; i<EQ; i++)
    Eq_Array[i].zeroOut();

  ZeroProperties.zeroOut();
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

  if(Eq_Array[over_ten_thousand].getProperties()->energy.getNumberSamples() == 10000)
    {
      ten_thousand_samples(over_ten_thousand,0) = Eq_Array[over_ten_thousand].getProperties()->energy.getAverage();
      ten_thousand_samples(over_ten_thousand,1) = Eq_Array[over_ten_thousand].getProperties()->energy.getBlockStandardDeviation(0);
      over_ten_thousand++;
    }

  // Activates the ith element on the 2^ith step.

  // I want to eliminate elements 1, 2, 3, 4, 5, and 6 from the array.  The
  // reason is that element 7 has 128 equilibration steps, so the previous
  // contain mostly redundant data. 

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
  if (index == -1)
    return ZeroProperties.getProperties();
  else
    return Eq_Array[index].getProperties();
}

// Gets the index of the element of the array with the lowest variance for the 
// energy.
int QMCEquilibrationArray::getDecorrObjectIndex()
{
  // We are basing equilibration on the first ten thousand samples.  We need
  // the first element to have at least twenty thousand samples before we can
  // tell if it is equilibrated.  Before then it will return zero samples.

  if (Eq_Array[0].getProperties()->energy.getNumberSamples() < 20000)
    return -1;

  double init_hi, init_lo, final_hi, final_lo, final_avg;
  int obj = -1;

  for (int i=0; i<over_ten_thousand; i++)
    {
      init_hi = ten_thousand_samples(i,0) + 2*ten_thousand_samples(i,1);
      init_lo = ten_thousand_samples(i,0) - 2*ten_thousand_samples(i,1);

      final_hi = Eq_Array[i].getProperties()->energy.getAverage() + 2*Eq_Array[i].getProperties()->energy.getBlockStandardDeviation(0);
      final_lo = Eq_Array[i].getProperties()->energy.getAverage() - 2*Eq_Array[i].getProperties()->energy.getBlockStandardDeviation(0);
      final_avg = Eq_Array[i].getProperties()->energy.getAverage();

      if ((final_hi > init_lo && final_hi < init_hi) || (final_lo > init_lo && final_lo < init_hi) || (final_avg > init_lo && final_avg < init_hi))
	{
	  obj = i;
	  break;
	}
    }

  if ((obj == -1) || (Eq_Array[obj].getProperties()->energy.getNumberSamples() < 20000))
    return -1;
  else
    return obj;
}

Stopwatch * QMCEquilibrationArray::getPropagationStopwatch()
{
  int index = getDecorrObjectIndex();
  if (index == -1)
    return ZeroProperties.getPropagationStopwatch();
  else
    return Eq_Array[index].getPropagationStopwatch();
}

Stopwatch * QMCEquilibrationArray::getEquilibrationStopwatch()
{
  int index = getDecorrObjectIndex();
  if (index == -1)
    return ZeroProperties.getEquilibrationStopwatch();
  else
    return Eq_Array[index].getEquilibrationStopwatch();
}

void QMCEquilibrationArray::startTimers()
{
  for (int i=0; i<decorr_objects; i++)
    Eq_Array[i].getPropagationStopwatch()->start();
  for (int j=decorr_objects; j<EQ; j++)
    Eq_Array[j].getEquilibrationStopwatch()->start();
  ZeroProperties.getEquilibrationStopwatch()->start();
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
  if (ZeroProperties.getEquilibrationStopwatch()->isRunning())
    ZeroProperties.getEquilibrationStopwatch()->stop();
}

void QMCEquilibrationArray::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCEquilibrationArray>" << endl;
  strm << "<decorr_objects>\n" << decorr_objects << "\n</decorr_objects>\n";
  strm << "<over_ten_thousand>\n" << over_ten_thousand << "\n</over_ten_thousand>\n";
  strm << "<ten_thousand_samples>\n";
  for (int i=0; i<over_ten_thousand; i++)
    strm << ten_thousand_samples(i,0) << "\t" << ten_thousand_samples(i,1) << "\n";
  strm << "</ten_thousand_samples>\n";
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

  strm >> temp;
  if (temp != "<over_ten_thousand>")
    return false;
  strm >> temp;
  over_ten_thousand = atoi(temp.c_str());
  strm >> temp;
  if (temp != "</over_ten_thousand>")
    return false;

  strm >> temp;
  if (temp != "<ten_thousand_samples>")
    return false;
  for (int i=0; i<over_ten_thousand; i++)
    {
      strm >> temp;
      ten_thousand_samples(i,0) = atof(temp.c_str());
      strm >> temp;
      ten_thousand_samples(i,1) = atof(temp.c_str());
    }
  strm >> temp;
  if (temp != "</ten_thousand_samples>")
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
  
  ZeroProperties.zeroOut();

  return true;
}

long QMCEquilibrationArray::power(int a, int b)
{
  long result = 1;
  for (int i=0; i<b; i++)
    result *= a;
  return result;
}


