// This class includes the data, and a timer and a counter about when the data
// was collected.

#include "QMCExtendedProperties.h"

QMCExtendedProperties::QMCExtendedProperties()
{
  zeroOut();
}

void QMCExtendedProperties::zeroOut()
{
  properties.zeroOut();
  startingStep = 0;
  equilibrationStopwatch.reset();
  propagationStopwatch.reset();
}


Stopwatch * QMCExtendedProperties::getEquilibrationStopwatch()
{
  return &(this->equilibrationStopwatch);
}

Stopwatch * QMCExtendedProperties::getPropagationStopwatch()
{
  return &(this->propagationStopwatch);
}

QMCProperties * QMCExtendedProperties::getProperties()
{
  return &properties;
}

long QMCExtendedProperties::getStartingStep()
{
  return startingStep;
}

void QMCExtendedProperties::setStartingStep(long i)
{
  startingStep = i;
}



