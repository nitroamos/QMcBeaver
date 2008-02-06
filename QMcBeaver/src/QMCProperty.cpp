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

#include "QMCProperty.h"

QMCProperty::QMCProperty()
{
  /*
   for(int i=0; i<numberParams; i++)
   {
     stdevFittingParameters[i] = 0.0;
   }
   */
  
  // initialize the rest of the calculation
  
  zeroOut();
  
#ifdef PARALLEL
  if( !mpiTypeCreated )
  {
    mpiTypeCreated = true;
    buildMpiType();
    buildMpiReduce();
  }
#endif
  
}

void QMCProperty::zeroOut()
{
  for(int i=0;i<DCL;i++)
  {
    DeCorr[i].zeroOut();
    //DeCorr_flags[i]  = 0;
    //DeCorr_sample[i] = 0.0;
    //DeCorr_weight[i] = 0.0;
  }
  memset(&DeCorr_flags,0,sizeof(int)*DCL);
  memset(&DeCorr_sample,0,sizeof(double)*DCL);
  memset(&DeCorr_weight,0,sizeof(double)*DCL);
}

unsigned long QMCProperty::getNumberSamples()
{
  return DeCorr[0].getNumberSamples();
}

double QMCProperty::getAverage()
{
  return DeCorr[0].getAverage();
}

double QMCProperty::getStandardDeviation()
{
  int decorr_depth = getDecorrDepth();
  if (decorr_depth == -2)
  {
    // return infinity
    return 1.0e300;
  }
  else if (decorr_depth == -1)
  {
    // DDDA has not converged.
    return 99;
  }
  else return getBlockStandardDeviation(decorr_depth);
}

double QMCProperty::getSkewness()
{
  int decorr_depth = getDecorrDepth();
  if (decorr_depth == -2)
  {
    // return infinity
    return 1.0e300;
  }
  else if (decorr_depth == -1)
  {
    // DDDA has not converged.
    return 99;
  }
  else return getBlockSkewness(decorr_depth);
}

double QMCProperty::getKurtosis()
{
  int decorr_depth = getDecorrDepth();
  if (decorr_depth == -2)
  {
    // return infinity
    return 1.0e300;
  }
  else if (decorr_depth == -1)
  {
    // DDDA has not converged.
    return 99;
  }
  else return getBlockKurtosis(decorr_depth);
}

int QMCProperty::getDecorrDepth()
{
  /*
   Fit the decorrelation standard deviation to 
   f(x;a) = a^2 * e^( c^2 * x + d^2 * x^2 )/(1 + b^2 * e^( c^2 * x + 
                                                           d^2 * x^2 ));
   x = log2(blocksize);
   The fit is done using the sum of the square of the residuals.
   
   100 samples must be accumulated before the fitting procedure is used.
   Otherwise, infinity is returned.  Before the fitting procedure is used,
   initial guess parameter sets are calculated from the accumulated data.
   After the fitting has begun, the previous optimal parameter set is the
   initial guess set of parameters.
   */
  
  // if too few data points have been accumulated, return an infinite
  // standard deviation
  
  if( getNumberSamples() < 100 )
  {
    return -2;
  }
  
  /******************** TEST **********************************/
  
  double value = getBlockStandardDeviation(0) - 
	getBlockStandardDeviationStandardDeviation(0);
  
  int block = 1;
  
  while(true)
  {
    if( DeCorr[block].getNumberSamples() <= 4 ) break;
    
    double temp = getBlockStandardDeviation(block) - 
      getBlockStandardDeviationStandardDeviation(block);
    
    if( temp <= value )
    {
      return block-1;          
    }
    else
    {
      value = temp;
    }
    
    block++;
  }
  
  return -1;
  
  /*********************************************************/
  
  
  /*
   
   // copy the parameters into an array to allow easier matrix-vector math
   
   Array1D<double> params(numberParams);
   
   for(int i=0;i<numberParams;i++)
   {
     params(i) = stdevFittingParameters[i];
   }
   
   // if the initial guess parameters have no length generate new ones
   
   generateInitialGuessFittingParameters();
   
   
   for(int i=0; i<numberParams; i++)
   {
     if( params(i)*params(i) < 1e-15 )
     {
       generateInitialGuessFittingParameters();
       
       // copy the new parameters into an array
       
       for(int i=0;i<numberParams;i++)
       {
         params(i) = stdevFittingParameters[i];
       }
     }
   }
   
   // find the number of blocks to use in the calculation
   
   int blocks = 0;
   
   while(true)
   {
     if( DeCorr[blocks].getNumberSamples() <= 4 ) break;
     blocks++;
   }
   
   
   
   // generate a vector of standard deviations and associated errors
   // this is a performance tuning adjustment
   
   Array1D<double> standardDeviations(blocks);
   Array1D<double> standardDeviationsErrors(blocks);
   
   for(int i=0; i<blocks; i++)
   {
     standardDeviations(i) = getBlockStandardDeviation(i);
     standardDeviationsErrors(i) = 
       getBlockStandardDeviationStandardDeviation(i);
   }
   
   
   // initialize the conjugate gradient optimization
   
   bool firstStep = true;
   
   Array1D<double> grad(numberParams);
   grad = 0.0;
   Array1D<double> gradOld(numberParams);
   gradOld = 0.0;
   Array1D<double> step(numberParams);
   step = 0.0;
   
   double functionValue = -1;
   double functionValueOld = -1;
   
   const double epsilon = 1e-8;
   
   
   // perform a Polak-Ribiere conjugate gradient optimization of the parameters
   
   while( firstStep || ( (grad*grad) > epsilon && fabs(functionValue - functionValueOld) > epsilon ) )
   //fabs(functionValue - functionValueOld) > epsilon ) 
   //(grad*grad) > epsilon )
   {
     //cout << functionValue << "\t" << functionValue - functionValueOld << "\t" << grad*grad << endl;
     //cout << "\t\tP: " << params << endl;
     //cout << "\t\t\tG: " << grad << endl;
     cout << ".";
     
     // calculate the gradient of the residual squared
     
     gradOld = grad;
     
     functionValueOld = functionValue;
     
     calculateObjectiveFunction(params, standardDeviations, 
                                standardDeviationsErrors, functionValue, 
                                grad);
     
     
     // if the grad comes across a NaN set the component to zero
     // this isn't rigorously correct, but the specific objective function
     // being optimized is plagued with overflow and underflow problems
     // from finite arithmatic.  This is an easy correction to most of these.
     
     for(int i=0; i<grad.dim1(); i++)
     {
       if( isnan( grad(i) ) != 0 )
       {
         grad(i) = 0.0;
       }
     }
     
     
     
     // calculate the next set of parameters using polak-ribiere 
     // conjugate gradient step
     
     double beta;
     
     if( firstStep )
     {
       beta = 0;
     }
     else
     {
       beta = (grad*grad - gradOld*grad)/(gradOld*gradOld);
     }
     
     for(int i=0; i<numberParams; i++)
     {
       step(i) = -grad(i) + beta*step(i);
     }
     
     
     // calculate the step length to use
     const double alpha_guess = 2.0;
     
     double alpha = wolfeStepLength(alpha_guess, params, 
                                    step, grad, functionValue, 
                                    standardDeviations, 
                                    standardDeviationsErrors);
     
     
     // calculate the new set of parameters
     
     for(int i=0; i<params.dim1(); i++)
     {
       params(i) = params(i) + 
       alpha * step(i);
     }
     
     firstStep = false;
   }
   
   for(int i=0;i<numberParams;i++)
   {
     stdevFittingParameters[i] = params(i);
   }
   
   cout << endl;
   
   return params(0)*params(0)/(params(1)*params(1));
   */
}

double QMCProperty::getSeriallyCorrelatedVariance()
{
  return DeCorr[0].getVariance();
}

double QMCProperty::getVariance()
{
  double std = getStandardDeviation();
  return std*std;
}

double QMCProperty::getSeriallyCorrelatedStandardDeviation()
{
  return sqrt(getSeriallyCorrelatedVariance());
}

double QMCProperty::getStandardDeviationStandardDeviation()
{
  int decorr_depth = getDecorrDepth();
  if (decorr_depth < 0) return 0.0;
  else return getBlockStandardDeviationStandardDeviation(decorr_depth);
}

void QMCProperty::newSample(double s, double weight)
{
  //This is a binary sum of the flags with the new sample bringing in a "1"
  //in the 2^0 spot. (This means that the binary representation in DeCorr_flags
  //starts at the DeCorr_flags[1] level
  
  int done=0;
  //always add the current sample to the "0"th DeCorr element
  DeCorr[0].newSample(s,weight);
  int inc=1;
  double pushed_sample=s; 
  double pushed_weight=weight; 
  
  while(!done)
  {
    if(inc>=DCL)
    {  //so we really don't care what happens with the pushed value
      done=1;
    }
    else
    {
      if(DeCorr_flags[inc]==0)//this means that this box is empty
	    {                     //and we won't need to push
                            //put the pushed up sample in the queueing box
	      DeCorr_sample[inc]=pushed_sample;
	      DeCorr_weight[inc]=pushed_weight;
	      DeCorr_flags[inc]=1; //it is now holding a sample
	      done=1;
	    }
      else // inc<DCL => this means that we need to push values up a box
	    {
	      //first add the pushed_sample to the old pushed_sample
	      //which will give us a new pushed_sample
        
	      pushed_sample=(pushed_sample*pushed_weight
                       +DeCorr_sample[inc]*DeCorr_weight[inc])
        /(pushed_weight+DeCorr_weight[inc]);
	      pushed_weight=(pushed_weight+DeCorr_weight[inc])/2.0;
	      
	      //then add the new pushed_sample to this DeCorr to be counted
	      DeCorr[inc].newSample(pushed_sample,pushed_weight);
	      
	      //this also means we must unflag this box
	      DeCorr_flags[inc]=0;  //it now holds nothing of value
	      
	      //now repeat pushing the new pushed_sample up another level
	      //(this is done keeping done=0 in the while loop)
	      inc++;
	    }
    }
  }
}

void QMCProperty::operator = ( const QMCProperty & rhs )
{
  memcpy(&DeCorr, &rhs.DeCorr, sizeof(QMCStatistic)*DCL);
  memcpy(&DeCorr_flags, &rhs.DeCorr_flags, sizeof(int)*DCL);
  memcpy(&DeCorr_sample, &rhs.DeCorr_sample, sizeof(double)*DCL);
  memcpy(&DeCorr_weight, &rhs.DeCorr_weight, sizeof(double)*DCL);
  
  // Don't include the following line because this data is not sent by MPI
  // Because of this, some MPI implementations cause segmentation faults.
  // stdevFittingParameters = rhs.stdevFittingParameters;
}

void QMCProperty::reWeight(double w)
{
  for(int i=0; i<DCL; i++)
  {
    DeCorr[i].reWeight(w);
    if(DeCorr_flags[i]==0)
    {
      DeCorr_weight[i] *= w;
    }
  }
}

void QMCProperty::operator += ( QMCProperty &rhs)
{
 *this = *this + rhs;
}

QMCProperty QMCProperty::operator + ( QMCProperty &rhs)
{
  QMCProperty result;
  
  for(int i=0;i<DCL;i++)
  {
    result.DeCorr[i]=DeCorr[i]+rhs.DeCorr[i];
  }
  
  int carry_up=0;
  double carry_up_sample=0.0;
  double carry_up_weight=0.0;
  int binary_sum=0;         //this is the binary sum for the DeCorr level i
  for(int i=0;i<DCL;i++)
  {
    binary_sum=carry_up+DeCorr_flags[i]+rhs.DeCorr_flags[i];
    if(binary_sum==0) //no new sample needs to get added in here
    {
      result.DeCorr_flags[i]=0;
      result.DeCorr_sample[i]=0.0;
      carry_up=0;
      carry_up_sample=0.0;
      carry_up_weight=0.0;
    }
    if(binary_sum==1) 
    {
      //we must put this one sample in a place holder
      
      result.DeCorr_flags[i]=1;
      if(DeCorr_flags[i]==1)
	    {
	      result.DeCorr_sample[i]=DeCorr_sample[i];
	      result.DeCorr_weight[i]=DeCorr_weight[i];
	    }
      else if(rhs.DeCorr_flags[i]==1)
	    {
	      result.DeCorr_sample[i]=rhs.DeCorr_sample[i];
	      result.DeCorr_weight[i]=rhs.DeCorr_weight[i];
	    }
      else
	    {
	      result.DeCorr_sample[i]=carry_up_sample;
	      result.DeCorr_weight[i]=carry_up_weight;
	    }
      carry_up=0;
      carry_up_sample=0.0;
      carry_up_weight=0.0;
    }
    if(binary_sum==2) 
    { 
      //we just made a whole new sample for this level
      //we will also push this sample up to the next higher level
      
      result.DeCorr_flags[i]=0;
      if(DeCorr_flags[i]==0)
	    {
	      // carry+rhs
	      carry_up_sample=(carry_up_sample*carry_up_weight
                         +rhs.DeCorr_sample[i]*rhs.DeCorr_weight[i])
        /(carry_up_weight+rhs.DeCorr_weight[i]);
	      carry_up_weight=(carry_up_weight+rhs.DeCorr_weight[i])/2.0;
	    }
      else if(rhs.DeCorr_flags[i]==0)
	    {
	      // carry+self
	      carry_up_sample=(carry_up_sample*carry_up_weight
                         +DeCorr_sample[i]*DeCorr_weight[i])
        /(carry_up_weight+DeCorr_weight[i]);
	      carry_up_weight=(carry_up_weight+DeCorr_weight[i])/2.0;
	    }
      else
	    {
	      // rhs+self
	      carry_up_sample=(rhs.DeCorr_sample[i]*rhs.DeCorr_weight[i]
                         +DeCorr_sample[i]*DeCorr_weight[i])
        /(rhs.DeCorr_weight[i]+DeCorr_weight[i]);
	      carry_up_weight=(rhs.DeCorr_weight[i]+DeCorr_weight[i])/2.0;
	    }
      carry_up=1;
      result.DeCorr[i].newSample(carry_up_sample,carry_up_weight);
    }
    if(binary_sum==3) 
    {
      //we made a sample for the next level plus have an extra
      //we will put in a place holder until later 
      
      result.DeCorr_flags[i]=1;
      result.DeCorr_sample[i]=carry_up_sample;
      result.DeCorr_weight[i]=carry_up_weight;
      carry_up=1;
      carry_up_sample=(rhs.DeCorr_sample[i]*rhs.DeCorr_weight[i]
                       +DeCorr_sample[i]*DeCorr_weight[i])
        /(rhs.DeCorr_weight[i]+DeCorr_weight[i]);
      carry_up_weight=(rhs.DeCorr_weight[i]+DeCorr_weight[i])/2.0;
      result.DeCorr[i].newSample(carry_up_sample,carry_up_weight);
    }
  }
  
  return  result;
}

void QMCProperty::toXML(ostream& strm)
{
  // Open XML
  strm << "<QMCProperty>" << endl;
  
  // DCL=30.  This means that 2^30=1e9 iterations are necessary for there to be
  // one sample in the last element of the array.  Most of our calculations
  // will run for a few million iterations at most, so most of the data we are
  // printing out in the checkpoint files will be zeros.  I am changing this
  // to only print out the elements of the array that have samples.  This will
  // make the checkpoint files smaller and faster to read and write.

  // Figure out which elements of the array have samples.

  int active_elements = 0;

  for (int i=0; i<DCL; i++)
    {
      if (DeCorr[i].getNumberSamples() == 0)
	break;
      else 
	active_elements++;
    }

  strm << "<NumberOfElements>" << endl;
  strm << active_elements << endl;
  strm << "</NumberOfElements>" << endl;

  // statistics
  strm << "<DeCorrStatistics>" << endl;
  for(int i=0; i<active_elements; i++)
    DeCorr[i].toXML(strm);
  strm << "</DeCorrStatistics>" << endl;
  
  // decorr flags
  strm << "<DeCorrFlags>" << endl;
  for(int i=0; i<active_elements; i++)
    strm << DeCorr_flags[i] << endl;
  strm << "</DeCorrFlags>" << endl;
  
  // decorr sample
  strm << "<DeCorrSamples>" << endl;
  for(int i=0; i<active_elements; i++)
    strm << DeCorr_sample[i] << endl;
  strm << "</DeCorrSamples>" << endl;
  
  // decorr weight
  strm << "<DeCorrWeights>" << endl;
  for(int i=0; i<active_elements; i++)
    strm << DeCorr_weight[i] << endl;
  strm << "</DeCorrWeights>" << endl;
  
  // Close XML
  strm << "</QMCProperty>" << endl;
}

bool QMCProperty::readXML(istream& strm)
{
  string temp;
  
  // Open XML
  strm >> temp;
  if (temp != "<QMCProperty>")
    return false;

  // Get the number of active elements.

  strm >> temp;
  if (temp != "<NumberOfElements>")
    return false;
  strm >> temp;
  int active_elements = atoi(temp.c_str());
  strm >> temp;
  if (temp != "</NumberOfElements>")
    return false;

  // statistics
  strm >> temp;
  if (temp != "<DeCorrStatistics>")
    return false;
  for(int i=0; i<active_elements; i++)
    {
      if (!DeCorr[i].readXML(strm))
        return false;
    }
  strm >> temp;
  if (temp != "</DeCorrStatistics>")
    return false;

  // decorr flags
  strm >> temp;
  if (temp != "<DeCorrFlags>")
    return false;
  for(int i=0; i<active_elements; i++)
    {
      strm >> temp;
      DeCorr_flags[i] = atoi(temp.c_str());
    }
  strm >> temp;
  if (temp != "</DeCorrFlags>")
    return false;

  // decorr sample
  strm >> temp;
  if (temp != "<DeCorrSamples>")
    return false;
  for(int i=0; i<active_elements; i++)
    {
      strm >> temp;
      DeCorr_sample[i] = atof(temp.c_str());
    }
  strm >> temp;
  if (temp != "</DeCorrSamples>")
    return false;

  // decorr weight
  strm >> temp;
  if (temp != "<DeCorrWeights>")
    return false;
  for(int i=0; i<active_elements; i++)
    {
      strm >> temp;
      DeCorr_weight[i] = atof(temp.c_str());
    }
  strm >> temp;
  if (temp != "</DeCorrWeights>")
    return false;
  
  for (int i=active_elements; i<DCL; i++)
    DeCorr[i].zeroOut();

  // Close XML
  strm >> temp;

  if (temp != "</QMCProperty>")
    return false;

  return true;
}
  
void QMCProperty::printAll(ostream & strm)
{
  strm << *this;
  int width = 14;
  strm.precision(6);

  strm << setw(width) << "DeCorr_depth"
       << setw(width) << "samples"
       << setw(width) << "Ave" 
       << setw(width) << "Std"
       << setw(width) << "StdStd"
       << setw(width) << "Var" 
       << setw(width) << "VarStd"
       << setw(width) << "Corr. Len."
       << setw(width) << "Variance"
       << setw(width) << "Skew"
       << setw(width) << "Kurtosis"
       << endl;
  
  for(int i=0;i<DCL;i++)
    {
      if( DeCorr[i].getNumberSamples() > 0 )
	{
	  strm.setf(ios::fixed);
	  strm.unsetf(ios::scientific);
	  strm << setw(width) << i << setw(width) << DeCorr[i].getNumberSamples() 
	       << setw(width) << DeCorr[i].getAverage();
	  
	  if( DeCorr[i].getNumberSamples() > 1 )
	    {
	      strm.unsetf(ios::fixed);
	      strm.setf(ios::scientific);
	      strm << setw(width) << getBlockStandardDeviation(i)
		   << setw(width) << getBlockStandardDeviationStandardDeviation(i)
		   << setw(width) << getBlockVariance(i)
		   << setw(width) << getBlockVarianceStandardDeviation(i);

	      strm.setf(ios::fixed);
	      strm.unsetf(ios::scientific);
	      strm << setw(width) << getBlockVariance(i) / getBlockVariance(0)
		   << setw(width) << DeCorr[i].getVariance()
		   << setw(width) << getBlockSkewness(i)
		   << setw(width) << getBlockKurtosis(i);
	    }
	  
	  strm << endl;
	}  
    }
}

ostream& operator <<(ostream& strm, QMCProperty &rhs)
{    
  if(!false)
  {
    strm.precision(12);
    strm.width(20);
    strm << scientific << rhs.getAverage() << " +/- ";
    if( fabs(rhs.getStandardDeviation()) > 1e-300 )
      if( log( fabs( rhs.getStandardDeviation() )) > 10.0)
	strm << scientific;
    strm.precision(12);
    strm.width(20);
    strm << scientific << rhs.getStandardDeviation() << " (" << rhs.getNumberSamples() << " samples, decorr depth " << rhs.getDecorrDepth() << ")" << endl;
  } else {
    printf("%20.12e +/- %20.12e (%li samples, decorr depth %i)\n",rhs.getAverage(),rhs.getStandardDeviation(),rhs.getNumberSamples(),rhs.getDecorrDepth());
  }

  /*
   double a = rhs.stdevFittingParameters[0] * rhs.stdevFittingParameters[0];
   double b = rhs.stdevFittingParameters[1] * rhs.stdevFittingParameters[1];
   double c = rhs.stdevFittingParameters[2] * rhs.stdevFittingParameters[2];
   double d = rhs.stdevFittingParameters[3] * rhs.stdevFittingParameters[3];
   
   
   strm << "StdDev Fit: stdev(x) = " << a << "*exp(" << c << "*x+" 
   << d << "*x**2)/(1.0+" << b << "*exp(" << c << "*x+" 
   << d << "*x**2))" << endl;
   
   strm << endl;
   */
  return strm;
}

double QMCProperty::getBlockVariance(int i)
{
  double v = DeCorr[i].getVariance();
  double n = DeCorr[i].getNumberSamples();

  return v/(n-1);
}

double QMCProperty::getBlockSkewness(int i)
{
  double v = DeCorr[i].getSkewness();
  //double n = DeCorr[i].getNumberSamples();

  return v;
}

double QMCProperty::getBlockKurtosis(int i)
{
  double v = DeCorr[i].getKurtosis();
  //double n = DeCorr[i].getNumberSamples();

  return v;
}

double QMCProperty::getBlockVarianceStandardDeviation(int i)
{
  double var = getBlockVariance(i);
  
  double n     = DeCorr[i].getNumberSamples();
  return var*sqrt(2.0/(n-1.0));  
}

double QMCProperty::getBlockStandardDeviation(int i)
{
  return sqrt( fabs(getBlockVariance(i)) );
}

double QMCProperty::getBlockStandardDeviationStandardDeviation(int i)
{
  double stdev = getBlockStandardDeviation(i);
  double n     = DeCorr[i].getNumberSamples();
  return stdev/sqrt(2.0*(n-1));
}


void QMCProperty::calculateObjectiveFunction(Array1D<double> & params, 
                                             Array1D<double> & standardDeviations, 
                                             Array1D<double> & standardDeviationsErrors,
                                             double & functionValue, 
                                             Array1D<double> & gradientValue)
{
  /*
   
   gradientValue.allocate(numberParams);
   
   gradientValue = 0.0;
   functionValue = 0.0;
   
   for(int logBlockSize=0; logBlockSize < standardDeviations.dim1(); 
       logBlockSize++)
   {
     double x = (double) logBlockSize;
     
     double weight = 1.0;//1.0/standardDeviationsErrors(logBlockSize);//1.0/(standardDeviationsErrors(logBlockSize) * 
                         // pow(2.0,-x/2) );
       
       double expFactor1 = exp( params(2)*params(2) * x + 
                                params(3)*params(3) * x * x ); 
       
       // use the cutoff in calculating expFactor2 to avoid overflow errors
       
       double expFactor2 = expFactor1 < 1e10 ?
         expFactor1/(1.0 + params(1)*params(1)*expFactor1) : 
         1.0/(params(1)*params(1));
       
       double functionValueLocal = params(0)*params(0) * expFactor2;
       
       double residual = standardDeviations(logBlockSize) - functionValueLocal;
       
       functionValue += residual*residual*weight*weight;
       
       double gradFactor = -2.0*residual*weight*weight;
       
       Array1D<double> gradF(numberParams);
       
       gradF(0) = 2.0 * params(0) * expFactor2;
       gradF(1) = -2.0 * params(0)*params(0) * params(1) * 
         expFactor2*expFactor2;
       gradF(2) = 2.0 * params(2) * x * functionValueLocal * 
         ( 1.0  - params(1)*params(1) * expFactor2 ); 
       gradF(3) = 2.0 * params(3) * x*x * functionValueLocal * 
         ( 1.0  - params(1)*params(1) * expFactor2 ); 
       
       for(int i=0; i<numberParams; i++)
       {
         gradientValue(i) = gradientValue(i) + gradFactor*gradF(i);
       }
   }
   
   // force the parameters to a smaller value
   
   //const double constant = 0.001;
   
   //functionValue += constant * ( pow(params(2),4) + pow(params(3),4) );
   
   //gradientValue(2) += constant * 4 * pow(params(2),3);
   //gradientValue(3) += constant * 4 * pow(params(3),3);
   
   
   const double constant0 = 0.0;
   const double constant1 = 1.0e-5;
   functionValue += constant0 * ( params(0)*params(0) + params(1)*params(1) ) + 
   constant1 * ( params(2)*params(2) + params(3)*params(3) );
   
   gradientValue(0) += constant0 * 2 * params(0);
   gradientValue(1) += constant0 * 2 * params(1);
   gradientValue(2) += constant1 * 2 * params(2);
   gradientValue(3) += constant1 * 2 * params(3);
   
   */
}      


void QMCProperty::calculateObjectiveFunction(Array1D<double> & params, 
                                             Array1D<double> & standardDeviations,
                                             Array1D<double> & standardDeviationsErrors,
                                             double & functionValue)
{
  /*
   
   functionValue = 0.0;
   
   for(int logBlockSize=0; logBlockSize < standardDeviations.dim1(); 
       logBlockSize++)
   {
     double x = (double) logBlockSize;
     
     double weight = 1.0;//1.0/standardDeviationsErrors(logBlockSize); //1.0/(standardDeviationsErrors(logBlockSize) * 
                         //pow(2.0,-x/2) );
       
       double expFactor1 = exp( params(2)*params(2) * x + 
                                params(3)*params(3) * x * x ); 
       
       // use the cutoff in calculating expFactor2 to avoid overflow errors
       
       double expFactor2 = expFactor1 < 1e10 ?
         expFactor1/(1.0 + params(1)*params(1)*expFactor1) : 
         1.0/(params(1)*params(1));
       
       
       double functionValueLocal = params(0)*params(0) * expFactor2;
       
       double residual = standardDeviations(logBlockSize) - functionValueLocal;
       
       functionValue += residual*residual*weight*weight;
   }
   
   
   const double constant0 = 0.0;
   const double constant1 = 1.0e-5;
   functionValue += constant0 * ( params(0)*params(0) + params(1)*params(1) ) + 
   constant1 * ( params(2)*params(2) + params(3)*params(3) );
   
   */
}      


double QMCProperty::calculateLineSearchObjectiveFunction(
                                                         Array1D<double> & params,
                                                         Array1D<double> & searchDirection,
                                                         double stepLength,
                                                         Array1D<double> & standardDeviations,
                                                         Array1D<double> & standardDeviationsErrors)
{
  Array1D<double> p;
  
  p = params + (searchDirection*stepLength);
  
  double functionValue = 1e50;
  
  calculateObjectiveFunction(p, standardDeviations, standardDeviationsErrors,
                             functionValue);
  
  return functionValue;
} 

double QMCProperty::calculateLineSearchObjectiveFunctionDerivative(
                                                                   Array1D<double> & params,
                                                                   Array1D<double> & searchDirection,
                                                                   double stepLength,
                                                                   Array1D<double> & standardDeviations,
                                                                   Array1D<double> & standardDeviationsErrors)
{
  Array1D<double> p;
  
  p = params + (searchDirection*stepLength);
  
  double functionValue;
  Array1D<double> grad;
  
  calculateObjectiveFunction(p, standardDeviations, standardDeviationsErrors,
                             functionValue, grad);
  
  return grad*searchDirection;
}


double QMCProperty::cubicInterpolateStep(double a_lo, double a_hi, 
                                         double phi_0, 
                                         double phi_a_lo, double phi_a_hi, 
                                         double phi_prime_0)
{
  // Here is a bunch of nonsense to find the cubic interpolated a_new
  // See the book for details
  
  // make sure the problem isn't singular!
  if( fabs(a_lo) < 1e-8 ) a_lo = 1e-8;
  if( fabs(a_hi) < 1e-8 ) a_hi = 1e-8;
  if( fabs(a_lo - a_hi) < 1e-8) return (a_lo+a_hi)/2.0; 
  
  Array2D<double> M(2,2);
  
  M(0,0) = a_lo * a_lo;
  M(0,1) = -a_hi * a_hi;
  M(1,0) = -a_lo * a_lo * a_lo;
  M(1,1) = a_hi * a_hi * a_hi;
  
  Array1D<double> v(2);
  
  v(0) = phi_a_hi - phi_0 - a_hi * phi_prime_0;
  v(1) = phi_a_lo - phi_0 - a_lo * phi_prime_0;
  
  double factor = a_lo * a_lo * a_hi * a_hi * ( a_hi - a_lo );
  
  Array1D<double> ab(2);
  
  for(int i=0; i<2; i++)
  {
    ab(i) = 0.0;
    
    for(int j=0; j<2; j++)
    {
      ab(i) = M(i,j) * v(j);
    }
  }
  
  ab /= factor;
  
  double desc  = ab(1)*ab(1)-3.0*ab(0)*phi_prime_0;
  
  if( desc < 0 ) return (a_hi+a_lo)/2.0;
  
  double a_new = (-ab(1)+sqrt(desc))/3.0/ab(0);
  
  // if the interpolated point isn't in the range, return the midpoint
  if( a_new < a_lo ) a_new = (a_lo + a_hi)/2.0;
  if( a_new > a_hi ) a_new = (a_lo + a_hi)/2.0;
  
  // if the new point is too close to either end scoot it over
  double interptol = 1e-6;
  
  if( fabs(a_new-a_lo) < interptol * fabs(a_hi-a_lo) )
  {
    a_new = a_lo + interptol;
  }
  
  if( fabs(a_new-a_hi) < interptol * fabs(a_hi-a_lo) )
  {
    a_new = a_hi - interptol;
  }
  
  return a_new;
}


double QMCProperty::zoom(double a_lo, double a_hi, double phi_0, 
                         double phi_a_lo, double phi_a_hi, double phi_prime_0,
                         Array1D<double> & params, 
                         Array1D<double> & searchDirection,
                         Array1D<double> & standardDeviations,
                         Array1D<double> & standardDeviationsErrors)
{
  const double zoomTol = 1e-7;   // numerical stability tolerance
  const int maxZoomSteps = 1000; // max steps in converging zoom
  const double c1 = 1.0e-4;      // sufficient decrease Wolfe parameter
  const double c2 = 0.1;         // curvature Wolfe parameter
  
  double a;
  
  for(int i=0; i<maxZoomSteps; i++)
  {
    // Check for neumerically unstable problem
    if( fabs(a_lo-a_hi) < zoomTol ) return (a_lo+a_hi)/2.0;
    
    // find a point in [a_lo,a_hi]
    a = cubicInterpolateStep( a_lo, a_hi, phi_0, phi_a_lo, 
                              phi_a_hi, phi_prime_0);
    
    // evaluate phi(a)
    double phi_a = 
      calculateLineSearchObjectiveFunction(params,searchDirection,a,  
                                           standardDeviations,standardDeviationsErrors);
    
    if( phi_a > phi_0 + c1 * a * phi_prime_0 || phi_a >= phi_a_lo )
    {
      a_hi     = a;
      phi_a_hi = phi_a;
    }
    else
    {
      double phi_prime_a = 
	    calculateLineSearchObjectiveFunctionDerivative(params, 
                                                     searchDirection,a,standardDeviations,standardDeviationsErrors);
      
      if( fabs(phi_prime_a) <= -c2 * phi_prime_0 )
	    {
	      return a;
	    }
      else if( phi_prime_a * (a_hi - a_lo) >= 0 )
	    {
	      a_hi     = a_lo;
	      phi_a_hi = phi_a_lo;
	    }
      
      a_lo     = a;
      phi_a_lo = phi_a;
    }
  }
  
  cerr << "WARNING: zoom(...) in QMCProperty did not converge!" << endl;
  
  return a;  
}


double QMCProperty::wolfeStepLength(double alpha_guess, 
                                    Array1D<double> & params,
                                    Array1D<double> & searchDirection,
                                    Array1D<double> & gradient,
                                    double functionValue,
                                    Array1D<double> & standardDeviations,
                                    Array1D<double> & standardDeviationsErrors)
{
  const double alpha_max = 2.0;    // maximum allowed step length
  const int MaxWolfeSteps = 1000;  // max steps to converge step length
  const double c1 = 1.0e-4;        // sufficient decrease Wolfe parameter
  const double c2 = 0.1;           // curvature Wolfe parameter
  
  
  
  // This is all right out of the book so look it up to make sense of it.
  
  double a_max = alpha_max;
  double phi_max = 
    calculateLineSearchObjectiveFunction(params,searchDirection,alpha_max,  
                                         standardDeviations,standardDeviationsErrors);
  
  double a_0 = 0.0;
  double a_1 = alpha_guess;
  
  double phi_0       = functionValue;
  double phi_prime_0 = gradient*searchDirection;
  
  double phi_a_0 = phi_0;
  double phi_a_1;
  double phi_prime_a_1;
  
  for( int i=0; i<MaxWolfeSteps; i++ )
  {
    phi_a_1 = 	
    calculateLineSearchObjectiveFunction(params,searchDirection,a_1,  
                                         standardDeviations,standardDeviationsErrors);
    
    // Check for sufficient descent
    if( phi_a_1 > phi_0 + c1 * a_1 * phi_prime_0 ||
        ( phi_a_1 >= phi_a_0 && i>0 ) )
    {
      return zoom(a_0,a_1,phi_0,phi_a_0,phi_a_1,phi_prime_0,params,
                  searchDirection,standardDeviations,
                  standardDeviationsErrors);
    }
    
    // Calculate the directional derivative at the trial point
    phi_prime_a_1 = 
	    calculateLineSearchObjectiveFunctionDerivative(params, 
                                                     searchDirection,a_1,standardDeviations,
                                                     standardDeviationsErrors);
    
    // Check the curvature condition
    if( fabs(phi_prime_a_1) <= -c2 * phi_prime_0 )
    {
      return a_1;
    }
    
    if( phi_prime_a_1 >= 0 )
    {
      return zoom(a_1,a_0,phi_0,phi_a_1,phi_a_0,phi_prime_0,params,
                  searchDirection,standardDeviations,
                  standardDeviationsErrors);
    }
    
    // move a(i) to a(i-1)
    a_0     = a_1;
    phi_a_0 = phi_a_1;
    
    // choose next a in (a_i, a_max)
    a_1 = cubicInterpolateStep(a_0,a_max,phi_0,phi_a_0,phi_max,
                               phi_prime_0);
  }
  
  cerr << "WARNING: wolfeStepLength(...) in QMCProperty did not converge!" 
    << endl;
  
  return a_1;  
}


void QMCProperty::generateInitialGuessFittingParameters()
{
  /*
   
   // use the data generated so far to guess the initial guess parameters
   // ***** these can be tuned to improve performance ****
   
   if( DeCorr[4].getNumberSamples() > 2 )
   {
     double temp_0 = getBlockStandardDeviation(0);
     double temp_N = getBlockStandardDeviation(4);
     
     temp_N = temp_0 > temp_N ? 2*temp_0 : temp_N;
     
     // deals with the case of a flat slope when estimating the parameters
     
     if( fabs(temp_N-temp_0) < 1.0e-6 )
     {
       stdevFittingParameters[1] = 100.0;
     }
     else
     {
       stdevFittingParameters[1] = sqrt(temp_0/(temp_N-temp_0));
     }
     
     stdevFittingParameters[0] = sqrt(temp_0*(1.0+stdevFittingParameters[1]*
                                              stdevFittingParameters[1]));
     stdevFittingParameters[2] = 1.0e-1;
     stdevFittingParameters[3] = 1.0e-2;
   }
   
   */
}


#ifdef PARALLEL

bool QMCProperty::mpiTypeCreated = false;

void QMCProperty::buildMpiReduce()
{
  MPI_Op_create((MPI_User_function*)Reduce_Function,
                true,&MPI_REDUCE);
}

void QMCProperty::buildMpiType()
{
  QMCProperty indata;
  
  int          block_lengths[4];
  MPI_Aint     displacements[4];
  MPI_Aint     addresses[5];
  MPI_Datatype typelist[4];
  
  typelist[0] = QMCStatistic::MPI_TYPE;     //this is the DeCorr[DCL]
  typelist[1] = MPI_INT;                    //this is the DeCorr_flags[DCL]
  typelist[2] = MPI_DOUBLE;                 //this is the DeCorr_sample[DCL]
  typelist[3] = MPI_DOUBLE;                 //this is the DeCorr_weight[DCL]
  
  block_lengths[0] = DCL;
  block_lengths[1] = DCL;
  block_lengths[2] = DCL;
  block_lengths[3] = DCL;
  
  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.DeCorr[0]), &addresses[1]);
  MPI_Address(&(indata.DeCorr_flags[0]), &addresses[2]);
  MPI_Address(&(indata.DeCorr_sample[0]), &addresses[3]);
  MPI_Address(&(indata.DeCorr_weight[0]), &addresses[4]);
  
  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];
  displacements[3] = addresses[4] - addresses[0];

  //MPI_Type_struct(4, block_lengths, displacements, typelist, &MPI_TYPE);

  //This is the same sort of fix that was applied to QMCProperties
  MPI_Datatype temp;
  MPI_Type_struct(4, block_lengths, displacements, typelist, &temp);
  MPI_Type_create_resized(temp,0,sizeof(QMCProperty),&MPI_TYPE);   

  MPI_Type_commit(&MPI_TYPE);
}


MPI_Datatype QMCProperty::MPI_TYPE;

MPI_Op QMCProperty::MPI_REDUCE;

void QMCProperty::Reduce_Function(QMCProperty *in, QMCProperty *inout, 
                                  int *len, MPI_Datatype *dptr)
{
  for(int i=0; i < *len; i++)
    inout[i] = inout[i] + in[i];
}


#endif

