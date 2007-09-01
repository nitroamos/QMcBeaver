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

#include "QMCJastrow.h"

//this will print out stopwatch counts of
//how much time Jastrow Electron-Electron took
static const bool showTimings  = false;

//This will print out the ln JEE result
//as well as its laplacian
static const bool printJastrow = false;

QMCJastrow::QMCJastrow()
{
}

QMCJastrow::~QMCJastrow()
{
  for(int i=0; i<grad_sum_U.dim1(); i++)
    grad_sum_U(i).deallocate();

  sum_U.deallocate();
  grad_sum_U.deallocate();
  laplacian_sum_U.deallocate();

  for(int i=0; i<p2_xa.dim1(); i++)
    for(int j=0; j<p2_xa.dim2(); j++)
      p2_xa(i,j).deallocate();

  p_a.deallocate();
  p2_xa.deallocate();
  p3_xxa.deallocate();
}

void QMCJastrow::initialize(QMCInput * input)
{
  Input = input;
  JastrowElectronNuclear.initialize(Input);
  JastrowElectronElectron.initialize(Input);
  
  if (Input->flags.use_three_body_jastrow == 1)
    ThreeBodyJastrow.initialize(Input);

  int walkersPerPass = Input->flags.walkers_per_pass;
  sum_U.allocate(walkersPerPass);
  grad_sum_U.allocate(walkersPerPass);
  laplacian_sum_U.allocate(walkersPerPass);

  int numJW = globalInput.JP.getNumberJWParameters();
  p_a.allocate(walkersPerPass,numJW);
  p2_xa.allocate(walkersPerPass,numJW);
  p3_xxa.allocate(walkersPerPass,numJW);

  for(int i=0; i<p2_xa.dim1(); i++)
    for(int j=0; j<p2_xa.dim2(); j++)
      p2_xa(i,j).allocate(Input->WF.getNumberElectrons(),3);

#ifdef QMC_GPU
  GPUQMCJastrowElectronElectron temp(JastrowElectronElectron, Input->flags.getNumGPUWalkers());
  gpuJEE = temp;
#endif
}

double QMCJastrow::getJastrow(int which)
{
  return exp(sum_U(which));
}

double QMCJastrow::get_p_a(int which, int ai)
{
  return getJastrow(which)*p_a(which,ai);
}

double QMCJastrow::getLnJastrow(int which)
{
  return sum_U(which);
}

double QMCJastrow::get_p_a_ln(int which, int ai)
{
  return p_a(which,ai);
}

Array2D<double> * QMCJastrow::getGradientLnJastrow(int which)
{
  return &grad_sum_U(which);
}

Array2D<double> * QMCJastrow::get_p2_xa_ln(int which, int ai)
{
  return &p2_xa(which,ai);
}

double QMCJastrow::getLaplacianLnJastrow(int which)
{
  return laplacian_sum_U(which);
}

double QMCJastrow::get_p3_xxa_ln(int which, int ai)
{
  return p3_xxa(which,ai);
}

void QMCJastrow::evaluate(Array1D<QMCWalkerData *> &walkerData,
			  Array1D<Array2D<double>*> &X, int num, int start)
{
  evaluate(Input->JP,walkerData,X,num,start);
}

void QMCJastrow::evaluate(Array2D<double> & R)
{
  Array1D< Array2D<double>* > temp(1);
  temp(0) = &R;

  Array1D<QMCWalkerData *> wdArray(1);
  wdArray(0) = 0;
  evaluate(Input->JP,wdArray,temp,1,0);
}

#ifdef QMC_GPU
void QMCJastrow::gpuEvaluate(Array1D<Array2D<double>*> &X, int num)
{
  Array2D<double> * grad_JEN;
  Array2D<double> * grad_3body;
  Array2D<double> * grad_JEE;
  gpuJEE.unloadResults();

  for(int walker = 0; walker < num; walker++)
  {
    JastrowElectronNuclear.evaluate(Input->JP,*X(walker));

    sum_U(walker) =
      JastrowElectronNuclear.getLnJastrow() + gpuJEE.getLnJastrow(walker);

    laplacian_sum_U(walker) =
      JastrowElectronNuclear.getLaplacianLnJastrow() +
      gpuJEE.getLaplacianLnJastrow(walker);

    grad_JEN = JastrowElectronNuclear.getGradientLnJastrow();

    grad_JEE = gpuJEE.getGradientLnJastrow(walker);

    grad_sum_U(walker).allocate(X(walker)->dim1(),3);

    for(int i=0; i<grad_JEE->dim1(); i++)
      for(int j=0; j<grad_JEE->dim2(); j++)
        (grad_sum_U(walker))(i,j) = (*grad_JEE)(i,j) + (*grad_JEN)(i,j);

    if (Input->flags.use_three_body_jastrow == 1)
      {
	ThreeBodyJastrow.evaluate(Input->JP,*X(walker));

	sum_U(walker) += ThreeBodyJastrow.getLnJastrow();
	laplacian_sum_U(walker) += ThreeBodyJastrow.getLaplacianLnJastrow();

	grad_3body = ThreeBodyJastrow.getGradientLnJastrow();
	for (int i=0; i<grad_3body->dim1(); i++)
	  for (int j=0; j<grad_3body->dim2(); j++)
	    (grad_sum_U(walker))(i,j) += (*grad_3body)(i,j);
      }

    if(printJastrow)
    {
      printf("%4d: ",walker);
      printf("lnJEE # %s%18.15e ",gpuJEE.getLnJastrow(walker)<0?"":" ",gpuJEE.getLnJastrow(walker) );
      printf("laplnJEE # %s%18.15g\n",gpuJEE.getLaplacianLnJastrow(walker)<0?"":" ",gpuJEE.getLaplacianLnJastrow(walker) );
    }
  }
}
void QMCJastrow::setUpGPU(GLuint aElectrons, GLuint bElectrons, int num)
{
  static double averageJ = 0, timeJ = 0;
  static double numT = -5;
  static int multiplicity = gpuJEE.getNumIterations();

  Stopwatch sw = Stopwatch();
  if(showTimings)
  {
      sw.reset(); sw.start();
  }
  gpuJEE.runCalculation(aElectrons, bElectrons, num);
  if(showTimings)
  {
    glFinish();
    sw.stop();
    timeJ = sw.timeUS();
    if(numT >= 0) averageJ += sw.timeUS();
    if( num > 1) numT++;
    cout << "gpu ee: " << (int)(timeJ/multiplicity+0.5) << " ( " << (int)(averageJ/(numT*multiplicity)+0.5) << ")\n";
  } 
}
#endif

void QMCJastrow::evaluate(QMCJastrowParameters & JP,
			  Array1D<QMCWalkerData *> &walkerData,
			  Array1D<Array2D<double>*> &X, int num, int start)
{
  static double averageJ = 0, timeJ = 0;
  static double numT = -5;
  static int multiplicity = 1;
  
  Stopwatch sw = Stopwatch();
  sw.reset();

  Array2D<double> * grad_JEN;
  Array2D<double> * grad_3body;
  Array2D<double> * grad_JEE;

  for(int walker = start; walker < start+num; walker++)
    {
      JastrowElectronNuclear.evaluate(JP,walkerData(walker),*X(walker));

      if(showTimings){sw.start();}
      JastrowElectronElectron.evaluate(JP,walkerData(walker),*X(walker));
      if(showTimings){sw.stop();}
      
      sum_U(walker) = JastrowElectronNuclear.getLnJastrow() +
	JastrowElectronElectron.getLnJastrow();

      laplacian_sum_U(walker) = JastrowElectronNuclear.getLaplacianLnJastrow()
	+ JastrowElectronElectron.getLaplacianLnJastrow();

      grad_JEN = JastrowElectronNuclear.getGradientLnJastrow();
      grad_JEE = JastrowElectronElectron.getGradientLnJastrow();

      grad_sum_U(walker).allocate(X(walker)->dim1(),3);

      for(int i=0; i<grad_JEE->dim1(); i++)
	for(int j=0; j<grad_JEE->dim2(); j++)
	  (grad_sum_U(walker))(i,j) = (*grad_JEE)(i,j) + (*grad_JEN)(i,j);

      if (Input->flags.use_three_body_jastrow == 1)
	{
	  ThreeBodyJastrow.evaluate(JP,walkerData(walker),*X(walker));
	  sum_U(walker) += ThreeBodyJastrow.getLnJastrow();
	  laplacian_sum_U(walker) += ThreeBodyJastrow.getLaplacianLnJastrow();

	  grad_3body = ThreeBodyJastrow.getGradientLnJastrow();

	  for (int i=0; i<grad_3body->dim1(); i++)
	    for (int j=0; j<grad_3body->dim2(); j++)
	      (grad_sum_U(walker))(i,j) += (*grad_3body)(i,j);
	}

      if (IeeeMath::isNaN( sum_U(walker) ))
	{
	  cerr << "WARNING: lnJastrow = " << sum_U(walker) << endl;
	  cerr << "lnJastrow is being set to 0 to avoid ruining the " << endl;
	  cerr << "calculation." << endl;
	  sum_U(walker) = 0.0;
	  laplacian_sum_U(walker) = 0.0;
	  grad_sum_U(walker) = 0.0;
	}

      if(Input->flags.calculate_Derivatives == 1)
	{
	  // We have not figured out the parameter derivatives with respect to
	  // the three body Jastrow functions yet.

	  int numEE = globalInput.JP.getNumberEEParameters();
	  for(int ai=0; ai<numEE; ai++)
	    {
	      p_a(walker,ai)    = JastrowElectronElectron.get_p_a_ln(ai);
	      p2_xa(walker,ai)  = *JastrowElectronElectron.get_p2_xa_ln(ai);
	      p3_xxa(walker,ai) = JastrowElectronElectron.get_p3_xxa_ln(ai);
	    }
	  int shift = numEE;
	  int numNE = globalInput.JP.getNumberNEParameters();
	  for(int ai=0; ai<numNE; ai++)
	    {
	      p_a(walker,ai+shift)    = JastrowElectronNuclear.get_p_a_ln(ai);
	      p2_xa(walker,ai+shift)  = *JastrowElectronNuclear.get_p2_xa_ln(ai);
	      p3_xxa(walker,ai+shift) = JastrowElectronNuclear.get_p3_xxa_ln(ai);
	    }
	  shift += numNE;
	  int numJW = globalInput.JP.getNumberJWParameters(); 
	  for(int ai=shift; ai<numJW; ai++)
	    {
	      p_a(walker,ai)    = ThreeBodyJastrow.get_p_a_ln(ai-shift);
	      p2_xa(walker,ai)  = *ThreeBodyJastrow.get_p2_xa_ln(ai-shift);
	      p3_xxa(walker,ai) = ThreeBodyJastrow.get_p3_xxa_ln(ai-shift);
	    }
	}

      if(printJastrow)
        {
          printf("%4d: ",walker);
          printf("lnJEE # %s%18.15e ",JastrowElectronElectron.getLnJastrow()<0?"":" ",JastrowElectronElectron.getLnJastrow() );
          printf("laplnJEE # %s%18.15g\n",JastrowElectronElectron.getLaplacianLnJastrow()<0?"":" ",JastrowElectronElectron.getLaplacianLnJastrow() );
        }
    }
  if(showTimings)
    {
      timeJ = sw.timeUS();
      if(numT >= 0) averageJ += sw.timeUS();
      if( num > 1) numT++;
      cout << "cpu ee: " << (int)(timeJ/multiplicity+0.5) << " ( " << (int)(averageJ/(numT*multiplicity)+0.5) << ")\n";
    }
}

void QMCJastrow::operator=(const QMCJastrow & rhs )
{
  Input = rhs.Input;

  sum_U = rhs.sum_U;
  grad_sum_U = rhs.grad_sum_U;
  laplacian_sum_U = rhs.laplacian_sum_U;
  JastrowElectronNuclear = rhs.JastrowElectronNuclear;

  JastrowElectronElectron = rhs.JastrowElectronElectron;

  if (Input->flags.use_three_body_jastrow == 1)
    ThreeBodyJastrow = rhs.ThreeBodyJastrow;

#ifdef QMC_GPU
  gpuJEE = rhs.gpuJEE;
#endif
}
