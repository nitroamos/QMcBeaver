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
//as well as it's laplacian
static const bool printJastrow = false;


void QMCJastrow::initialize(QMCInput * input)
{
  Input = input;

  int walkersPerPass = Input->flags.walkers_per_pass;
  sum_U.allocate(walkersPerPass);
  grad_sum_U.allocate(walkersPerPass);
  laplacian_sum_U.allocate(walkersPerPass);

  JastrowElectronNuclear.initialize(Input);
  JastrowElectronElectron.initialize(Input);

#ifdef QMC_GPU
  GPUQMCJastrowElectronElectron temp(JastrowElectronElectron, Input->flags.getNumGPUWalkers());
  gpuJEE = temp;
#endif
}

double QMCJastrow::getJastrow(int which)
{
  return exp(sum_U(which));
}

double QMCJastrow::getLnJastrow(int which)
{
  return sum_U(which);
}

Array2D<double> * QMCJastrow::getGradientLnJastrow(int which)
{
  return &grad_sum_U(which);
}

double QMCJastrow::getLaplacianLnJastrow(int which)
{
  return laplacian_sum_U(which);
}

void QMCJastrow::evaluate(Array1D<Array2D<double>*> &X, int num, int start)
{
  evaluate(Input->JP,X,num,start);
}

#ifdef QMC_GPU
void QMCJastrow::gpuEvaluate(Array1D<Array2D<double>*> &X, int num)
{
  Array2D<double> * grad_JEN;
  Array2D<double> * grad_JEE;
  gpuJEE.unloadResults();

  for(int walker = 0; walker < num; walker++)
  {
    JastrowElectronNuclear.evaluate(Input->JP,*X(walker));

    sum_U(walker) =
      JastrowElectronNuclear.getLnJastrow() +
      gpuJEE.getLnJastrow(walker);

    laplacian_sum_U(walker) =
      JastrowElectronNuclear.getLaplacianLnJastrow() +
      gpuJEE.getLaplacianLnJastrow(walker);

    grad_JEN = JastrowElectronNuclear.getGradientLnJastrow();
    grad_JEE = gpuJEE.getGradientLnJastrow(walker);

    grad_sum_U(walker).allocate(X(walker)->dim1(),3);

    for(int i=0; i<grad_JEE->dim1(); i++)
    {
      for(int j=0; j<grad_JEE->dim2(); j++)
      {
        (grad_sum_U(walker))(i,j) =
          (*grad_JEE)(i,j) + (*grad_JEN)(i,j);
      }
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
    timeJ = sw.timeMS();
    if(numT >= 0) averageJ += sw.timeMS();
    if( num > 1) numT++;
    cout << "gpu ee: " << (int)(timeJ/multiplicity+0.5) << " ( " << (int)(averageJ/(numT*multiplicity)+0.5) << ")\n";
  } 
}
#endif

void QMCJastrow::evaluate(QMCJastrowParameters & JP, Array1D<Array2D<double>*> &X, int num, int start)
{
  static double averageJ = 0, timeJ = 0;
  static double numT = -5;
  static int multiplicity = 1;
  
  Stopwatch sw = Stopwatch();
  sw.reset();

  Array2D<double> * grad_JEN;
  Array2D<double> * grad_JEE;

  for(int walker = start; walker < start+num; walker++)
    {
      JastrowElectronNuclear.evaluate(JP,*X(walker));

      if(showTimings){sw.start();}
      JastrowElectronElectron.evaluate(JP,*X(walker));
      if(showTimings){sw.stop();}

      sum_U(walker) = JastrowElectronNuclear.getLnJastrow() +
                      JastrowElectronElectron.getLnJastrow();

      laplacian_sum_U(walker) = 
        JastrowElectronNuclear.getLaplacianLnJastrow() +
        JastrowElectronElectron.getLaplacianLnJastrow();

      grad_JEN = JastrowElectronNuclear.getGradientLnJastrow();
      grad_JEE = JastrowElectronElectron.getGradientLnJastrow();

      grad_sum_U(walker).allocate(X(walker)->dim1(),3);

      for(int i=0; i<grad_JEE->dim1(); i++)
        {
          for(int j=0; j<grad_JEE->dim2(); j++)
            {
              (grad_sum_U(walker))(i,j) =
                (*grad_JEE)(i,j) + (*grad_JEN)(i,j);
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
    timeJ = sw.timeMS();
    if(numT >= 0) averageJ += sw.timeMS();
    if( num > 1) numT++;
    cout << "cpu ee: " << (int)(timeJ/multiplicity+0.5) << " ( " << (int)(averageJ/(numT*multiplicity)+0.5) << ")\n";
  }
}

void QMCJastrow::operator=(const QMCJastrow & rhs ){
  sum_U = rhs.sum_U;
  grad_sum_U = rhs.grad_sum_U;
  laplacian_sum_U = rhs.laplacian_sum_U;
  JastrowElectronNuclear = rhs.JastrowElectronNuclear;
  JastrowElectronElectron = rhs.JastrowElectronElectron;

#ifdef QMC_GPU
  gpuJEE = rhs.gpuJEE;
#endif

  Input = rhs.Input;
}
