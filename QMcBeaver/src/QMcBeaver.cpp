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

#include <iostream>
#include <string>

#include "QMCManager.h"

using namespace std;

void qmcbeaver(int argc, char **argv);

int main(int argc, char *argv[])
{
#ifdef PARALLEL

  // MPI initilizations
  if( MPI_Init(&argc,&argv) )
    {
      cerr << "ERROR: MPI_Init Error" << endl;
      exit(1);
    }

#endif

  qmcbeaver(argc,argv);

#ifdef PARALLEL

  if( MPI_Finalize() )
    {
      cerr << "ERROR: MPI_Finalize Error" << endl;
      exit(1);
    }

#endif

  // Added to make the code totally ansi compliant
  return 0;
}
 
void qmcbeaver(int argc, char ** argv)
{
  Stopwatch timer;
  timer.start();
  int width = 19;

  if(argc < 2)
    {
      cerr << "ERROR: No input file given" << endl;
      exit(1);
    }

  QMCManager TheMan;
  TheMan.initialize( argc, argv);

  if( TheMan.getInputData()->flags.zero_out_checkpoint_statistics )
    {
      TheMan.zeroOut();
    }

  if( TheMan.getInputData()->flags.my_rank == 0 )
    {
      cout << "***************  TheMan.run();" << endl;
      cout << setw(10) << "Iteration" << setw(width) << "Eavg" << setw(width) << "Estd" << setw(width)
	<< "Eavg-Estd" << setw(width) << "Eavg+Estd" << setw(width) << "Num. Walkers" 
	<< setw(width) << "Trial Energy" << setw(width) << "dt_effective"
	<< setw(width) << "Weights" << setw(width) << "Total Samples" << endl;
    }

  TheMan.run();

  if( TheMan.getInputData()->flags.my_rank == 0 )
    {
      TheMan.writeRestart();

      cout << "***************  TheMan.print_results();" << endl;

      cout << TheMan;
      *TheMan.getResultsOutputStream() << TheMan;
    }

  int optloops = 0;
  while( TheMan.getInputData()->flags.optimize_Psi && 
	 optloops < TheMan.getInputData()->flags.max_optimize_Psi_steps )
    {
      if( TheMan.getInputData()->flags.my_rank == 0 )
	{
	  cout << "***************  TheMan.optimize(), iteration: "  << (optloops+1) << ";" << endl;
	}

      TheMan.optimize();

      if( TheMan.getInputData()->flags.my_rank == 0 )
	{
	  TheMan.writeRestart();
	  
	  cout << "***************  TheMan.run();" << endl;
	  cout << setw(10) << "iteration" << setw(width) << "Eavg" << setw(width) << "Estd" << setw(width)
	    << "Eavg-Estd" << setw(width) << "Eavg+Estd" << setw(width) << "Num. Walkers" 
	    << setw(width) << "Trial Energy" << setw(width) << "Eff. dt"
	    << setw(width) << "Weights" << endl;
	}

      TheMan.zeroOut();
      TheMan.run();

      if( TheMan.getInputData()->flags.my_rank == 0 )
	{
	  TheMan.writeRestart();

	  cout << "***************  TheMan.print_results();" << endl;

	  cout << TheMan;
 	  *TheMan.getResultsOutputStream() << TheMan;
	}

      optloops++;
    }

  if( TheMan.getInputData()->flags.my_rank == 0 )
    {
      cout << "***************  TheMan.finalize();" << endl;
    }
  TheMan.finalize();

  timer.stop();
  if( TheMan.getInputData()->flags.my_rank == 0 )
    {
      cout << "Run Time: " << timer << endl;
      TheMan.writeTimingData(cout);
 
      *TheMan.getResultsOutputStream() << "RunTime: " << timer << endl;
      TheMan.writeTimingData( *TheMan.getResultsOutputStream() );
   }
}












