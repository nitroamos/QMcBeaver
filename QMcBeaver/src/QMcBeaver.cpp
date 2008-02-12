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
#include <signal.h>

#include "QMCManager.h"

static const bool showExtraHeaders = !false;

QMCInput globalInput;

#ifdef QMC_GPU
#include "GPUGlobals.h"

CGcontext g_cgContext;
CGprofile g_cgProfile;

void idle();
void openGLBootStrap();
void reshape(int w, int h);
void cgErrorCallback();

#endif

int argc;
char ** argv;
void qmcbeaver();
void atExitCallback();
void atSignalCallback(int sig);

//If someone else is trapping signals (e.g. by default LAM traps USR2)
//change these to a free signal number. This is not portable! But
//this is just a feature anyway.
//enum signalChannels { CHANNEL1 = SIGUSR1, CHANNEL2 = SIGUSR2 };
//enum signalChannels { CHANNEL1 = SIGUSR1, CHANNEL2 = 40 };
//this is for use on QSC. LSF on QSC does not force the program to die using SIGURG
enum signalChannels { CHANNEL1 = SIGURG, CHANNEL2 = 40 };

using namespace std;

void printCompileInfo(ostream & strm)
{
  strm << "QMcBeaver version " << VERSION << " was compiled in ";

#ifdef SINGLEPRECISION
  strm << "single";
#else
  strm << "double";
#endif
  strm << " precision." << endl;

#ifdef __DATE__
  strm << "Compiled on " << __DATE__ << " " << __TIME__ << endl;
#endif
#ifdef __VERSION__
  strm << "Compiler version " << __VERSION__ << endl;
#endif

  strm << "With the libraries:";
#ifdef PARALLEL
  strm << " MPI";
#endif
#ifdef USEBLAS
  strm << " BLAS";
#endif
#ifdef USELAPACK
  strm << " LAPACK"; 
#endif
#ifdef USESPRNG
  strm << " SPRNG";
#endif
#ifdef USEHDF5
  strm << " HDF5";
#endif
#ifdef QMC_GPU
  strm << endl;
  GPUGlobals::printVersions(strm);
#endif
  strm << endl;
  strm << "Signaling: CHANNEL1 = " << CHANNEL1 << " CHANNEL2 = " << CHANNEL2 << endl;
}

int main(int argcTemp, char *argvTemp[])
{
  argc = argcTemp;
  argv = argvTemp;
  atexit(atExitCallback);

#ifndef QMC_DEBUG
  signal(CHANNEL1,atSignalCallback);
  signal(CHANNEL2,atSignalCallback);
  signal(SIGTERM,atSignalCallback);
  signal(SIGHUP,SIG_IGN);
#endif
 
#ifdef PARALLEL

  // MPI initilizations
  if( MPI_Init(&argc,&argv) )
    {
      cerr << "ERROR: MPI_Init Error" << endl;
      exit(1);
    }

  int rank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if( rank != 0 )
    {
      //The idea is that we need a stream to which warnings
      //can be printed, but only one processor actually needs
      //do to the warning.

      ofstream dump("/dev/null");
      if(dump.good())
	clog.rdbuf(dump.rdbuf());
    } else {
      //point them to the same stream
      //by default, they're different; which means unix pipe operator >
      //won't catch everything
      clog.rdbuf(cout.rdbuf());
    }
#else
  clog.rdbuf(cout.rdbuf());
#endif

  if(argcTemp == 1 || showExtraHeaders)
    {
      printCompileInfo(clog);
      clog << endl << endl;
    }

#ifdef QMC_GPU
  openGLBootStrap();
#else
  qmcbeaver();
#endif

#ifdef PARALLEL

  if( MPI_Finalize() )
    {
      cerr << "ERROR: MPI_Finalize Error" << endl;
      exit(1);
    }

  //coyote will crash for some reason unless we explicitly exit:
  exit(0);
#endif


  // Added to make the code totally ansi compliant
  return 0;
}

void qmcbeaver()
{
  Stopwatch timer;
  timer.start();

  if(argc < 2)
    {
      cerr << "ERROR: No input file given" << endl;
      exit(1);
    }
    
  globalInput.read( string(argv[ 1 ]) );

  QMCManager TheMan;
  TheMan.initialize( argc, argv);

  if( globalInput.flags.zero_out_checkpoint_statistics )
    {
      TheMan.zeroOut();
    }

  bool ok = TheMan.run(false);
  TheMan.writeRestart();
  
  if( globalInput.flags.my_rank == 0 )
    {
      clog << "***************  Print Results (ok = " << ok << ")" << endl;
      clog << TheMan;
      *TheMan.getResultsOutputStream() << TheMan;
    }

  int optloops = 1;
  while( globalInput.flags.optimize_Psi &&
         optloops <= globalInput.flags.max_optimize_Psi_steps &&
	 ok)
    {
      TheMan.writeTimingData(clog);
      clog << "***************  Optimize iteration: "  << optloops << ";" << endl;

      if(optloops > 1) TheMan.resetTimers();

      TheMan.optimize();
      TheMan.zeroOut();
      ok = TheMan.run(globalInput.flags.equilibrate_every_opt_step);

      stringstream save_opt;
      save_opt << globalInput.flags.checkout_file_name << "/"
	       << globalInput.flags.base_file_name << ".opt" << optloops << ".ckmf";

      if(globalInput.flags.a_diag < 0)
	{
	  /*
	    This method requires 2 optloops to complete one optimization cycle.

	    When optloops % 2 == 1, we have a new wavefunction.
	    When optloops % 2 == 0, we will have a VMC result to attach to that WF.

	    So instead of writing it out every odd step, lets write the WF when
	    we have a measurement of how good it is. If we really want the new WF,
	    we still have the other restart file.
	  */
	  if(optloops % 2 == 0)
	    {
	      TheMan.writeRestart(save_opt.str());
	      TheMan.writeRestart();
	    }
	}
      else
	{
	  TheMan.writeRestart(save_opt.str());
	  TheMan.writeRestart();
	}

      if( globalInput.flags.my_rank == 0 )
        {
	  clog << "***************  Print Results (ok = " << ok << ")" << endl;

          clog << TheMan;
          *TheMan.getResultsOutputStream() << TheMan;
        }

      optloops++;
    }

  if(!ok)
    {
      clog << "Error: the calculation was terminated prematurely.\n";
      /*
      clog << "Some possible problems:\n";
      clog << "   1) You didn't request enough time from your queuing system.\n";
      clog << "   2) Your wavefunction has a problem.\n";
      */
    }

  clog << "***************  Finalize" << endl;
  TheMan.finalize();

  timer.stop();

  if( globalInput.flags.my_rank == 0 &&
      globalInput.flags.calculate_bf_density == 1)
    {
      TheMan.writeBFDensity();
    }

  if( globalInput.flags.my_rank == 0 &&
      globalInput.flags.nuclear_derivatives != "none")
    {
      TheMan.writeForces();
    }
    
  if( globalInput.flags.my_rank == 0 )
    {
      TheMan.writeTimingData(clog);
      clog << "Wallclock Time:         " << timer << endl;
      
      TheMan.writeTimingData( *TheMan.getResultsOutputStream() );
      *TheMan.getResultsOutputStream()
	<< "Wallclock Time:         " << timer << endl;
    }
}

void atExitCallback()
{
#if defined(_WIN32) && !defined(__CYGWIN__)
  clog << "Press any key...\n";
  getchar();
#endif
}

void atSignalCallback(int signal)
{
  signalType choice;

  switch(signal){
  case CHANNEL1:
    choice = SIG_REDUCE;
    break;
  case CHANNEL2:
    choice = SIG_INCREASE;
    break;
  case SIGTERM:
    //cout << "SIGTERM";
    choice = SIG_QUIT;
    break;
  default:
    cout << "No response available for trapped signal " << signal << "\n.";
    return;
    break;
  }

  QMCManager::receiveSignal(choice);
}

#ifdef QMC_GPU
void openGLBootStrap()
{
  //first we need to make sure that we have support for required extensions
  GPUGlobals::checkExtensions();

  ////////Initialize GLUT and OpenGL render context
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
  glutInitWindowSize(200, 1);
  glutInitWindowPosition(500, 500);
  glutCreateWindow("Hello, QMC-GPU!");
  glutHideWindow();

  glDrawBuffer(GL_BACK);
  glReadBuffer(GL_BACK);
  glDisable (GL_DEPTH_TEST);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glDrawBuffer(GL_FRONT);
  glReadBuffer(GL_FRONT);
  glDisable(GL_DEPTH_TEST);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  int err = glewInit();
  if (GLEW_OK != err)
    {
      // problem: glewInit failed, something is seriously wrong
      fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
      exit(-1);
    }

  /*The QMcBeaver program is set up as the OpenGL idle callback function
  As soon as OpenGL is finished setting up (finished processing events in
  the queue) then it will idle and call the qmcbeaver function. At this
  point, because there are no keyboard or mouse callbacks, there shouldn't
  be any interuptions...*/
  glutIdleFunc(idle);
  glutReshapeFunc(reshape);
  GPUGlobals::getOpenGLError(0);
  ////////End initialize GLUT and OpenGL

  ////////Initialize Cg
  cgSetErrorCallback(cgErrorCallback);
  g_cgContext = cgCreateContext();

  // get the best profile for this hardware
  g_cgProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  assert(g_cgProfile != CG_PROFILE_UNKNOWN);
  cgGLSetOptimalOptions(g_cgProfile);
  ////////End Initialize Cg

  /*Once this is called, the whole program is placed on the GLUT "poll-for-events"
  stack. It will not return from this function call maybe unless there is an error*/
  glutMainLoop();
}

// Called when Cg detects an error
void cgErrorCallback()
{
  CGerror lastError = cgGetError();

  if(lastError)
    {
      printf("%s\n\n", cgGetErrorString(lastError));
      printf("%s\n", cgGetLastListing(g_cgContext));
      printf("Cg error!\n");
    }
}

//This function is wired to the glutIdleFunc hook. It will be called when glut finishes setting up.
void idle()
{
  //The entire application is run in one idle function call.
  //idle will only be called when OpenGL is finished setting up.
  qmcbeaver();
#ifdef PARALLEL
  if( MPI_Finalize() )
    {
      cerr << "ERROR: MPI_Finalize Error" << endl;
      exit(1);
    }
#endif
  exit(0);
}

//must call this function if the size of the screen has been changed
void reshape(int w, int h)
{
  glViewport(0, 0, (GLsizei)(w), (GLsizei)(h));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
#endif

