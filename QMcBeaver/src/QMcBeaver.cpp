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

static const bool showExtraHeaders = false;

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

using namespace std;

int main(int argcTemp, char *argvTemp[])
{
  argc = argcTemp;
  argv = argvTemp;
  atexit(atExitCallback);
#ifdef PARALLEL

  // MPI initilizations
  if( MPI_Init(&argc,&argv) )
    {
      cerr << "ERROR: MPI_Init Error" << endl;
      exit(1);
    }

#endif

#ifdef QMC_GPU
  if(showExtraHeaders) cout << "GPU mode\n";
  openGLBootStrap();
#else

  string precision;
#ifdef SINGLEPRECISION
  precision = "single";
#else
  precision = "double";
#endif

#ifdef USEATLAS
  if(showExtraHeaders) cout << "CPU (" << precision << ") mode, using ATLAS\n";
#else
  if(showExtraHeaders) cout << "CPU (" << precision << ") mode, not using ATLAS\n";
#endif

  qmcbeaver();
#endif

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

void qmcbeaver()
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

  if( TheMan.getInputData()->flags.my_rank == 0 &&
      TheMan.getInputData()->flags.calculate_bf_density == 1)
    {
      TheMan.writeBFDensity();
    }

  if( TheMan.getInputData()->flags.my_rank == 0 )
    {
      cout << "Run Time: " << timer << endl;
      TheMan.writeTimingData(cout);

      *TheMan.getResultsOutputStream() << "RunTime: " << timer << endl;
      TheMan.writeTimingData( *TheMan.getResultsOutputStream() );
    }
}

void atExitCallback()
{
#if defined(_WIN32) && !defined(__CYGWIN__)
  cout << "Press any key...\n";
  getchar();
#endif
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
      getchar();
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

  //some feedback for what has been set up
  GPUGlobals::printVersions();

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

