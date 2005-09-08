/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#ifndef GPUGLOBALS
#define GPUGLOBALS

#ifdef QMC_GPU

#include <windows.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <Cg/cgGL.h>
#include "Array2D.h"
#include <sstream>
#include <fstream>

typedef Array1D< Array2D<qmcfloat> > ArrayGPU;

//this is a modified version of something originally one of the OpenGL examples
#define GET_GLERROR(message)                                        \
  {                                                                 \
    string temp = __FILE__;                                         \
    temp = temp.substr(temp.find_last_of('\\')+1);                  \
    GLenum err = glGetError();                                      \
    if (err != GL_NO_ERROR) {                                       \
        fprintf(stderr, "\n[%s line %d] GL Error: %s\n",            \
                temp.c_str(), __LINE__, gluErrorString(err));       \
        fprintf(stderr, "       %s\n\n", message);                  \
        fflush(stderr);                                             \
      }                                                             \
  }

extern CGcontext g_cgContext;
extern CGprofile g_cgProfile;

static int gpu_make;

const static int nvidia_gpu   = 1;
const static int ati_gpu      = 2;
const static int labs3d_gpu   = 3;
const static int unknown_gpu  = 4;

class GPUGlobals
  {
    public:
      /**
       If there is an error, msg will be printed
       I still haven't figured out which I like better. macro or function...
       I'm leaning towards macro cause it can incorporate __LINE__ and because
       I can turn of all macros with a compiler-time switch, so this
       function might be deleted
      */
      static void getOpenGLError(char *msg);
      
      /**
       This will print pixels in a box form.
       r g
       b a
       uses the 'f' format specifier in printf.
       @param pix an array of pixels
       @param w internal width of pix
       @param h internal height of pix
      */
      static void PrintRGBAPixelsBoxF(float* pix, int w, int h);
      
      /**
       This will print pixels in a box form.
       r g
       b a
       uses the 'e' format specifier in printf.
       @param pix an array of pixels
       @param w internal width of pix
       @param h internal height of pix
       @param wBorder print only the first wBorder columns and the last wBorder columns
       @param hBorder print only the first hBorder rows and the last wBorder rows
       @param maxJ print no more than maxJ rows
       @param maxI print no more than maxI columns
      */
      static void PrintRGBAPixelsBoxE(float* pix, int w, int h,
                                      int wBorder, int hBorder,
                                      int maxJ, int maxI,
                                      bool printRC);
                                      
      /**
       This function will print pixels to a terminal screen in a column format.
       this representation isn't as wide in a terminal screen.
       @param pix an array of pixels
       @param w internal width of pix
       @param h internal height of pix
      */
      static void PrintRGBAPixelsColumn(float* pix, int w, int h);
      
      /**
       This function will print pixels to a terminal screen in a column format.
       this representation isn't as wide in a terminal screen.
       this prints only 3 components!
       @param pix an array of pixels
       @param w internal width of pix
       @param h internal height of pix
      */
      static void PrintRGBPixelsColumn(float* pix, int w, int h);
      
      /**
       A quick helper function for Array2D to help it print stuff in the same
       format as my Print*Pixels* functions.
      */
      static void PrintMatrix(Array2D<GLfloat> matrix);
      
      /**
       ld couldn't find this function if it was static and defined in the .cpp file...
       This will scan through a string and replace all instances of one string with another
       templated thing.
       @param source the string to be modified
       @param find what to look for
       @param replace what to replace all instances with. can be an integer, string, etc...
      */
      template<class T> void findandreplace( string& source, const string& find, T replace )
      {
        size_t j;
        stringstream strm;
        strm << replace;
        for (;(j = source.find( find )) != string::npos;)
          {
            source.replace( j, find.length(), strm.str() );
          }
      }
      
      /**
       This will find the two factors closest to each other,
       the most 'square' way to divide the number of walkers into rows and columns.
       for hmx, cols (fact2) can be as high as 7, rows (fact1) can be as high as 21
       @param num how many walkers to divide
       @param fact1 is output, factor number 1
       @param fact2 is output, factor number 1
      */
      static void getFactors(int num, int & fact1, int & fact2);
      
      /**
       Taken from the gpubench source code available from sourceforge.net. Basically helps identify the
       version of the driver and which company made the GPU.
      */
      static const char *gpubench_getDriverVersion(void);
      
      /**
       Prints all the GPU specific info to cout.
      */
      static void printVersions();
      
      /**
       Just dumps a char * into a file.
       This will clobber that file if it already exists.
       @param shader the shader
       @param filename the name of the file to create.
      */
      static void writeShader(const char * shader, const char * filename);
      
      /**
       Makes sure that we're up to date with respect to needed functions.
       Actually, this might not compile if we weren't ok.
      */
      static void checkExtensions();
      
      /**
       This method will print out the limitations placed on shader
       programs by the profile currently loaded
      */
      static void specTest();
  };
  
  
#endif
#endif
  
  