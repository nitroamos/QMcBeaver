/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#ifndef GPU_BASISFUNCTIONS
#define GPU_BASISFUNCTIONS

#ifdef QMC_GPU

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <string>
#include "GPUGlobals.h"
#include "GPUQMCFramebuffer.h"
#include "Array1D.h"
#include "Array2D.h"
#include "Stopwatch.h"
#include "QMCBasisFunction.h"

/**
 The texture id for that texture which holds all the electron positions.
 Is it ok to have this as a static variable?
*/
static GLuint electronsTexID;

/**
 The texture id for that texture which holds all the other random parameters
 needed for a basis function calculation. (e.g. the a's, b's, klm's, Rc)
*/
static GLuint bfParametersTexID;

/**
 This class is designed to calculate basis functions for several walkers simultaneously.
*/
class GPUQMCBasisFunction : public QMCBasisFunction, GPUGlobals
  {
    public:
    
      GPUQMCBasisFunction();
      
      /**
       This class lets QMCBasisFunction do some of the initial set up work, and then
       GPUQMCBasisFunction sets up the rest.
       @param bf is the QMCBasisFunction to get the details from
       @param numElectrons is the number of electrons for the specific slater determinant
        it was set up for.
       @param max_calcs is the number of walkers_per_pass
      */
      GPUQMCBasisFunction(QMCBasisFunction bf, int numElectrons, int max_calcs);
      
      /**
       Copy constructor.
       Because the GPU data (like the input coefficients) are common to the calculation,
       each GPUQMCBasisFunction does not need it's own copy of the texture ids.
      */
      void operator=(GPUQMCBasisFunction &rhs);
      
      /**
       After being set up, the class can accept sets of electronic positions.
       @param X the set of all the electronic positions, both alpha and beta.
       @param num is how many walkers to calculate for
       @param start where in the array the electrons for this spin start
       @param stop where in the array the electrons for this spin stop
      */
      GLuint runCalculation(Array1D<Array2D<double>*> &X, int num, int start, int stop);
      
      /**
       Because timing shaders seems to be system dependent... This repeats the calculation
       several times to get averaged numbers.
       @return how many times the calculation was repeated
      */
      int getNumIterations();
      
      /**
       The Destroyer.
      */
      ~GPUQMCBasisFunction();
      
    private:
    
      /**
       This function will set up the input textures as needed. After this function is called,
       no more references to QMCBasisFunction elements should be required...
      
       This function will set up the input textures as desired. After this function is called,
       no more references to QMCBasisFunction elements should be required...
       we need to load into bfParametersTexID the following parameters
       n is the length of the contraction
      
       2*n  --- contraction coefficients, double  --- (a, b)
       3    --- momentum indices, int             --- (k, l, m)
       3    --- center of the nucleus, double     --- (Rc)
       = 6 + 2n (where n is typically 1, 3, or 6)
       = 8, 12, 18
      
       for each run, electronsTexID needs to have
       3    --- electron position, double         --- (r)
      
       one of the dimensions of the cpuData structure is shared
      */
      void setUpInputs();
      
      /**
       This function will generate a different shader for each type of calculation
       (psi, grad psi, lap psi).
       @param the id for which calc we want
       @return the shader as a string
      */
      string generateShader(int which);
      
      /**
       A separate shader is used to translate the GPU output into matrix multiplication
       input.
       @param is12 indicates whether the number of electrons mod 4 is 1 or 2.
       @return the shader as a string
      */
      string generateTranslationShader(bool is12);
      
      /**
       This will compile the shaders either from generateShader or from file.
       It will write out the shaders for private viewing later.
      */
      void loadShaders();
      
      /**
       This function oversees the GPU translation of the basis function output
       into matrix multiplication input.
      */
      void translate();
      
      /**
       The abstraction of this to it's own function makes it easy to do all
       the shifting of coordinates that the basisfunction calculations needed.
       @param maxs the width of the primative
       @param maxt the height of the primative
       @param sShift how far the coordinates need to be shifted along x
       @param tShift how far the coordinates need to be shifted along y
      */
      void drawPrimative(GLfloat maxs, GLfloat maxt, GLfloat sShift, GLfloat tShift);
      
      /**
       This will draw all of the electron positions into texture data.
       @param X the electronic positions for all electrons
       @param start where in the array to find our electrons
       @param stop where in the array to find our electrons
      */
      void loadElectronPositions(Array1D<Array2D<double>*> &X, int start, int stop);
      
      /**
       Used by setUpInputs to position the data in the texture.
       @param i selects row
       @param j selects column
       @param h is the height of the data (currently unused)
       @param w is the width of the data
      */
      int mapping(int i, int j, int h, int w);
      
      /**
       This function only has debugging purposes.
       It will read the frame from a rendertexture
       and spit out the results to the terminal.
       @param which the framebuffer we want to investigate
       @param w how much of the framebuffer we want to look at
       @param h how much of the framebuffer we want to look at
      */
      void unloadData(GPUQMCFramebuffer & which, int w, int h);
      
      /**
       Originally, I had been thinking that perhaps I would want
       to run fewer than all 5 on the GPU...
      */
      const static int nMats = 5;
      
      /**
       Some IDs to make some of the code readable.
      */
      const static int psi = 0;
      const static int grx = 1;
      const static int gry = 2;
      const static int grz = 3;
      const static int lap = 4;
      
      /**
       Number of simultaneous walkers currently being processed
       (changed everytime runCalculation is called)
      */
      int nRows, nCols;
      
      /**
       Number of simultaneous walkers the calculation was
       originally allocated for.
      */
      int allocatedRows, allocatedCols;
      
      /**
       There is no bijection between the dimensions of CPU data and texture dimensions
       So these need to be stored.
      */
      int nElectrons, nBasisF;
      
      /**
       These are the sizes of the GPU bf data when
       it's in four by one configuration.
      */
      int fxo_deltaOE, fxo_deltaBF;
      
      /**
       These are the sizes of the GPU bf data when
       it's in two by two configuration.
      */
      int txt_deltaOE, txt_deltaBF;
      
      /**
       These are the dimensions of the GPU electron data
      */
      int elecW, elecH;
      
      /**
       The number of gaussians required for this molecule.
      */
      int maxGaussians;
      
      /**
       This is how large the coefficient texture needs to be
       to accomodate all the a's and b's. It depends upon:
       basisfunctionParamsH = 2 + (int)(maxGaussians/2.0 + 0.5);
      */
      int basisfunctionParamsH;
      
      /**there is one for each type of calculation (psi, grad, lap)*/
      static vector<CGprogram> fragProg;
      static vector<CGparameter> electronsCGP, paramsCGP;
      
      /**for the translated data to give to the matrix multiply*/
      static CGprogram fxo_to_txt_CG;
      static CGparameter fxo_to_txt_CGP;
      
      GPUQMCFramebuffer basisFunctionsFB;
      GPUQMCFramebuffer outputFB;
      
      /**
       Some data storage space used by several functions.
      */
      GLfloat * cpuData;
  };
#endif
#endif
  
  