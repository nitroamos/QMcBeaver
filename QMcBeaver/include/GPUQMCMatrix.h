#ifndef GPUMATRIX
#define GPUMATRIX
#ifdef QMC_GPU

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>

#include "GPUGlobals.h"
#include "GPUQMCFramebuffer.h"
#include "Array1D.h"
#include "Array2D.h"
#include "Stopwatch.h"

/**
 This class multiplies QMcBeaver matrices on the GPU.
 
 The Coeff and basisFunctions submatricies have the same shape because the
 Coeff matrix is passed in transposed.
*/
class GPUQMCMatrix : public GPUGlobals
  {
    public:
      GPUQMCMatrix();
      
      /**
       @param Coeffs is the matrix of coefficients
       @param numCalcs is the walkers_per_pass parameter.
      */
      GPUQMCMatrix(Array1D< Array2D<qmcfloat> > & Coeffs, int numCalcs);
      
      /**
       This method was created instead of a regular deconstructor because...
       Was I having trouble deconstructing when I wanted?
      */
      void destroy();
      
      /**
       This method is called after the GPUQMCMatrix object is created to actually run
       the matrix multiplication. It can be called multiple times with different
       basisfunction data.
       @param numCalcs is the number of simulataneous calculations we want.
        During initialization, this is 1, afterwards it is walkers_per_pass.
        This function will split numCalcs into common factors. e.g. 6 goes to 2,3
        The common factors created here must be smaller than the common factors
        it would get with the numCalcs in the constructor. E.g. if original was 6 (2,3)
        and runCalculation is given 5 (1,5) something will crash.
       @param bfInput is the texture ID number containting the basisfunction data
        produced by GPUQMCBasisFunction
       @return (not implemented) an error code... actually, an error will probably be
        given somewhere else.
      */
      int runCalculation(int numCalcs, GLuint bfInput);
      
      /**
       After runCalculation has been called, you are free to call this function to
       collect your results.
       @param results this is the data structure the function will use to put numbers
        into. The top dimensions are arranged like: results(iDet,c*nRows*nMats + r)
      */
      void getResults(Array2D< Array2D<qmcfloat>* > & results);
      
      /**
       When timing OpenGL operations, to eliminate system dependent effects on the timings,
       it can be useful to iterate a couple times. This tells how many times it
       iterates -- useful in gflops calculations.
       @return how many times runCalculation and getResults iterated
      */
      int getNumIterations();
      
    private:
      /**
       This function will create a shader for a specific pass.
       @param start where in the dot-product to start
       @param stop where in the dot-product to stop
       @param isFirstPass not implemented.
      */
      string generateShader(int start, int stop, bool isFirstPass);
      
      /*
       This function write/reads shaders from file
      */
      void loadShaders();
      
      /**
       Currently only used to load the coefficient data to the GPU.
       @param data the coefficient data
       @param textureID is the ID to load the data into
       @param numRows for the Coeff matrix, this is always 1
       @param numCols for the Coeff matrix, this is always 1
      */
      void loadData(ArrayGPU & data, GLuint & textureID, int numRows, int numCols);
      
      /**
       This is a helper function for the mapData function.
       r and c indicate which submatrix we're in
       i and j indicate which value in that submatrix
      */
      inline int operandMapping(int r, int c, int i, int j, int numCols);
      
      /**
       This function helps loadData by mapping the data from the 2D
       cpu form into a 1D texture format.
      
       This function will write into pixelData, so to prevent the need to
       be continually reallocating it, it will be made as large as the
       largest necessary when this class is initialized.
      */
      void mapData(ArrayGPU & data, int numRows, int numCols);
      
      /**
       The texture id's of the coefficient data. One is needed for each determinant
       in a multiconfiguration calculation
      */
      GLuint * coTexID;
      
      /**
       Originally, I had been thinking that I may not want the GPU
       to handle all of psi, grad psi, and lap psi.
      */
      const static int nMats = 5;
      
      /**
       nRows*nCols equals the number of walkers to simultaneously
       process. These are defined with the numCalcs passed to the
       constructor, and it is based on these that memory is allocated.
      */
      int nRows, nCols;
      
      /**
       There is no bijection between the dimensions of CPU data and texture dimensions*/
      int nDets, nOrbitals, nBasisF;
      
      /**
       Internal dimensions of gpuData.
      */
      int deltaOE, deltaBF;
      
      /**
       When using the multiple rendertargets compiler choice, these define
       the new internal dimensions of gpuData. If MRT_W = MRT_H = 1, then
       these are both equal to deltaOE.
      */
      int deltaW, deltaH;
      
      /**
       Compiled version of the shader.
       They are indexed according to which pass it is used for.
      */
      static vector<CGprogram> fp;
      
      /**
       The basisfunction and coefficient input texture parameters to the shader.
       They are indexed according to which pass it is used for.
      */
      static vector<CGparameter> tELxBF, tORxBF;
      
      /**
       The rendertextures used for multipass and gpu output.
       First, they are indexed according to which pass it is used for,
       and second according to which rendertexture.
      */
      static vector< vector<CGparameter> > tELxOR;
      
      /**
       The problem dependent length of each dot-product.
      */
      int numLoops;
      
      /**
       Depends upon the parameter LOOPS_PER_PASS to determine
       how many passes the shader should use to calculate the
       results.
      */
      int numPasses;
      
      /**
       Records which framebuffer runCalculation stored the results
       in. If negative, runCalculation has not been run. Otherwise,
       resultsIn = (numPasses+1)%2;
      */
      int resultsIn;
      
      /**
       This is the GPU data structure meant to hold the result of the calculation
      */
      GPUQMCFramebuffer * gpuDataFB;
      
      /**
       This is the CPU data structure meant to hold the result of the calculation.
       Only one of these is needed because (for now) the data for each
       determinant configuration is downloaded independently.
      */
      GLfloat * cpuData;
      
      /**
       This is the CPU data structure meant to hold the input used in creating the texture.
       For now, it is only used when uploading coefficient data.
      */
      GLfloat * pixelData;
  };
#endif
#endif
  