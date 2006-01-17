/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#ifndef GPU_JASTROWELECTRONELECTRON
#define GPU_JASTROWELECTRONELECTRON

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
#include "QMCJastrowElectronElectron.h"

class GPUQMCJastrowElectronElectron : public QMCJastrowElectronElectron, GPUGlobals
  {
    public:
    
      GPUQMCJastrowElectronElectron();
      
      GPUQMCJastrowElectronElectron(QMCJastrowElectronElectron jee, int max_calcs);
      
      void operator=(const GPUQMCJastrowElectronElectron &rhs);
      
      GLuint runCalculation(GLuint aElectronsTexID, GLuint bElectronsTexID, int num);
      
      void unloadResults();

      /**
      Gets the value of the natural log of the electron-electron Jastrow 
      function for the last 
      evaluated electronic configuration and parameter set.  
      \f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$

      @return natural log of the electron-electron Jastrow function 
      (\f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$)
      */
      double getLnJastrow(int which);

      /**
      Gets the gradient of the natural log of the electron-electron 
      Jastrow function with
      respect to the cartesian electronic coordinates for the last 
      evaluated electronic configuration and parameter set.  
      \f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$

      @return gradient natural log of the electron-electron Jastrow function 
      (\f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$)
      */
      Array2D<double> * getGradientLnJastrow(int which);

      /**
      Gets the laplacian of the natural log of the electron-electron 
      Jastrow function with
      respect to the cartesian electronic coordinates for the last 
      evaluated electronic configuration and parameter set.  
      \f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$

      @return gradient natural log of the electron-electron Jastrow function 
      (\f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$)
      */
      double getLaplacianLnJastrow(int which);

      int getNumIterations();

      ~GPUQMCJastrowElectronElectron();
      
    private:
    
      void setUpInputs();
      
      string generatePolynomialShader(int which);

      string generateTranslationShader(int which);

      string generateReductionShader();

      string generateGradientReductionShader(int which);

      void translateElectronPositions(GLuint aElectronsTexID, GLuint bElectronsTexID);

      void sumAllJastrowValues();

      void sumGradJastrowValues();

      void unloadData(GPUQMCFramebuffer & fb, int w, int h);

      string coeffToCgString(Array1D<double> & input, string name);
      
      void loadShaders();
      
      /**
       The abstraction of this to it's own function makes it easy to do all
       the shifting of coordinates that the basisfunction calculations needed.
       @param maxs the width of the primative
       @param maxt the height of the primative
       @param sShift how far the coordinates need to be shifted along x
       @param tShift how far the coordinates need to be shifted along y
      */
      void drawPrimative(GLfloat maxs, GLfloat maxt, GLfloat sShift, GLfloat tShift);

      void drawTriangles(GLfloat maxs, GLfloat maxt, int nCols, int nRows);

      static const int aa = 0;
      static const int bb = 1;
      static const int ab = 2;

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
      int numE, numA, numB, numLarger;
      
      /**
       These are the dimensions of the GPU electron data
      */
      int elecW, elecH;
      
      Array1D<double> array_sum;
      Array1D< Array2D<double> > array_grad_sum;
      Array1D<double> array_lap_sum;

      GLfloat * cpuData;

      static vector<CGprogram> mapElectronsCG;
      static CGparameter mixedInputCGP;
      static vector<CGparameter> inputCGP;

      static vector<CGprogram> polynomialCG;
      static vector<CGparameter> polyInputCGP;

      static CGprogram sumReductionCG;
      static vector<CGparameter> sumReductionCGP;

      static vector<CGprogram> gradientReductionCG;
      static vector<CGparameter> gradientReductionCGP;

      GPUQMCFramebuffer * r1r2FB;
      GPUQMCFramebuffer * polynomialFB;
      GPUQMCFramebuffer * finalUandLapUFB;
      GPUQMCFramebuffer * finalGradUFB;
  };
#endif
#endif
