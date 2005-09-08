/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "GPUQMCMatrix.h"

#ifdef QMC_GPU

/**********DEBUGGING PARAMETERS**********/
//#define PRINT_TIMINGS
static const bool PRINT_SHADER     = false;
static const  int TIMING_REPS      =  10;
static const bool INT_FINISHES     = false;

/**********SHADER PARAMETERS*************/
/**
if unrolled, max allowed for LOOPS_PER_PASS:
(this list needs to be filled in)
SEPARATE_MAD    = 144
TOGETHER_MAD    > 144
KAHAN_SUMS      ~ 55
MRT 2x2 SepM    ~
MRT 1x2 SepM    ~
MRT 2x1 SepM    ~
MRT 2x2 KS      ~
MRT 1x2 KS      ~
MRT 2x1 KS      ~
if rolled must be less than 256 (it won't give an error if too large -- just the wrong answer)
*/
static const  int LOOPS_PER_PASS   = 144;
static const bool UNROLL_LOOP      = !true; //this penalizes the constructor
static const bool USE_CHEAPER      = true;  //changes the precision of the coordinate vars. with this, 4 registers only?
static const bool USE_TRIANGLES    = true;  //Slightly improves performance
static const bool CHEAPER_FIRST    = false; //something is wrong with this.
static const bool ALL_HALF         = false; //faster, but precision is dead. error seems to be 0.1

static const  int SEPARATE_MAD     = 1; //Compiles to 2 mads in shader
static const  int TOGETHER_MAD     = 2; //Compiles to 1 mad, 1 add, 1 mul
static const  int KAHAN_SUMS       = 3; //Keeps summation error down to machine error
static const  int EXPERIMENTAL     = 4; //Something else
static const  int HOW_TO_ADD       = TOGETHER_MAD;

/**
1 <= MRT_W*MRT_H <= 4 (max 4 render textures)
For matrices that don't evenly divide by 3 or 4, fringe cases are difficult,
so they haven't been done. So these two must be either 1 or 2.
According to cjiang, the important case is 2x2 anyway...
*/
static const  int MRT_W            = 2;
static const  int MRT_H            = 2;

static const bool REUSE_SHADERS    = false;
static       bool shadersCreated   = false; //this must be false


/**********OpenGL/Cg PARAMETERS**********/
/**This doesn't seem to make a difference*/
#if UNROLL_LOOP
static const CGprofile matrixProfile = CG_PROFILE_FP30;
#else
static const CGprofile matrixProfile = CG_PROFILE_FP40;
#endif

/**
The GL_TEXTURE_RECTANGLE_NV must be used for non-power-of-2 textures. correspondingly,
one must use samplerRECT instead of sampler2D. makes one wonder if blocking optimizations
are available here... i don't know what the performance difference is.
 
GL_FLOAT_RGBA16_NV vs GL_FLOAT_RGBA32_NV? for some reason, using GL_FLOAT_RGBA16_NV is a bit slower...
*/
#if ALL_HALF
#define TEXTURE_INTERNAL_FORMAT    GL_FLOAT_RGBA32_NV
#define TEXTURE_TARGET             GL_TEXTURE_RECTANGLE_NV
#else
#define TEXTURE_INTERNAL_FORMAT    GL_FLOAT_RGBA32_NV
#define TEXTURE_TARGET             GL_TEXTURE_RECTANGLE_NV
#endif

/**********END PARAMETERS****************/

/*class static variable's definition*/
vector<CGprogram> GPUQMCMatrix::fp;
vector<CGparameter> GPUQMCMatrix::tELxBF;
vector<CGparameter> GPUQMCMatrix::tORxBF;
vector< vector<CGparameter> > GPUQMCMatrix::tELxOR;

GPUQMCMatrix::GPUQMCMatrix()
{
  nOrbitals = 0;
  nBasisF = 0;
  nRows = 0;
  nCols = 0;
  resultsIn = -1;
}

GPUQMCMatrix::GPUQMCMatrix(Array1D< Array2D<qmcfloat> > & Coeffs, int numCalcs)
{
  /*>>>>>>>>>>>>>>>>>>>> SET UP CONSTANTS <<<<<<<<<<<<<<<<<<<<*/
  getFactors(numCalcs,nRows,nCols);
  nDets = Coeffs.dim1();
  nOrbitals = Coeffs(0).dim1();
  nBasisF = Coeffs(0).dim2();
  
  deltaBF = (int)(nBasisF/2.0);
  deltaOE = (int)(nOrbitals/2.0);
  if(nBasisF%2 != 0)   deltaBF += 1;
  if(nOrbitals%2 != 0) deltaOE += 1;
  
  numLoops = deltaBF;
  numPasses = (int)(ceil((double)numLoops/LOOPS_PER_PASS));
  
  deltaW = (int)(deltaOE/MRT_W);
  deltaH = (int)(deltaOE/MRT_H);
  if(deltaOE%MRT_W != 0 || deltaOE%MRT_H != 0)
    cout << "Error: MRT fringe cases not implemented\n";
  if(MRT_H*MRT_W < 1 || MRT_H*MRT_W > 4)
    cout << "Error: 1 <= MRT_W*MRT_H <= 4 (max 4 render textures)\n";
    
  /*>>>>>>>>>>>>>>>>>>>> SET UP DATA STRUCTURES <<<<<<<<<<<<<<<<<<<<*/
  //CPU DATA
  //storage used to download data from GPU
  cpuData = new GLfloat[ nCols*deltaOE * nRows*nMats*deltaOE * 4 ];
  //storage used to upload data to GPU
  pixelData = (GLfloat *) calloc( nCols*deltaOE * nRows*nMats*deltaOE * 4 , sizeof(GLfloat) );
  
  //GPU DATA
  gpuDataFB = new GPUQMCFramebuffer[nDets];
  for(int i=0; i<nDets; i++)
    {
      gpuDataFB[i].initialize(nCols*deltaW, nRows*nMats*deltaH, 2, MRT_H*MRT_W);
    }
  getOpenGLError("Error creating framebuffers");
  
  //GPU COEFFICIENT TEXTURE
  coTexID = new GLuint[nDets];
  glGenTextures(nDets, coTexID);
  
  ArrayGPU coeffA = ArrayGPU(1);
  for(int i=0; i<nDets; i++)
    {
      coeffA(0) = Coeffs(i);
      loadData(coeffA, coTexID[i], 1, 1);
    }
  getOpenGLError("Error loading Coeffs");
  
  /*>>>>>>>>>>>>>>>>>>>> SET UP SHADERS <<<<<<<<<<<<<<<<<<<<*/
  if(fp.empty())
    loadShaders();
    
  //this guy records which framebuffer the result is in
  resultsIn = -1;
}

int GPUQMCMatrix::runCalculation(int numCalcs, GLuint bfInput)
{
  if(numCalcs == 0) cout << "ERROR: numCalcs can not be zero\n";
  resultsIn = -1;
  getFactors(numCalcs,nRows,nCols);
  
#ifdef PRINT_TIMINGS
  Stopwatch sw = Stopwatch();
  double temp;
  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
    {
#endif
    
      for(int iDet=0; iDet<nDets; iDet++)
        {
          gpuDataFB[iDet].cleanAllBuffers();
          gpuDataFB[iDet].drawTo(0);
          
          cgGLEnableProfile(matrixProfile);// <-- this line MUST be within a Begin/End Capture pair
          for(int passIndex=0; passIndex<numPasses; passIndex++)
            {
              cgGLSetTextureParameter(tELxBF[passIndex], bfInput);
              cgGLSetTextureParameter(tORxBF[passIndex], coTexID[iDet]);
              
              cgGLEnableTextureParameter(tELxBF[passIndex]);
              cgGLEnableTextureParameter(tORxBF[passIndex]);
            }
            
          int maxs = gpuDataFB[iDet].getWidth();
          int maxt = gpuDataFB[iDet].getHeight();
          
          glEnable(GL_SCISSOR_TEST);
          glScissor( 0, 0, maxs, maxt);
          
          for(int i=0; i<numPasses; i++)
            {
              gpuDataFB[iDet].drawTo(i%2);
              
              int whichRT;
              for(int mrth=0; mrth<MRT_H; mrth++)
                {
                  for(int mrtw=0; mrtw<MRT_W; mrtw++)
                    {
                      whichRT = mrth*MRT_W + mrtw;
                      cgGLSetTextureParameter(tELxOR[i][whichRT], gpuDataFB[iDet].getTextureID((i+1)%2,whichRT));
                      cgGLEnableTextureParameter(tELxOR[i][whichRT]);
                    }
                }
                
              cgGLBindProgram(fp[i]);
              
              if(USE_TRIANGLES)
                {
                  glBegin(GL_TRIANGLES);
                  glTexCoord2f(0.0f, 0.0f);    glVertex2f(-1.0f, -1.0f);
                  glTexCoord2f(0.0f, maxt*2);  glVertex2f(-1.0f, 3.0f);
                  glTexCoord2f(maxs*2, 0.0f);  glVertex2f(3.0f, -1.0f);
                  glEnd();
                }
              else
                {
                  glBegin(GL_QUADS);
                  glTexCoord2f(0.0,  0.0 );  glVertex2f(-1.0, -1.0);
                  glTexCoord2f(0.0,  maxt);  glVertex2f(-1.0,  1.0);
                  glTexCoord2f(maxs, maxt);  glVertex2f( 1.0,  1.0);
                  glTexCoord2f(maxs, 0.0 );  glVertex2f( 1.0, -1.0);
                  glEnd();
                }
            }
            
          gpuDataFB[iDet].drawTo(0);
          
          cgGLDisableProfile(matrixProfile);
          for(int passIndex=0; passIndex<numPasses; passIndex++)
            {
              cgGLDisableTextureParameter(tELxBF[passIndex]);
              cgGLDisableTextureParameter(tORxBF[passIndex]);
              for(int mrt=0; mrt<MRT_H*MRT_W; mrt++)
                {
                  cgGLDisableTextureParameter(tELxOR[passIndex][mrt]);
                }
            }
        }
        
      if(INT_FINISHES)
        {
          glFinish();
          glFlush();
        }
        
      glDisable(GL_SCISSOR_TEST);
      
#ifdef PRINT_TIMINGS
      
    }
  sw.stop();
  temp = (double)sw.timeMS()/TIMING_REPS;
  printf(" mm_cg: %7.2f", temp );
#endif
  
  resultsIn = (numPasses+1)%2;
  
  GET_GLERROR("Error in QMC matrix multiply");
  return 0;
}

void GPUQMCMatrix::loadData(Array1D< Array2D<qmcfloat> > & data, GLuint & textureID, int numRows, int numCols)
{
  glBindTexture(TEXTURE_TARGET, textureID);
  mapData(data,numRows,numCols);
  glTexImage2D(TEXTURE_TARGET, 0, TEXTURE_INTERNAL_FORMAT,
               numCols*deltaBF, numRows*deltaOE, 0, GL_RGBA, GL_FLOAT,
               pixelData);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

inline int GPUQMCMatrix::operandMapping(int r, int c, int i, int j, int numCols)
{
  /**
  r and c indicate which submatrix we're in
  i and j indicate which value in that submatrix
  4*( (row selection)*(total width) + (column selection) );
  */
  return 4*( (r*deltaOE + i)*numCols*deltaBF + (c*deltaBF + j) );
}

void GPUQMCMatrix::mapData(ArrayGPU & data, int numRows, int numCols)
{
  int h = deltaOE;
  int w = deltaBF;
  
  GLfloat fringe = 0.0;
  Array2D<qmcfloat> * matrix;
  
  //we iterate over all the matricies
  for(int c = 0; c < numCols; c++)
    {
      for(int r = 0; r < numRows; r++)
        {
          matrix = & data( c*numRows + r);
          
          int i, j, index;
          for(i=0; i<h-1; i++)
            {
              for(j=0; j<w-1; j++)
                {
                  index = operandMapping(r,c,i,j,numCols);
                  pixelData[index   ] = matrix->get( i*2   , j*2   );
                  pixelData[index +1] = matrix->get( i*2   , j*2+1 );
                  pixelData[index +2] = matrix->get( i*2+1 , j*2   );
                  pixelData[index +3] = matrix->get( i*2+1 , j*2+1 );
                }
            }
            
          j = w-1;
          if(nBasisF%2 != 0)
            {
              for(int i=0; i<h; i++)
                {
                  index = operandMapping(r,c,i,j,numCols);
                  pixelData[index   ] = matrix->get( i*2   , j*2   );
                  pixelData[index +1] = fringe;
                  if(nOrbitals%2 == 0 || i<h-1)
                    pixelData[index +2] = matrix->get( i*2+1 , j*2   );
                  else
                    pixelData[index +2] = fringe;
                  pixelData[index +3] = fringe;
                }
            }
          else
            {
              for(int i=0; i<h; i++)
                {
                  index = operandMapping(r,c,i,j,numCols);
                  pixelData[index   ] = matrix->get( i*2   , j*2   );
                  pixelData[index +1] = matrix->get( i*2   , j*2+1 );
                  if(nOrbitals%2 == 0 || i<h-1)
                    {
                      pixelData[index +2] = matrix->get( i*2+1 , j*2   );
                      pixelData[index +3] = matrix->get( i*2+1 , j*2+1 );
                    }
                  else
                    {
                      pixelData[index +2] = fringe;
                      pixelData[index +3] = fringe;
                    }
                }
            }
            
          i = h-1;
          if(nOrbitals%2 != 0)
            {
              for(int j=0; j<w; j++)
                {
                  index = operandMapping(r,c,i,j,numCols);
                  pixelData[index   ] = matrix->get( i*2   , j*2   );
                  if(nBasisF%2 == 0 || j<w-1)
                    pixelData[index +1] = matrix->get( i*2   , j*2+1 );
                  else
                    pixelData[index +1] = fringe;
                  pixelData[index +2] = fringe;
                  pixelData[index +3] = fringe;
                }
            }
          else
            {
              for(int j=0; j<w; j++)
                {
                  index = operandMapping(r,c,i,j,numCols);
                  pixelData[index   ] = matrix->get( i*2   , j*2   );
                  pixelData[index +2] = matrix->get( i*2+1 , j*2   );
                  if(nBasisF%2 == 0 || j<w-1)
                    {
                      pixelData[index +1] = matrix->get( i*2   , j*2+1 );
                      pixelData[index +3] = matrix->get( i*2+1 , j*2+1 );
                    }
                  else
                    {
                      pixelData[index +1] = fringe;
                      pixelData[index +3] = fringe;
                    }
                }
            }
        }
    }
  //*
  if(numCols*numRows > 0 && false)
    {
      cout << "loaded by matrix multiply\n";
      PrintRGBAPixelsBoxE(pixelData,numCols*deltaBF,numRows*deltaOE,20,20,-1,-1,true);
    }
  //cout << "loaded from\n";
  //for(int i=0; i<numRows*numCols; i++)
  // cout << data(i);
  //cout << endl;
  //*/
}

void GPUQMCMatrix::getResults(Array2D< Array2D<qmcfloat>* > & results)
{
  if(resultsIn < 0)
    cout << "Error: call runCalculation first!\n";
    
  Array2D<qmcfloat> * fetch;
  int subSize = 4 * nCols*deltaW * nRows*nMats*deltaH;
  
  if(results.dim1() != nDets)
    {
      cout << "results data structure has incorrect dimensions";
      return;
    }
    
  for(int iDet=0; iDet<nDets; iDet++)
    {
    
#ifdef PRINT_TIMINGS
      Stopwatch sw = Stopwatch();
      sw.reset(); sw.start();
      for(int numReps=0; numReps<TIMING_REPS; numReps++)
        {
#endif
        
          for(int mrth=0; mrth<MRT_H; mrth++)
            {
              for(int mrtw=0; mrtw<MRT_W; mrtw++)
                {
                  int whichRT = mrth*MRT_W + mrtw;
                  gpuDataFB[iDet].readFrom(resultsIn, whichRT );
                  glReadPixels(0, 0, nCols*deltaW, nRows*nMats*deltaH, GL_RGBA, GL_FLOAT,
                               (GLfloat *)(cpuData + whichRT*subSize));
                }
            }
            
          for(int whichRT=0; whichRT<4; whichRT++)
            {
              if(whichRT == 0 && !true)
                {
                  cout << "this rt is " << whichRT << endl;
                  if(!true)
                    {
                      PrintRGBAPixelsBoxE( (GLfloat *)(cpuData + whichRT*subSize),
                                           nCols*deltaW, nRows*nMats*deltaH, -1, 46, -1, 46, true);
                    }
                  else
                    {
                      PrintRGBAPixelsBoxE( (GLfloat *)(cpuData + whichRT*subSize),
                                           nCols*deltaW, nRows*nMats*deltaH, 5, 5, -1, -1, true);
                    }
                }
            }
            
#ifdef PRINT_TIMINGS
            
        }
      sw.stop();
      double temp = (double)sw.timeMS()/TIMING_REPS;
      printf(" mm_reading: %7.2f", temp );
      
      sw.reset(); sw.start();
      for(int numReps=0; numReps<TIMING_REPS; numReps++)
        {
#endif
        
          int rstart, cstart, index;
          for(int r=0; r<nRows*nMats; r++)
            {
              rstart = r*deltaH;
              for(int c=0; c<nCols; c++)
                {
                  cstart = c*deltaW;
                  
                  //begin mess
                  fetch = results(iDet,c*nRows*nMats + r);
                  
                  /*This whole mess is to avoid having the if statements
                  inside both the i and j loops. Having made the change, it
                  doesn't seem to have made much difference
                  used to be:
                  index = 4*( (rstart + i)*nCols*deltaOE + (cstart + j) );
                  */
                  int i, j;
                  for(i=0; i<deltaOE-1; i++)
                    {
                      for(j=0; j<deltaOE-1; j++)
                        {
                          index = 4*( (rstart + i%deltaH)*nCols*deltaW + (cstart + j%deltaW) );
                          index += subSize*( (int)(i/deltaH)*MRT_W + (int)(j/deltaW) );
                          (*fetch)(i*2,j*2)       = cpuData[index   ];
                          (*fetch)(i*2,j*2+1)     = cpuData[index +1];
                          (*fetch)(i*2+1,j*2)     = cpuData[index +2];
                          (*fetch)(i*2+1,j*2+1)   = cpuData[index +3];
                        }
                    }
                    
                  //scanning the right edge
                  j = deltaOE-1;
                  if(nOrbitals%2 != 0)
                    {
                      for(i=0; i<deltaOE; i++)
                        {
                          index = 4*( (rstart + i%deltaH)*nCols*deltaW + (cstart + j%deltaW) );
                          index += subSize*( (int)(i/deltaH)*MRT_W + (int)(j/deltaW) );
                          (*fetch)(i*2,j*2)       = cpuData[index   ];
                          if(i*2+1 < nOrbitals)
                            (*fetch)(i*2+1,j*2)   = cpuData[index +2];
                        }
                    }
                  else
                    {
                      for(i=0; i<deltaOE; i++)
                        {
                          index = 4*( (rstart + i%deltaH)*nCols*deltaW + (cstart + j%deltaW) );
                          index += subSize*( (int)(i/deltaH)*MRT_W + (int)(j/deltaW) );
                          (*fetch)(i*2,j*2)       = cpuData[index   ];
                          (*fetch)(i*2,j*2+1)     = cpuData[index +1];
                          if(i*2+1 < nOrbitals)
                            {
                              (*fetch)(i*2+1,j*2)   = cpuData[index +2];
                              (*fetch)(i*2+1,j*2+1) = cpuData[index +3];
                            }
                        }
                    }
                    
                  //scanning the bottom edge
                  i = deltaOE-1;
                  if(nOrbitals%2 != 0)
                    {
                      for(j=0; j<deltaOE; j++)
                        {
                          index = 4*( (rstart + i%deltaH)*nCols*deltaW + (cstart + j%deltaW) );
                          index += subSize*( (int)(i/deltaH)*MRT_W + (int)(j/deltaW) );
                          (*fetch)(i*2,j*2)       = cpuData[index   ];
                          if(j*2+1 < nOrbitals)
                            (*fetch)(i*2,j*2+1)   = cpuData[index +1];
                        }
                    }
                  else
                    {
                      for(j=0; j<deltaOE; j++)
                        {
                          index = 4*( (rstart + i%deltaH)*nCols*deltaW + (cstart + j%deltaW) );
                          index += subSize*( (int)(i/deltaH)*MRT_W + (int)(j/deltaW) );
                          (*fetch)(i*2,j*2)       = cpuData[index   ];
                          (*fetch)(i*2+1,j*2)     = cpuData[index +2];
                          if(j*2+1 < nOrbitals)
                            {
                              (*fetch)(i*2,j*2+1)   = cpuData[index +1];
                              (*fetch)(i*2+1,j*2+1) = cpuData[index +3];
                            }
                        }
                    }
                  //end mess
                  
                  //this line lowers the final calculated the result.
                  //the numbers had been elevated both in QMCWavefunction.cpp
                  //(coefficients) and in the GPUQMCBasisFunction.cpp
                  //fetch->operator *= ( pow(2.0,0.0));
                  
                }//end c loop
            }//end r loop
            
#ifdef PRINT_TIMINGS
            
        }
      sw.stop();
      temp = (double)sw.timeMS()/TIMING_REPS;
      printf(" mm_copying: %7.2f\n", temp );
#endif
      
    }
  /*
  if(nCols*nRows >= 1 && false){
   cout << "\nunloaded by matrix multiply\n";
   PrintRGBAPixelsBoxE(cpuData,nCols*deltaOE,nRows*nMats*deltaOE);
  }
  cout << "unloaded to\n";
  cout << endl << *results(0,0) << endl;
  //*/
}

string GPUQMCMatrix::generateShader(int start, int stop, bool isFirstPass)
{
  /**
   in this shader:
   a refers to the basis function texture
   b refers to the coefficient texture
   c selects which column we're in (when there multiple multiplications happening at once)
   coord.x is where this shader needs to start for its pass
   coord.w is which column within our matrix we're looking at
  
   could come up with a new version for pass #1 -- e.g. accum not passed in
  */
  
  /**
   This first part creates several multple
   rendertexture dependent portions of the shader.
  */
  
  //This is the struct defining the output data type
  string outputAccumulator =
    "struct outputType {                                                 \n";
  //This portion helps define the inputs to the shader
  string inputAccumulator;
  //Each rendertexture output gets loaded with results from the previous pass
  string loadAccumulator;
  //This is where all the a's and b's are defined for the shader
  string abAccumulator;
  //Each rendertexture needs help adjusting coordinates
  string coordAccumulator;
  
  /*These guys all need their double-letters replaced with an actual index*/
  string outputSingle =
    "  float4 oIIJJ : COLORKK;                                           \n";
  string inputSingle =
    "    uniform samplerRECT  accumIIJJ,                                 \n";
  string loadSingle =
    "    output.oIIJJ = texRECT(accumIIJJ, position);                    \n";
  string abSingle =
    "    float4 LETTERII;                                                \n";
  string coordSingle =
    "    temptype2 coordKK = temptype2(LL*DELTAW,LL*DELTAH);             \n";
    
  for(int mrth=0; mrth<MRT_H; mrth++)
    {
      abAccumulator += abSingle;
      findandreplace(abAccumulator,"LETTER",  "a");
      findandreplace(abAccumulator,"II",  mrth+1);
    }
  for(int mrtw=0; mrtw<MRT_W; mrtw++)
    {
      abAccumulator += abSingle;
      findandreplace(abAccumulator,"LETTER",  "b");
      findandreplace(abAccumulator,"II",  mrtw+1);
    }
    
  for(int mrth=0; mrth<MRT_H; mrth++)
    {
      for(int mrtw=0; mrtw<MRT_W; mrtw++)
        {
          loadAccumulator += loadSingle;
          inputAccumulator += inputSingle;
          outputAccumulator += outputSingle;
          coordAccumulator += coordSingle;
          findandreplace(loadAccumulator,"II",  mrth+1);
          findandreplace(inputAccumulator,"II",  mrth+1);
          findandreplace(outputAccumulator,"II",  mrth+1);
          findandreplace(outputAccumulator,"JJ",  mrtw+1);
          findandreplace(inputAccumulator,"JJ",  mrtw+1);
          findandreplace(loadAccumulator,"JJ",  mrtw+1);
          findandreplace(outputAccumulator,"KK",  mrth*MRT_W + mrtw);
          findandreplace(coordAccumulator,"KK",  mrth*MRT_W + mrtw + 1);
          findandreplace(coordAccumulator,"LL",  mrth*MRT_W + mrtw);
        }
    }
  outputAccumulator +=
    "};                                                                \n\n";
    
  /*Now we're finally ready to start constructing the actual shader*/
  string shader;
  shader += outputAccumulator;
  shader +=
    "outputType main(                                                    \n"
    "    uniform samplerRECT  elbf,                                      \n"
    "    uniform samplerRECT  oebf,                                      \n";
  shader += inputAccumulator;
  shader +=
    "    in temptype2 position : TEX0)                                   \n"
    "{                                                                   \n"
    "    int c = position.x/DELTAW;                                      \n"
    "    ROW_SHIFTER1                                                    \n"
    "    temptype4 coord =                                               \n"
    "        {index + c*deltaBF, position.y ROW_SHIFTER2,                \n"
    "         index,             fmod(position.x, DELTAW)};              \n"
    "    outputType output;                                              \n";
  shader += coordAccumulator;
  shader += abAccumulator;
  shader += loadAccumulator;
  
  //If the basisfunction matrix is larger than our rendertextures, then
  //each row needs to be shifted to make up the difference
  if(MRT_H != 1)
    {
      findandreplace(shader,"ROW_SHIFTER1", "int r = position.y/DELTAH;");
      findandreplace(shader,"ROW_SHIFTER2", "+ r*DELTAH");
    }
  else
    {
      findandreplace(shader,"ROW_SHIFTER1", "");
      findandreplace(shader,"ROW_SHIFTER2", "");
    }
    
  //A few extra variables are needed when using the KAHAN_SUMS method
  if(HOW_TO_ADD == KAHAN_SUMS)
    {
      shader +=
        "    float4 T=0, C=0, Y=0;                                           \n";
    }
    
  string mad;
  switch(HOW_TO_ADD)
    {
        //With Cg1.3 at least, this will compile to 2 mads
        case SEPARATE_MAD:
        {
          mad =
            "      output.oIIJJ += aII.xxzz*bJJ.xzxz;                          \n"
            "      output.oIIJJ += aII.yyww*bJJ.ywyw;                          \n";
          break;
        }
        //With Cg1.3 at least, this will compile to 1 mad, 1 add, 1 mul
        case TOGETHER_MAD:
        {
          mad =
            "      output.oIIJJ += aII.xxzz*bJJ.xzxz + aII.yyww*bJJ.ywyw;      \n";
          break;
        }
        //With this scheme, the error does not scale with the number of loops
        case KAHAN_SUMS:
        {
          mad =
            "      Y = aII.xxzz*bJJ.xzxz - C;                                   \n"
            "      T = output.oIIJJ + Y;                                        \n"
            "      C = (T - output.oIIJJ) - Y;                                  \n"
            "      output.oIIJJ = T;                                            \n"
            "      Y = aII.yyww*bJJ.ywyw - C;                                   \n"
            "      T = output.oIIJJ + Y;                                        \n"
            "      C = (T - output.oIIJJ) - Y;                                  \n"
            "      output.oIIJJ = T;                                            \n";
          break;
        }
        //Playpen
        case EXPERIMENTAL:
        {
          mad =
            "      output.oIIJJ += aII.xxzz*bJJ.xzxz;                           \n"
            "      output.oIIJJ += aII.yyww*bJJ.ywyw;                           \n";
          break;
        }
    }
    
  string aLoader =
    "      aII   = texRECT(elbf, coord.xy + float2(0,coordII.y));           \n";
  string bLoader =
    "      bJJ   = texRECT(oebf, coord.zw + float2(0,coordJJ.x));           \n";
  string innerLooper =
    "      coord.xz++;                                                      \n";
    
  /*This section creates the multiplication and addition part of the shader
  the best way to understand what it is doing is to look at the shader produced*/
  bool bNeeded = true;
  if(UNROLL_LOOP)
    {
      for(int i=start; i<stop; i++)
        {
          for(int mrth=0; mrth<MRT_H; mrth++)
            {
              for(int mrtw=0; mrtw<MRT_W; mrtw++)
                {
                  shader += innerLooper;
                  shader += mad;
                  findandreplace(shader,"II",  mrth+1);
                  findandreplace(shader,"JJ",  mrtw+1);
                }
            }
        }
    }
  else
    {
      shader +=
        "    while(coord.z < stopOps) {                                      \n";
      for(int mrth=0; mrth<MRT_H; mrth++)
        {
          shader += aLoader;
          for(int mrtw=0; mrtw<MRT_W; mrtw++)
            {
              if(bNeeded)
                {
                  shader += bLoader;
                }
              shader += mad;
              findandreplace(shader,"II",  mrth+1);
              findandreplace(shader,"JJ",  mrtw+1);
            }
          bNeeded = false;
        }
      shader += innerLooper;
      shader +=
        "    }                                                               \n";
      findandreplace(shader,"startOps", start);
      findandreplace(shader,"stopOps",  stop);
    }
    
  //The end of the shader
  shader +=
    "    return output;                                                  \n"
    "}                                                                   \n";
    
  findandreplace(shader,"index",    start);
  findandreplace(shader,"deltaOE", deltaOE);
  findandreplace(shader,"deltaBF", deltaBF);
  findandreplace(shader,"DELTAW", deltaW);
  findandreplace(shader,"DELTAH", deltaH);
  findandreplace(shader,"MRTH", MRT_H);
  findandreplace(shader,"MRTW", MRT_W);
  
  //This slightly reduces the memory requirements and doesn't seem to
  //affect the precision
  if(USE_CHEAPER)
    {
      findandreplace(shader,"temptype", "half");
    }
  else
    {
      findandreplace(shader,"temptype", "float");
    }
    
  //This improves performance, but totally kills precision...
  if(ALL_HALF)
    {
      findandreplace(shader,"float", "half");
    }
    
  if(PRINT_SHADER)
    {
      cout << "shader is:\n" << shader << endl;
      getchar();
    }
  return shader;
}

/*This function's purpose is to aid in gflops calculations by QMCSlater*/
int GPUQMCMatrix::getNumIterations()
{
#ifdef PRINT_TIMINGS
  return TIMING_REPS;
#else
return 1;
#endif
}

void GPUQMCMatrix::loadShaders()
{
  tELxOR.resize(numPasses);
  
  for(int i=0; i<numPasses; i++)
    {
      char shaderName[256];
      sprintf(shaderName,"shader_matmul.%d.%d.%d_pass%d.cg-asm",LOOPS_PER_PASS,deltaOE,deltaBF,i);
      
      if(!shadersCreated && !REUSE_SHADERS)
        {
        
          string generatedShader = generateShader(i*LOOPS_PER_PASS,
                                                  i<numPasses-1?(i+1)*LOOPS_PER_PASS:numLoops,
                                                  i==0?true:false);
                                                  
          //This dumps the raw Cg code into a file for viewing it later
          char cgName[256];
          sprintf(cgName,"shader_matmul.%d.%d.%d_pass%d.cg",LOOPS_PER_PASS,deltaOE,deltaBF,i);
          writeShader(generatedShader.c_str(),cgName);
          
          //This compiles the shader
          fp.push_back(cgCreateProgram(g_cgContext, CG_SOURCE,
                                       generatedShader.c_str(),
                                       matrixProfile, "main", NULL));
                                       
          //This dumps the compiled Cg code into a file for viewing it later
          if(fp[i])
            writeShader(cgGetProgramString(fp[i],CG_COMPILED_PROGRAM),shaderName);
            
        }
      else
        {
        
          //This can load compiled Cg code from a file
          fp.push_back(cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
                                               shaderName,
                                               matrixProfile, "main", NULL));
                                               
        }
        
      if(!fp[i])
        {
          cout << "ERROR: Matrix multiply shader did not compile.\n";
          cout << "     : LOOPS_PER_PASS is set to " << LOOPS_PER_PASS << endl;
          exit(1);
        }
        
      //allegedly, if accum doesn't exist, cgGetNamedParameter will return 0
      //This section generates the names of the inputs from the shader
      //so that they can be used to create the cg parameters
      string accumIIJJ = "accumIIJJ";
      string temp;
      for(int mrth=0; mrth<MRT_H; mrth++)
        {
          for(int mrtw=0; mrtw<MRT_W; mrtw++)
            {
              temp = accumIIJJ;
              findandreplace(temp,"II", mrth+1);
              findandreplace(temp,"JJ", mrtw+1);
              tELxOR[i].push_back(cgGetNamedParameter(fp[i], temp.c_str()));
            }
        }
      tELxBF.push_back(cgGetNamedParameter(fp[i], "elbf"));
      tORxBF.push_back(cgGetNamedParameter(fp[i], "oebf"));
      
      cgGLLoadProgram(fp[i]);
    }
    
  shadersCreated = true;
}

void GPUQMCMatrix::destroy()
{
  for(int i=0; i<nDets; i++)
    {
      //gpuDataFB[i].destroy();
    }
    
  glDeleteTextures(1, coTexID);
  
  delete [] cpuData;
  delete [] pixelData;
}
#endif
