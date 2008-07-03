/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "GPUQMCJastrowElectronElectron.h"
#ifdef QMC_GPU

/**********DEBUGGING PARAMETERS**********/
//#define PRINT_TIMINGS
static const bool PRINT_SHADER     = false;
static const  int TIMING_REPS      =  10;
static const bool INT_FINISHES     = false;

static const bool USE_TRIANGLES    = true;   //guess this should be true
static const bool REUSE_SHADERS    = false;  //to use shaders compiled from a previous run
static       bool shadersCreated   = false;  //don't change this

/*Parameters used by all textures*/
#define TEXTURE_INTERNAL_FORMAT    GL_FLOAT_RGBA32_NV
#define TEXTURE_TARGET             GL_TEXTURE_RECTANGLE_NV

/*Declarations for some class static variables*/
vector<CGprogram> GPUQMCJastrowElectronElectron::mapElectronsCG;
vector<CGparameter> GPUQMCJastrowElectronElectron::inputCGP;
CGparameter GPUQMCJastrowElectronElectron::mixedInputCGP;

vector<CGprogram> GPUQMCJastrowElectronElectron::polynomialCG;
vector<CGparameter> GPUQMCJastrowElectronElectron::polyInputCGP;

CGprogram GPUQMCJastrowElectronElectron::sumReductionCG;
vector<CGparameter> GPUQMCJastrowElectronElectron::sumReductionCGP;

vector<CGprogram> GPUQMCJastrowElectronElectron::gradientReductionCG;
vector<CGparameter> GPUQMCJastrowElectronElectron::gradientReductionCGP;

GPUQMCJastrowElectronElectron::GPUQMCJastrowElectronElectron()
{
  nCols = 0; nRows = 0;
}

GPUQMCJastrowElectronElectron::GPUQMCJastrowElectronElectron(
  QMCJastrowElectronElectron jee, int max_calcs) : QMCJastrowElectronElectron(jee)
{
  getFactors(max_calcs,nRows,nCols);
  allocatedRows = nRows;
  allocatedCols = nCols;

  numA = Input->WF.getNumberElectrons(true);
  numB = Input->WF.getNumberElectrons(false);
  numE = numA + numB;
  if(numA > numB) numLarger = numA;
  else numLarger = numB;

  r1r2FB = new GPUQMCFramebuffer[3];
  r1r2FB[aa].initialize(  nCols*numA, nRows*numA, 1, 1);
  r1r2FB[bb].initialize(  nCols*numB, nRows*numB, 1, 1);
  r1r2FB[ab].initialize(  nCols*numA, nRows*numB, 1, 1);

  polynomialFB = new GPUQMCFramebuffer[3];
  polynomialFB[aa].initialize(  nCols*numA, nRows*numA, 1, 2);
  polynomialFB[bb].initialize(  nCols*numB, nRows*numB, 1, 2);
  polynomialFB[ab].initialize(  nCols*numA, nRows*numB, 1, 2);

  finalUandLapUFB = new GPUQMCFramebuffer( nCols*1, nRows*numLarger, 1, 1);
  finalGradUFB = new GPUQMCFramebuffer( nCols*1, nRows*numE, 1, 1);

  array_sum.allocate(max_calcs);
  array_grad_sum.allocate(max_calcs);
  array_lap_sum.allocate(max_calcs);
  for(int i=0; i<max_calcs; i++)
  {
    array_grad_sum(i).allocate(numE,3);
  }

  GET_GLERROR("Error setting up framebuffer");
  
  cpuData = (GLfloat *) calloc( nCols*numLarger * nRows*numLarger * 4 , sizeof(GLfloat) );

  if(mapElectronsCG.empty())
    loadShaders();
}

GPUQMCJastrowElectronElectron::~GPUQMCJastrowElectronElectron()
{
  delete [] cpuData;
  
  for(int i=0; i<array_grad_sum.dim1(); i++)
  {
    array_grad_sum(i).deallocate();
  }
  array_sum.deallocate();
  array_grad_sum.deallocate();
  array_lap_sum.deallocate();

  delete finalUandLapUFB;
  delete finalGradUFB;
  delete [] r1r2FB;
  delete [] polynomialFB;
}

GLuint GPUQMCJastrowElectronElectron::runCalculation(GLuint aElectronsTexID, GLuint bElectronsTexID, int num)
{
  getFactors(num,nRows,nCols);

#ifdef PRINT_TIMINGS
  Stopwatch sw = Stopwatch();
  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
#endif

  translateElectronPositions(aElectronsTexID, bElectronsTexID);

#ifdef PRINT_TIMINGS
  sw.stop();
  double temp = (double)sw.timeUS()/TIMING_REPS;
  printf("   jee_trans: %7.2f", temp );

  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++){
#endif

  for(int i=0; i<3; i++){
    int maxs = polynomialFB[i].getWidth();
    int maxt = polynomialFB[i].getHeight();
    polynomialFB[i].cleanAllBuffers();
    polynomialFB[i].drawTo(0);

    cgGLSetTextureParameter(polyInputCGP[i] , r1r2FB[i].getTextureID(0,0));

    cgGLEnableProfile(g_cgProfile);

    cgGLEnableTextureParameter(polyInputCGP[i]);
    cgGLBindProgram(polynomialCG[i]);

    if(i == ab){
      drawPrimative(maxs,maxt,0,0);
    } else {
      drawTriangles(maxs,maxt,nCols,nRows);
    }

    cgGLDisableProfile(g_cgProfile);

    cgGLDisableTextureParameter(polyInputCGP[i]);
  }

  if(INT_FINISHES) glFinish();

    //for(int i=0; i<3; i++)
    //unloadData(polynomialFB[i],polynomialFB[i].getWidth(),polynomialFB[i].getHeight());
    /*
    polynomialFB[aa].readFrom(0,0);
    glReadPixels(0,0,polynomialFB[aa].getWidth(),polynomialFB[aa].getHeight(),GL_RGBA,GL_FLOAT,cpuData);
    cout << "unloaded by jee poly\n";
    PrintRGBAPixelsBoxE(cpuData,polynomialFB[aa].getWidth(),polynomialFB[aa].getHeight(),-1,-1,-1,-1,true);
    //*/

#ifdef PRINT_TIMINGS
  }
  sw.stop();
  temp = (double)sw.timeUS()/TIMING_REPS;
  printf("   jee_poly: %7.2f\n", temp );

  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
#endif

  sumAllJastrowValues();

#ifdef PRINT_TIMINGS
  sw.stop();
  temp = (double)sw.timeUS()/TIMING_REPS;
  printf("   jee_sum:   %7.2f", temp );

  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
#endif
    
    sumGradJastrowValues();

#ifdef PRINT_TIMINGS
    sw.stop();
    temp = (double)sw.timeUS()/TIMING_REPS;
    printf("   jee_grad: %7.2f\n", temp );
#endif

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glFlush();
    getOpenGLError("Error in Jastrow Electron-Electron calculation");
    return 0;
}

void GPUQMCJastrowElectronElectron::unloadResults()
{
  int index = 0;
  finalUandLapUFB->readFrom(0,0);
  glReadPixels(0,0,finalUandLapUFB->getWidth(),finalUandLapUFB->getHeight(),GL_RGB,GL_FLOAT,cpuData);

  array_sum = 0;
  array_lap_sum = 0;
  for(int r=0; r<nRows; r++)
  {
    for(int c=0; c<nCols; c++)
    {
      for(int i=0; i<numLarger; i++)
      {
        index = 3*( (r*numLarger + i)*nCols + c );
        array_sum( c*nRows + r ) += cpuData[index];
        array_lap_sum( c*nRows + r ) += cpuData[index+1];
      }
    }
  }
  
  finalGradUFB->readFrom(0,0);
  glReadPixels(0,0,finalGradUFB->getWidth(),finalGradUFB->getHeight(),GL_RGB,GL_FLOAT,cpuData);
  for(int r=0; r<nRows; r++)
  {
    for(int c=0; c<nCols; c++)
    {
      for(int i=0; i<numE; i++)
      {
        index = 3*( (r*numE + i)*nCols + c );
        (array_grad_sum(c*nRows + r))(i,0) = cpuData[index];
        (array_grad_sum(c*nRows + r))(i,1) = cpuData[index+1];
        (array_grad_sum(c*nRows + r))(i,2) = cpuData[index+2];
      }
    }
  }
}

void GPUQMCJastrowElectronElectron::translateElectronPositions(GLuint aElectronsTexID, GLuint bElectronsTexID)
{
  for(int i=0; i<3; i++){
    int maxs = r1r2FB[i].getWidth();
    int maxt = r1r2FB[i].getHeight();
    r1r2FB[i].cleanBuffer(0);

    if(i == aa){
      cgGLSetTextureParameter(inputCGP[i] , aElectronsTexID);
    } else if(i == bb){
      cgGLSetTextureParameter(inputCGP[i], bElectronsTexID);
    } else {
      cgGLSetTextureParameter(inputCGP[i], aElectronsTexID);
      cgGLSetTextureParameter(mixedInputCGP, bElectronsTexID);
      cgGLEnableTextureParameter(mixedInputCGP);
    }

    cgGLEnableProfile(g_cgProfile);

    cgGLEnableTextureParameter(inputCGP[i]);
    cgGLBindProgram(mapElectronsCG[i]);

    if(i == ab){
      drawPrimative(maxs,maxt,0,0);
    } else {
      drawTriangles(maxs,maxt,nCols,nRows);
    }

    cgGLDisableProfile(g_cgProfile);

    cgGLDisableTextureParameter(inputCGP[i]);
    if(i==ab)
      cgGLDisableTextureParameter(mixedInputCGP);
  }
      
  if(INT_FINISHES) glFinish();
  glFlush();

  //for(int i=0; i<1; i++)
  //unloadData(r1r2FB[i],r1r2FB[i].getWidth(),r1r2FB[i].getHeight());

  getOpenGLError("Error in Jastrow Electron-Electron translation");
}

void GPUQMCJastrowElectronElectron::sumAllJastrowValues()
{
  int maxs = finalUandLapUFB->getWidth();
  int maxt = finalUandLapUFB->getHeight();
  finalUandLapUFB->cleanBuffer(0);

  for(int i=0; i<3; i++){
    cgGLSetTextureParameter(sumReductionCGP[i], polynomialFB[i].getTextureID(0,0));
    cgGLEnableTextureParameter(sumReductionCGP[i]);
  }

  cgGLEnableProfile(g_cgProfile);
  cgGLBindProgram(sumReductionCG);
  drawPrimative(maxs,maxt,0,0);
  cgGLDisableProfile(g_cgProfile);

  for(int i=0; i<3; i++){
    cgGLDisableTextureParameter(sumReductionCGP[i]);
  }

  if(INT_FINISHES) glFinish();
  glFlush();

  //unloadData(*finalUandLapUFB,finalUandLapUFB->getWidth(),finalUandLapUFB->getHeight());

  getOpenGLError("Error in Jastrow Electron-Electron translation");
}

void GPUQMCJastrowElectronElectron::sumGradJastrowValues()
{
  int maxs = finalGradUFB->getWidth();
  int maxt = finalGradUFB->getHeight();
  finalGradUFB->cleanBuffer(0);

  int tShift;
  for(int i=0; i<2; i++){
    cgGLSetTextureParameter(gradientReductionCGP[2*i], polynomialFB[i==0?aa:bb].getTextureID(0,1));
    cgGLSetTextureParameter(gradientReductionCGP[2*i+1], polynomialFB[ab].getTextureID(0,1));
    cgGLEnableTextureParameter(gradientReductionCGP[2*i]);
    cgGLEnableTextureParameter(gradientReductionCGP[2*i+1]);

    cgGLEnableProfile(g_cgProfile);
    cgGLBindProgram(gradientReductionCG[i]);
    glEnable(GL_SCISSOR_TEST);
    int opp = (i+1)%2;
    for(int r=0; r<nRows; r++){
      tShift = r*numE + i*numA;
      glScissor(0,tShift,maxs,i==0?numA:numB);
      //we need to shift the texture coordinates of the
      //primative because the fbo target is larger than
      //any of the textures we need to lookup from
      drawPrimative(maxs,maxt,0,-numA*(r+i)+r*opp*(numA-numB));
    }

    glDisable(GL_SCISSOR_TEST);
    cgGLDisableProfile(g_cgProfile);

    for(int j=0; j<2; j++){
      cgGLDisableTextureParameter(gradientReductionCGP[2*i+j]);
    }
  }

  if(INT_FINISHES) glFinish();
  glFlush();

  //unloadData(*finalGradUFB,maxs,maxt);

  getOpenGLError("Error in Jastrow Electron-Electron translation");
}

void GPUQMCJastrowElectronElectron::unloadData(GPUQMCFramebuffer & fb, int w, int h)
{
  int num = 10;
  fb.readFrom(0,0);
  glReadPixels(0,0,w,h,GL_RGBA,GL_FLOAT,cpuData);
  cout << "unloaded by jee\n";
  PrintRGBAPixelsBoxE(cpuData,w,h,num,num,-1,-1,true);
}

void GPUQMCJastrowElectronElectron::drawPrimative(GLfloat maxs, GLfloat maxt, GLfloat sShift, GLfloat tShift)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if(USE_TRIANGLES)
    {
      glBegin(GL_TRIANGLES);
      glTexCoord2f(    sShift    , tShift        );  glVertex2f(-1.0f, -1.0f);
      glTexCoord2f(    sShift    , maxt*2+tShift );  glVertex2f(-1.0f,  3.0f);
      glTexCoord2f(maxs*2+sShift , tShift        );  glVertex2f( 3.0f, -1.0f);
      glEnd();
    }
  else
    {
      glBegin(GL_QUADS);
      glTexCoord2f(  sShift    , tShift      );  glVertex2f(-1.0, -1.0);
      glTexCoord2f(  sShift    , maxt+tShift );  glVertex2f(-1.0,  1.0);
      glTexCoord2f(maxs+sShift , maxt+tShift );  glVertex2f( 1.0,  1.0);
      glTexCoord2f(maxs+sShift , tShift      );  glVertex2f( 1.0, -1.0);
      glEnd();
    }
}

void GPUQMCJastrowElectronElectron::drawTriangles(GLfloat maxs, GLfloat maxt, int nCols, int nRows)
{
  int tShift, sShift;
  int deltaW = (int)(maxs/nCols);
  int deltaH = (int)(maxt/nRows);


  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, maxs, 0, maxt, 0, 100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glEnable(GL_SCISSOR_TEST);
  for(int r=0; r<nRows; r++){
    for(int c=0; c<nCols; c++){
      tShift = r*deltaH;
      sShift = c*deltaW;
      glScissor( sShift, tShift, deltaW, deltaH);

      glBegin(GL_TRIANGLES);
      glTexCoord2f( sShift, tShift);
      glVertex2f  ( sShift, tShift);
      glTexCoord2f( sShift, tShift + deltaH );
      glVertex2f  ( sShift, tShift + deltaH );
      glTexCoord2f( sShift + deltaW, tShift + deltaH );
      glVertex2f  ( sShift + deltaW, tShift + deltaH );
      glEnd();
    }
  }
  glDisable(GL_SCISSOR_TEST);
}


string GPUQMCJastrowElectronElectron::generateTranslationShader(int which)
{
  string shader;
  shader +=
    "float4 main(in float2 coords : TEX0,                                   \n"
    "            MIXED_PARAM                                                \n"
    "            uniform samplerRECT  input                                 \n"
    "           ) : COLOR                                                   \n"
    "{                                                                      \n"
    "   int2   rc = coords/float2(WIDTH,HEIGHT);                            \n"
    "   float2 pos = fmod(coords,float2(WIDTH,HEIGHT));                     \n"
    "   int2   eposY = pos/TEXW;                                            \n"
    "   float2 eposX = fmod(pos,TEXW);                                      \n"
    "   float3 electron1 = texRECT(input,float2(TEXW*rc.x + eposX.x,TEXH_I*rc.y + eposY.x)).xyz;        \n"
    "   float3 electron2 = texRECT(MIXIN,float2(TEXW*rc.x + eposX.y,TEXH_M*rc.y + eposY.y)).xyz;        \n"
    "   float4 output;                                                      \n"
    "   output.xyz = normalize(electron2 - electron1);                      \n"
    "   output.w   = length(electron2 - electron1);                         \n"
    "   return output;                                                      \n"
    //"   return float4(coords,pos);                                              \n"
    //"   return float4(eposX,eposY);                                              \n" 
    //"   return float4(TEXW*rc.x + eposX.x, TEXH_I*rc.y + eposY.x,           \n" 
    //"                 TEXW*rc.x + eposX.y, TEXH_M*rc.y + eposY.y);          \n"
    "}                                                                      \n";

  int elecW = 4;
  int elecHA = (int)(numA/4.0);
  if(numA%4 != 0) elecHA += 1;  
  int elecHB = (int)(numB/4.0);
  if(numB%4 != 0) elecHB += 1;  
  findandreplace(shader,"TEXW",elecW);

  switch(which)
  {
  case aa:
    {
      findandreplace(shader,"MIXIN","input");
      findandreplace(shader,"MIXED_PARAM", "");
      findandreplace(shader,"WIDTH",numA);
      findandreplace(shader,"HEIGHT",numA);
      findandreplace(shader,"TEXH_I",elecHA);
      findandreplace(shader,"TEXH_M",elecHA);
      break;
    }
  case bb:
    {
      findandreplace(shader,"MIXIN","input");
      findandreplace(shader,"MIXED_PARAM", "");
      findandreplace(shader,"WIDTH",numB);
      findandreplace(shader,"HEIGHT",numB);
      findandreplace(shader,"TEXH_I",elecHB);
      findandreplace(shader,"TEXH_M",elecHB);
      break;
    }
  case ab:
    {
      findandreplace(shader,"MIXIN","inputMix");
      findandreplace(shader,"MIXED_PARAM", "uniform samplerRECT  inputMix,");
      findandreplace(shader,"WIDTH",numA);
      findandreplace(shader,"HEIGHT",numB);
      findandreplace(shader,"TEXH_I",elecHA);
      findandreplace(shader,"TEXH_M",elecHB);
      break;
    }
  }
  return shader;
}

string GPUQMCJastrowElectronElectron::coeffToCgString(Array1D<double> & input, string name)
{
  string lines =
    "   float NAME[NUM] = {";
  findandreplace(lines,"NAME",name);
  findandreplace(lines,"NUM",input.dim1());
  for(int i=0; i<input.dim1(); i++){
    if(i==0){
      lines += "ENTRY";
    } else {
      lines += ", ENTRY";
    }
    findandreplace(lines,"ENTRY",input(i));
  }
  lines += "};\n";  
  return lines;
}

string GPUQMCJastrowElectronElectron::generatePolynomialShader(int which)
{
  Array1D<double> num, den;
  string shader;
  shader +=
    "struct outputType {                                                    \n"
    "  float4 o1 : COLOR0;                                                  \n"
    "  float4 o2 : COLOR1;                                                  \n"
    "};                                                                     \n"
    "outputType main(in float2 coords : TEX0,                               \n"
    "            uniform samplerRECT  input                                 \n"
    "           ) : COLOR                                                   \n"
    "{                                                                      \n"
    "   float4  r = texRECT(input,coords);                                  \n"
    "   outputType output;                                                  \n";

  switch(which)
  {
  case aa:
    {
      num = Input->JP.getElectronUpElectronUpParameters()->getCorrelationFunction()->getNumeratorCoeffs();
      den = Input->JP.getElectronUpElectronUpParameters()->getCorrelationFunction()->getDenominatorCoeffs();
      break;
    }
  case bb:
    {
      num = Input->JP.getElectronDownElectronDownParameters()->getCorrelationFunction()->getNumeratorCoeffs();
      den = Input->JP.getElectronDownElectronDownParameters()->getCorrelationFunction()->getDenominatorCoeffs();
      break;
    }
  case ab:
    {
      num = Input->JP.getElectronUpElectronDownParameters()->getCorrelationFunction()->getNumeratorCoeffs();
      den = Input->JP.getElectronUpElectronDownParameters()->getCorrelationFunction()->getDenominatorCoeffs();
      break;
    }
  }

  shader += coeffToCgString(num,"num");
  shader += coeffToCgString(den,"den");

  /*
  int n = coeffs.dim1()-1;
  if( n < 0 ) return;
  f = coeffs(n);
  df = 0.0;
  d2f = 0.0;
  for(int i=n-1; i >= 0; i--){
      d2f = d2f*x + df;
      df  = df*x  + f;
      f   = f*x   + coeffs(i);
  }
  d2f *= 2;
  */
  shader +=
    "   float4 dfd2f = 0;                                                   \n"
    "   float2 f = float2(num[NN],den[NN]);                                 \n"
    "   for(int i=NN-1; i >= 0; i--){                                       \n"
    "      dfd2f = dfd2f*r.w + float4(f.x, dfd2f.x, f.y, dfd2f.z);          \n"
    "      f = f*r.w + float2(num[i],den[i]);                               \n"
    "   }                                                                   \n"
    "   dfd2f.yw *= 2.0f;                                                   \n";

  //assumption: num and den are the same length. they don't have to be...
  //it's just more programming than it's worth at the moment
  assert(num.dim1() == den.dim1());
  findandreplace(shader,"NN",num.dim1()-1);

  //                x   y    z   w
  //dfd2f = float4( ap, app, bp, bpp )
  shader +=
    "   dfd2f /= f.y;                                                       \n"
    "   output.o1.x = f.x/f.y;                                              \n"
    "   output.o1.y = dfd2f.x - dfd2f.z*output.o1.x;                        \n"
    "   output.o1.z = dfd2f.y - dfd2f.w*output.o1.x - 2*dfd2f.z*output.o1.y;\n"
    "   output.o1.z = 2.0f*(2.0f * output.o1.y / r.w + output.o1.z);        \n"
    //"   output.o2.xyz = 1;                         \n";
    "   output.o2.xyz = r.xyz * output.o1.y;                                \n";

  /*
  // p   = a/b
  // p'  = a'/b - a b'/b^2
  // p'' = a"/b - 2 a' b'/b^2 + 2a (b')^2 /b^3 -a b"/b^2

  double aDIVb = a/b;
  double bpDIVb = bp/b;
  double apDIVb = ap/b;

  FunctionValue = aDIVb;
  dFunctionValue = apDIVb - bpDIVb*aDIVb;
  d2FunctionValue = (app - bpp*aDIVb)/b + 2*bpDIVb*(bpDIVb*aDIVb - apDIVb);

  laplacian_sum_U += 2.0*(2.0/r *
                     firstDeriv +
                     U_Function->getSecondDerivativeValue());
*/

  shader +=
    "   return output;                                                      \n"
    "}                                                                      \n";
  return shader;
}

string GPUQMCJastrowElectronElectron::generateReductionShader()
{
  string shader;
  shader +=
    "float4 main(in float2 coords : TEX0,                                   \n"
    "            uniform samplerRECT  inputAA,                              \n"
    "            uniform samplerRECT  inputAB,                              \n"
    "            uniform samplerRECT  inputBB                               \n"
    "           ) : COLOR                                                   \n"
    "{                                                                      \n"
    "   float2 sum = 0;                                                     \n"
    "   int2   rc  = coords/float2(WIDTH,HEIGHT);                           \n"
    "   int2   pos = fmod(coords,float2(WIDTH,HEIGHT));                     \n";

    if(numA > numB)
      shader +=
    "   if(pos.y < NUMB)                                                    \n";

    shader +=
    "   for(int i=0; i<NUMA; i++){                                          \n"
    "      sum += texRECT(inputAB,float2(i+rc.x*NUMA,coords.y-rc.y*ALARGER)).xz;   \n"
    "   }                                                                   \n";
    if(numA < numB)
      shader +=
    "   if(pos.y < NUMA)                                                    \n";

    shader +=
    "   for(int i=0; i<pos.y; i++){                                         \n"
    "      sum += texRECT(inputAA,float2(i+rc.x*NUMA,coords.y-rc.y*BLARGER)).xz;   \n"
    "   }                                                                   \n";

    if(numA > numB)
      shader +=
    "   if(pos.y < NUMB)                                                    \n";

    shader +=
    "   for(int i=0; i<pos.y; i++){                                         \n"
    "      sum += texRECT(inputBB,float2(i+rc.x*NUMB,coords.y-rc.y*ALARGER)).xz;   \n"
    "   }                                                                   \n"
    "   return float4(sum,0,0);                                             \n"
    "}                                                                      \n";

  findandreplace(shader,"NUMA",numA);
  findandreplace(shader,"NUMB",numB);
  findandreplace(shader,"WIDTH",1);
  findandreplace(shader,"HEIGHT",numLarger);

  /* When numA != numB, the difference in height between
     the input needs to be taken account. Specifically, if
     numA > numB, then each row will be offset by numA-numB.
     This code reverses the offset.
  */
  findandreplace(shader,"ALARGER",numLarger-numB);
  findandreplace(shader,"BLARGER",numLarger-numA);  
  return shader;
}

string GPUQMCJastrowElectronElectron::generateGradientReductionShader(int which)
{
  string shader;
  shader +=
    "float4 main(in float2 coords : TEX0,                                   \n"
    "            uniform samplerRECT  inputParallel,                        \n"
    "            uniform samplerRECT  inputOpposite                         \n"
    "           ) : COLOR                                                   \n"
    "{                                                                      \n"
    "   float3 sum = 0;                                                     \n"
    "   int2 rc  = coords/float2(WIDTH,HEIGHT);                             \n"
    "   int2 pos = fmod(coords,float2(WIDTH,HEIGHT));                       \n"
    "   for(float i=0; i<NUMO; i++){                                        \n"
    "      LINE                                                             \n"
    "   }                                                                   \n"
    "   for(int i=0; i < pos.y; i++){                                       \n"
    "      sum += texRECT(inputParallel,float2(i+rc.x*NUMP,coords.y)).xyz;  \n"
    "   }                                                                   \n"
    "   for(int i=(rc.y)*NUMP+pos.y+1; i <(rc.y+1)*NUMP; i++){              \n"
    "      sum -= texRECT(inputParallel,float2(pos.y+rc.x*NUMP,i)).xyz;     \n"
    "   }                                                                   \n"
    "   return float4(sum,0);                                               \n"
    //"   return texRECT(inputParallel,float2(pos.y+rc.x*NUMP,(rc.y+1)*NUMP));\n"
    "}                                                                      \n";

  if(which == 0){
    findandreplace(shader,"WIDTH",1);
    findandreplace(shader,"HEIGHT",numA);
    findandreplace(shader,"NUMP",numA);
    findandreplace(shader,"NUMO",numB);
    findandreplace(shader,"LINE","sum -= texRECT(inputOpposite,float2(pos.y+rc.x*NUMA,i+rc.y*NUMB)).xyz;");
  } else {
    findandreplace(shader,"WIDTH",1);
    findandreplace(shader,"HEIGHT",numB);
    findandreplace(shader,"NUMP",numB);
    findandreplace(shader,"NUMO",numA);
    findandreplace(shader,"LINE","sum += texRECT(inputOpposite,float2(i+rc.x*NUMA,pos.y+rc.y*NUMB)).xyz;");
  }
  findandreplace(shader,"NUMA",numA);
  findandreplace(shader,"NUMB",numB);

  return shader;
}

int GPUQMCJastrowElectronElectron::getNumIterations()
{
#ifdef PRINT_TIMINGS
  return TIMING_REPS;
#else
return 1;
#endif
}

void GPUQMCJastrowElectronElectron::operator=(const GPUQMCJastrowElectronElectron & rhs)
{
  nRows = rhs.nRows;
  nCols = rhs.nCols;
  allocatedRows = rhs.allocatedRows;
  allocatedCols = rhs.allocatedCols;

  numE = rhs.numE;
  numA = rhs.numA;
  numB = rhs.numB;
  numLarger = rhs.numLarger;

  elecW = rhs.elecW;
  elecH = rhs.elecH;

  //r1r2FB = rhs.r1r2FB;
  //polynomialFB = rhs.polynomialFB;
  r1r2FB = new GPUQMCFramebuffer[3];
  r1r2FB[aa].initialize(  nCols*numA, nRows*numA, 1, 1);
  r1r2FB[bb].initialize(  nCols*numB, nRows*numB, 1, 1);
  r1r2FB[ab].initialize(  nCols*numA, nRows*numB, 1, 1);

  polynomialFB = new GPUQMCFramebuffer[3];
  polynomialFB[aa].initialize(  nCols*numA, nRows*numA, 1, 2);
  polynomialFB[bb].initialize(  nCols*numB, nRows*numB, 1, 2);
  polynomialFB[ab].initialize(  nCols*numA, nRows*numB, 1, 2);

  finalUandLapUFB = new GPUQMCFramebuffer( nCols*1, nRows*numLarger, 1, 1);
  finalGradUFB = new GPUQMCFramebuffer( nCols*1, nRows*numE, 1, 1);
  
  //finalUandLapUFB = rhs.finalUandLapUFB;
  //finalGradUFB = rhs.finalGradUFB;
  
  array_sum = rhs.array_sum;
  array_grad_sum = rhs.array_grad_sum;
  array_lap_sum = rhs.array_lap_sum;

  cpuData = (GLfloat *) calloc( nCols*numLarger * nRows*numLarger * 4 , sizeof(GLfloat) );
}

double GPUQMCJastrowElectronElectron::getLaplacianLnJastrow(int which)
{
  return array_lap_sum(which);
}

Array2D<double> * GPUQMCJastrowElectronElectron::getGradientLnJastrow(int which)
{
  return &array_grad_sum(which);
}

double GPUQMCJastrowElectronElectron::getLnJastrow(int which)
{
  return array_sum(which);
}

void GPUQMCJastrowElectronElectron::loadShaders()
{
  for(int i=0; i<3; i++)
  {
    char shaderName[256];
    sprintf(shaderName,"shader_jeet.%d.type%d.cg-asm",numE,i);
    if(!shadersCreated && !REUSE_SHADERS)
    {
      //create the shader and write it to a file
      string generatedShader = generateTranslationShader(i);

      char cgName[256];
      sprintf(cgName,"shader_jeet.%d.type%d.cg",numE,i);
      writeShader(generatedShader.c_str(),cgName);

      mapElectronsCG.push_back(cgCreateProgram(g_cgContext, CG_SOURCE,
        generatedShader.c_str(),
        g_cgProfile, "main", NULL));

      if(mapElectronsCG[i])
        writeShader(cgGetProgramString(mapElectronsCG[i],
        CG_COMPILED_PROGRAM),shaderName);
    }
    else
    {
      //create the shader from precompiled text in a file
      mapElectronsCG.push_back(cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
        shaderName,
        g_cgProfile, "main", NULL));
    }

    if(!mapElectronsCG[i])
    {
      cerr << "ERROR: Jastrow Electron-Electron translation shader " << i << " did not compile.\n";
      exit(1);
    }

    inputCGP.push_back(cgGetNamedParameter(mapElectronsCG[i], "input"));
    if(i == ab)
      mixedInputCGP = cgGetNamedParameter(mapElectronsCG[ab], "inputMix");

    cgGLLoadProgram(mapElectronsCG[i]);
  }

  for(int i=0; i<3; i++)
  {
    char shaderName[256];
    sprintf(shaderName,"shader_jee.%d.type%d.cg-asm",numE,i);
    if(!shadersCreated && !REUSE_SHADERS)
    {
      //create the shader and write it to a file
      string generatedShader = generatePolynomialShader(i);

      char cgName[256];
      sprintf(cgName,"shader_jee.%d.type%d.cg",numE,i);
      writeShader(generatedShader.c_str(),cgName);

      polynomialCG.push_back(cgCreateProgram(g_cgContext, CG_SOURCE,
        generatedShader.c_str(),
        g_cgProfile, "main", NULL));

      if(polynomialCG[i])
        writeShader(cgGetProgramString(polynomialCG[i],
        CG_COMPILED_PROGRAM),shaderName);
    }
    else
    {
      //create the shader from precompiled text in a file
      polynomialCG.push_back(cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
        shaderName,
        g_cgProfile, "main", NULL));
    }

    if(!polynomialCG[i])
    {
      cerr << "ERROR: Jastrow Electron-Electron translation shader " << i << " did not compile.\n";
      exit(1);
    }

    polyInputCGP.push_back(cgGetNamedParameter(polynomialCG[i], "input"));

    cgGLLoadProgram(polynomialCG[i]);
  }

  //blocked off. i guess this sort of thing should be functionalized...
  {
    char shaderName[256];
    sprintf(shaderName,"shader_jeer.%d.cg-asm",numLarger);
    if(!shadersCreated && !REUSE_SHADERS)
    {
      //create the shader and write it to a file
      string generatedShader = generateReductionShader();

      char cgName[256];
      sprintf(cgName,"shader_jeer.%d.cg",numLarger);
      writeShader(generatedShader.c_str(),cgName);

      sumReductionCG = cgCreateProgram(g_cgContext, CG_SOURCE,
        generatedShader.c_str(),
        g_cgProfile, "main", NULL);

      if(sumReductionCG)
        writeShader(cgGetProgramString(sumReductionCG,
        CG_COMPILED_PROGRAM),shaderName);
    }
    else
    {
      //create the shader from precompiled text in a file
      sumReductionCG = cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
        shaderName, g_cgProfile, "main", NULL);
    }

    if(!sumReductionCG)
    {
      cerr << "ERROR: Jastrow Electron-Electron reduction shader did not compile.\n";
      exit(1);
    }

    sumReductionCGP.push_back(cgGetNamedParameter(sumReductionCG, "inputAA"));
    sumReductionCGP.push_back(cgGetNamedParameter(sumReductionCG, "inputBB"));
    sumReductionCGP.push_back(cgGetNamedParameter(sumReductionCG, "inputAB"));

    cgGLLoadProgram(sumReductionCG);
  }

  for(int i=0; i<2; i++)
  {
    char shaderName[256];
    sprintf(shaderName,"shader_jeeg.%d.type%d.cg-asm",numE,i);
    if(!shadersCreated && !REUSE_SHADERS)
    {
      //create the shader and write it to a file
      string generatedShader = generateGradientReductionShader(i);

      char cgName[256];
      sprintf(cgName,"shader_jeeg.%d.type%d.cg",numE,i);
      writeShader(generatedShader.c_str(),cgName);

      gradientReductionCG.push_back(cgCreateProgram(g_cgContext, CG_SOURCE,
        generatedShader.c_str(),
        g_cgProfile, "main", NULL));

      if(gradientReductionCG[i])
        writeShader(cgGetProgramString(gradientReductionCG[i],
        CG_COMPILED_PROGRAM),shaderName);
    }
    else
    {
      //create the shader from precompiled text in a file
      gradientReductionCG.push_back(cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
        shaderName,
        g_cgProfile, "main", NULL));
    }

    if(!gradientReductionCG[i])
    {
      cerr << "ERROR: Jastrow Electron-Electron gradient reduction shader " << i << " did not compile.\n";
      exit(1);
    }

    gradientReductionCGP.push_back(cgGetNamedParameter(gradientReductionCG[i], "inputParallel"));
    gradientReductionCGP.push_back(cgGetNamedParameter(gradientReductionCG[i], "inputOpposite"));

    cgGLLoadProgram(gradientReductionCG[i]);
  }

  //we mark the shader as compiled so that (within the same QMcBeaver run) we don't
  //have to recompile/load it.
  shadersCreated = true;
}
#endif
