#include "GPUQMCBasisFunction.h"
#ifdef QMC_GPU

/**********DEBUGGING PARAMETERS**********/
//#define PRINT_TIMINGS
static const bool PRINT_SHADER     = false;
static const  int TIMING_REPS      =  10;
static const bool INT_FINISHES     = false;

/**********SHADER PARAMETERS*************/
/** NOTE:
 Raising EXP_SHIFT will increase the accuracy because it helps
 the calculation avoid underflows (which the GPU is particularly
 sensitive to). However, it if is too high, then overflows can
 become a problem; an error about rejecting a NaN kinetic energy
 will appear in the output... A balance between accuracy and
 fraction of rejects is essential.
 
 It's possibly that this trick will only work for light atoms.
*/
static const  int EXP_SHIFT        = 70;     //if base 2, less than 100, if base e, less than 70
static const  int MULTIPLIER       = 0;      //a change to this requires an additional change in GPUQMCMatrix
static const bool USE_BASE_2       = true;   //it's possible that 2^x is faster than e^x
static const bool NAN_CHECKS       = false;  //these might be necessary for some systems like He, and not necessary for others like HMX
static const bool USE_TRIANGLES    = true;   //guess this should be true
static const bool REUSE_SHADERS    = false;  //to use shaders compiled from a previous run
static       bool shadersCreated   = false;  //don't change this

/*Parameters used by all textures*/
#define TEXTURE_INTERNAL_FORMAT    GL_FLOAT_RGBA32_NV
#define TEXTURE_TARGET             GL_TEXTURE_RECTANGLE_NV

/*Declarations for some class static variables*/
vector<CGprogram> GPUQMCBasisFunction::fragProg;
vector<CGparameter> GPUQMCBasisFunction::electronsCGP;
vector<CGparameter> GPUQMCBasisFunction::paramsCGP;
CGprogram GPUQMCBasisFunction::fxo_to_txt_CG;
CGparameter GPUQMCBasisFunction::fxo_to_txt_CGP;

GPUQMCBasisFunction::GPUQMCBasisFunction(QMCBasisFunction bf, int numElectrons, int max_calcs) : QMCBasisFunction(bf)
{
  getFactors(max_calcs,nRows,nCols);
  allocatedRows = nRows;
  allocatedCols = nCols;
  
  nElectrons = numElectrons;
  nBasisF = N_BasisFunctions;
  
  fxo_deltaBF = nBasisF;
  fxo_deltaOE = (int)(nElectrons/4.0);
  if(nElectrons%4 != 0) fxo_deltaOE += 1;
  
  elecW = 4;
  elecH = (int)(nElectrons/4.0);
  if(nElectrons%4 != 0) elecH += 1;
  
  txt_deltaBF = (int)(nBasisF/2.0);
  txt_deltaOE = (int)(nElectrons/2.0);
  if(nBasisF%2 != 0)   txt_deltaBF += 1;
  if(nElectrons%2 != 0) txt_deltaOE += 1;
  
  //we want to know how many gaussians we're going to need to allocate for
  maxGaussians = 0;
  for (int atom=0; atom<flags->Natoms; atom++)
    {
      for (int j=0; j<BFCoeffs(atom).getNumberBasisFunctions(); j++)
        {
          if(BFCoeffs(atom).N_Gauss(j) > maxGaussians) maxGaussians = BFCoeffs(atom).N_Gauss(j);
        }
    }
  basisfunctionParamsH = 2 + (int)(maxGaussians/2.0 + 0.5);
  
  basisFunctionsFB.initialize(nCols*fxo_deltaBF,nRows*nMats*fxo_deltaOE,1,1);
  outputFB.initialize(nCols*txt_deltaBF,nRows*nMats*txt_deltaOE,1,1);
  basisFunctionsFB.checkFramebufferStatus();
  outputFB.checkFramebufferStatus();
  GET_GLERROR("Error setting up framebuffer");
  
  //gpu input
  glGenTextures(1, &electronsTexID);
  glGenTextures(1, &bfParametersTexID);
  //temp storage for loading the gpu input
  cpuData = (GLfloat *) calloc( nCols*fxo_deltaBF * nRows*nMats*max(fxo_deltaOE,basisfunctionParamsH) * 4 , sizeof(GLfloat) );
  
  setUpInputs();
  
  if(fragProg.empty())
    loadShaders();
}

void GPUQMCBasisFunction::loadShaders()
{
  //all psi, grx, gry, grz, and lap are compiled
  for(int i=0; i<nMats; i++)
    {
      char shaderName[256];
      sprintf(shaderName,"shader_bf.%d.%d_type%d.cg-asm",nElectrons,nBasisF,i);
      if(!shadersCreated && !REUSE_SHADERS)
        {
          //create the shader and write it to a file
          string generatedShader = generateShader(i);
          
          char cgName[256];
          sprintf(cgName,"shader_bf.%d.%d_type%d.cg",nElectrons,nBasisF,i);
          writeShader(generatedShader.c_str(),cgName);
          
          fragProg.push_back(cgCreateProgram(g_cgContext, CG_SOURCE,
                                             generatedShader.c_str(),
                                             g_cgProfile, "main", NULL));
                                             
          if(fragProg[i])
            writeShader(cgGetProgramString(fragProg[i],
                                           CG_COMPILED_PROGRAM),shaderName);
                                           
        }
      else
        {
          //create the shader from precompiled text in a file
          fragProg.push_back(cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
                             shaderName,
                             g_cgProfile, "main", NULL));
        }
        
      if(!fragProg[i])
        {
          cout << "ERROR: Basisfunction shader did not compile.\n";
          exit(1);
        }
        
      electronsCGP.push_back(cgGetNamedParameter(fragProg[i], "epos"));
      paramsCGP.push_back(cgGetNamedParameter(fragProg[i], "params"));
      
      cgGLLoadProgram(fragProg[i]);
    }
    
  //this performs the same process with the translation shader
  char shaderName[256];
  sprintf(shaderName,"shader_trans.%d.%d.cg-asm",nElectrons,nBasisF);
  if(!shadersCreated && !REUSE_SHADERS)
    {
      string generatedShader = generateTranslationShader(nElectrons%4==1 || nElectrons%4 ==2);
      
      char cgName[256];
      sprintf(cgName,"shader_trans.%d.%d.cg",nElectrons,nBasisF);
      writeShader(generatedShader.c_str(),cgName);
      
      fxo_to_txt_CG = cgCreateProgram(g_cgContext, CG_SOURCE,
                                      generatedShader.c_str(),
                                      g_cgProfile, "main", NULL);
                                      
      writeShader(cgGetProgramString(fxo_to_txt_CG,CG_COMPILED_PROGRAM),shaderName);
      
    }
  else
    {
      fxo_to_txt_CG = cgCreateProgramFromFile(g_cgContext, CG_OBJECT,
                                              shaderName,
                                              g_cgProfile, "main", NULL);
    }
    
  if(!fxo_to_txt_CG)
    {
      cout << "ERROR: Translation shader did not compile.\n";
      exit(1);
    }
    
  fxo_to_txt_CGP = cgGetNamedParameter(fxo_to_txt_CG, "input");
  cgGLLoadProgram(fxo_to_txt_CG);
  
  //we mark the shader as compiled so that (within the same QMcBeaver run) we don't
  //have to recompile/load it.
  shadersCreated = true;
}

GPUQMCBasisFunction::~GPUQMCBasisFunction()
{
  //basisFunctionsFB.destroy();
  
  delete [] cpuData;
  //glDeleteTextures(1, electronsTexID);
  //glDeleteTextures(1, bfParametersTexID);
}

GLuint GPUQMCBasisFunction::runCalculation(Array1D<Array2D<double>*> &X, int num, int start, int stop)
{
  getFactors(num,nRows,nCols);
  if(nRows > allocatedRows || nCols > allocatedCols)
    cout << "Error: remainder walkers chose bad dimensions.\n";
    
#ifdef PRINT_TIMINGS
  Stopwatch sw = Stopwatch();
  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
#endif
  
    loadElectronPositions(X,start,stop);
    
#ifdef PRINT_TIMINGS
  sw.stop();
  double temp = (double)sw.timeMS()/TIMING_REPS;
  printf(" bf_loading: %7.2f", temp );
  
  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
    {
#endif
    
      int tShift;
      
      int maxs = basisFunctionsFB.getWidth();
      int maxt = basisFunctionsFB.getHeight();
      basisFunctionsFB.cleanBuffer(0);
      
      cgGLEnableProfile(g_cgProfile);// <-- this line MUST be within a Begin/End Capture pair
      
      glEnable(GL_SCISSOR_TEST);
      for(int program=0; program<nMats; program++)
        {
          cgGLSetTextureParameter(electronsCGP[program], electronsTexID);
          cgGLSetTextureParameter(paramsCGP[program], bfParametersTexID);
          cgGLEnableTextureParameter(electronsCGP[program]);
          cgGLEnableTextureParameter(paramsCGP[program]);
          
          cgGLBindProgram(fragProg[program]);
          
          for(int row=0; row<nRows; row++)
            {
              tShift = (program + row*nMats)*fxo_deltaOE;
              glScissor( 0, tShift, maxs, fxo_deltaOE);
              tShift -= row*fxo_deltaOE;
              drawPrimative(maxs,maxt,0,-tShift);
            }
            
        }
      glDisable(GL_SCISSOR_TEST);
      
      cgGLDisableProfile(g_cgProfile);
      for(int whichType=0; whichType<nMats; whichType++)
        {
          cgGLDisableTextureParameter(electronsCGP[whichType]);
          cgGLDisableTextureParameter(paramsCGP[whichType]);
        }
        
      if(INT_FINISHES)
        {
          glFinish();
        }
        
      GET_GLERROR("Error in QMC basis function calculation");
      
#ifdef PRINT_TIMINGS
      
    }
  sw.stop();
  temp = (double)sw.timeMS()/TIMING_REPS;
  printf(" bf_cg: %7.2f", temp );
  
  sw.reset(); sw.start();
  for(int numReps=0; numReps<TIMING_REPS; numReps++)
    {
#endif
    
      translate();
      
#ifdef PRINT_TIMINGS
      
    }
  sw.stop();
  temp = (double)sw.timeMS()/TIMING_REPS;
  printf(" bf_translate: %7.2f\n", temp );
#endif
  
  if(num > 1 && !true)
    {
      unloadData(basisFunctionsFB, nCols*fxo_deltaBF, nRows*nMats*fxo_deltaOE);
      unloadData(outputFB, nCols*txt_deltaBF, nRows*nMats*txt_deltaOE);
    }
    
  return outputFB.getTextureID(0,0);
}

void GPUQMCBasisFunction::loadElectronPositions(Array1D<Array2D<double>*> &X, int start, int stop)
{
  int index, i, j;
  for(int c = 0; c < nCols; c++)
    {
      for(int r = 0; r < nRows; r++)
        {
          for(int electron=0; electron<nElectrons; electron++)
            {
              i = electron/4;
              j = electron%4;
              //      3*( (row selection)*(total width) + (column selection) );
              index = 3*( (r*elecH + i)*nCols*elecW + (c*elecW + j) );
              cpuData[index    ] = (GLfloat) X(c*nRows + r)->get(electron+start, 0);
              cpuData[index + 1] = (GLfloat) X(c*nRows + r)->get(electron+start, 1);
              cpuData[index + 2] = (GLfloat) X(c*nRows + r)->get(electron+start, 2);
            }
        }
    }
  glBindTexture(TEXTURE_TARGET, electronsTexID);
  glTexImage2D(TEXTURE_TARGET, 0, TEXTURE_INTERNAL_FORMAT,
               elecW*nCols, elecH*nRows, 0, GL_RGB, GL_FLOAT, cpuData);
}

void GPUQMCBasisFunction::translate()
{
  //now we have to translate those calculations from 4x1 format into 2x2 format
  int maxs = outputFB.getWidth();
  int maxt = outputFB.getHeight();
  outputFB.cleanBuffer(0);
  cgGLSetTextureParameter(fxo_to_txt_CGP, basisFunctionsFB.getTextureID(0,0));
  
  cgGLEnableProfile(g_cgProfile);// <-- this line MUST be within a Begin/End Capture pair
  
  cgGLEnableTextureParameter(fxo_to_txt_CGP);
  cgGLBindProgram(fxo_to_txt_CG);
  
  /*if there are a couple columns being calculated simultaneously, then they must
  stay separate. if there are an odd number of basisfunctions, then each successive
  column will start in the last pixel of the preceeding column. to prevent this, each
  column is given a new set of shifted texture coordinates to rasterize over.*/
  if(fxo_deltaBF%2==0)
    {
      drawPrimative(maxs,maxt,0,0);
    }
  else
    {
      glEnable(GL_SCISSOR_TEST);
      for(int col=0; col<nCols; col++)
        {
          glScissor( col*txt_deltaBF, 0, txt_deltaBF, maxt);
          drawPrimative(maxs,maxt,-0.5*col,0);
        }
      glDisable(GL_SCISSOR_TEST);
    }
  cgGLDisableProfile(g_cgProfile);
  
  cgGLDisableTextureParameter(fxo_to_txt_CGP);
  
  if(INT_FINISHES)
    {
      glFinish();
      glFlush();
    }
    
  getOpenGLError("Error in QMC basis function translation");
}

void GPUQMCBasisFunction::drawPrimative(GLfloat maxs, GLfloat maxt, GLfloat sShift, GLfloat tShift)
{
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

string GPUQMCBasisFunction::generateTranslationShader(bool is12)
{
  string testVarIsNaN;
  
  /*at one point in debugging, i was getting errors because of NaN's showing up
  in the calculated basis functions. this segment of code can remove them during 
  the translation process... but i'm not sure that it's necessary because all the
  basisfunction data that isn't supposed to be there will be multiplied by a zero
  in the coefficient matrix anyway. the question is whether NaN*0 = 0 the way i want
  it to be (if NaNs still do show up)*/
  testVarIsNaN +=
    "   VAR.x = isnan(VAR.x)? 0 : VAR.x;                                    \n"
    "   VAR.y = isnan(VAR.y)? 0 : VAR.y;                                    \n"
    "   VAR.z = isnan(VAR.z)? 0 : VAR.z;                                    \n"
    "   VAR.w = isnan(VAR.w)? 0 : VAR.w;                                    \n";
    
  string shader;
  shader +=
    "float4 main(in float2 coords : TEX0,                                   \n"
    "            uniform samplerRECT  input) : COLOR                        \n"
    "{                                                                      \n";
    
  /*if the number of electrons % 4 is 1 or 2, then that means that an extra
  row will be generated in the process of translating from 4x1 to 2x2. to 
  prevent this, each set of 5 rows is shifted up one pixel, or s pixels for
  the s-th row.*/
  if(is12)
    {
      shader +=
        "   int s = coords.y/ROW_HEIGHT;                                        \n"
        "   int t = fmod(coords.y,ROW_HEIGHT);                                  \n"
        "   coords.y += s;                                                      \n";
    }
    
  /*the multiply and divide by two manifest the switch from 4x1 to 2x2.
  originally, left and right were shifted by -1 and 0 respectively, but
  some of the texture coordinates were skipped that way.*/
  /* //stuff that didn't work
  //"   float switcher1 = frac(coords.y*0.5);                            \n"
  //"   int switcher = 2*frac(coords.y);                                  \n"
  //"   coords.x *= 2.0;                                                 \n"
  */
  shader +=
    "   float4 left  = texRECT(input, float2(2*coords.x-0.99, coords.y/2)); \n"
    "   float4 right = texRECT(input, float2(2*coords.x+0.01, coords.y/2)); \n";
    
  if(NAN_CHECKS)
    {
      shader += testVarIsNaN;
      findandreplace(shader,"VAR", "right");
      shader += testVarIsNaN;
      findandreplace(shader,"VAR", "left");
    }
    
  /*the switch here is to choose between whether we are looking at the top
  or bottom of a 2x2 pair originating from a (left and right) 4x1 pair.*/
  shader +=
    "   int switcher = fmod(coords.y,2);                                    \n"
    "   return (switcher==0  EXTRA_CHECK) ?                                 \n"
    "          float4(left.x, right.x, left.y, right.y):                    \n"
    "          float4(left.z, right.z, left.w, right.w);                    \n"
    "}                                                                      \n";
  findandreplace(shader,"ROW_HEIGHT", txt_deltaOE);
  if(is12)
    {
      findandreplace(shader,"EXTRA_CHECK", "|| t==0");
    }
  else
    {
      findandreplace(shader,"EXTRA_CHECK", "");
    }
  if(PRINT_SHADER)
    {
      cout << shader << endl;
      getchar();
    }
  return shader;
}

string GPUQMCBasisFunction::generateShader(int which)
{
  string shader;
  shader +=
    "float4 main(in float2 coords : TEX0,                                   \n"
    "            uniform samplerRECT  params,                               \n"
    "            uniform samplerRECT  epos) : COLOR                         \n"
    "{                                                                      \n"
    "   int bfX = fmod(coords.x,WIDTH);                                     \n"
    "   int eposX = coords.x/WIDTH;                                         \n"
    "   int eposY = coords.y;                                               \n"
    "   eposX *= 4;                                                         \n"
    "   float3 n_center = texRECT(params,float2(bfX,0)).xyz;                \n"
    "   float4 klm = texRECT(params,float2(bfX,1));                         \n"
    "   float ntexs = klm.w;                                                \n"
    "   float4 output = 0;                                                  \n"
    "   float4 r_sq = 0;                                                    \n"
    "   float4 xyz_term = 1;                                                \n"
    "   float3 r = 0;                                                       \n";
    
  if(which != psi)
    {
      shader +=
        "   float4 r_extra = 0;                                                 \n";
    }
    
  //the base is not allowed to be negative for a GPU pow command,
  //even if the exponent is an integer. This ?: mess is the best
  //alternative i've found...
  shader +=
    "                                                                       \n"
    "   r = texRECT(epos, float2(eposX,eposY)).xyz;                         \n"
    "   r = r + n_center;                                                   \n"
    "   r_sq.x = dot(r,r);                                                  \n"
    "   xyz_term.x = klm.x>0 ? xyz_term.x*r.x : xyz_term.x;                 \n"
    "   xyz_term.x = klm.x>1 ? xyz_term.x*r.x : xyz_term.x;                 \n"
    "   xyz_term.x = klm.y>0 ? xyz_term.x*r.y : xyz_term.x;                 \n"
    "   xyz_term.x = klm.y>1 ? xyz_term.x*r.y : xyz_term.x;                 \n"
    "   xyz_term.x = klm.z>0 ? xyz_term.x*r.z : xyz_term.x;                 \n"
    "   xyz_term.x = klm.z>1 ? xyz_term.x*r.z : xyz_term.x;                 \n"
    "   VARIABLE.x = EQUATION;                                              \n"
    "                                                                       \n"
    "   r = texRECT(epos, float2(eposX + 1,eposY)).xyz;                     \n"
    "   r = r + n_center;                                                   \n"
    "   r_sq.y = dot(r,r);                                                  \n"
    "   xyz_term.y = klm.x>0 ? xyz_term.y*r.x : xyz_term.y;                 \n"
    "   xyz_term.y = klm.x>1 ? xyz_term.y*r.x : xyz_term.y;                 \n"
    "   xyz_term.y = klm.y>0 ? xyz_term.y*r.y : xyz_term.y;                 \n"
    "   xyz_term.y = klm.y>1 ? xyz_term.y*r.y : xyz_term.y;                 \n"
    "   xyz_term.y = klm.z>0 ? xyz_term.y*r.z : xyz_term.y;                 \n"
    "   xyz_term.y = klm.z>1 ? xyz_term.y*r.z : xyz_term.y;                 \n"
    "   VARIABLE.y = EQUATION;                                              \n"
    "                                                                       \n"
    "   r = texRECT(epos, float2(eposX + 2,eposY)).xyz;                     \n"
    "   r = r + n_center;                                                   \n"
    "   r_sq.z = dot(r,r);                                                  \n"
    "   xyz_term.z = klm.x>0 ? xyz_term.z*r.x : xyz_term.z;                 \n"
    "   xyz_term.z = klm.x>1 ? xyz_term.z*r.x : xyz_term.z;                 \n"
    "   xyz_term.z = klm.y>0 ? xyz_term.z*r.y : xyz_term.z;                 \n"
    "   xyz_term.z = klm.y>1 ? xyz_term.z*r.y : xyz_term.z;                 \n"
    "   xyz_term.z = klm.z>0 ? xyz_term.z*r.z : xyz_term.z;                 \n"
    "   xyz_term.z = klm.z>1 ? xyz_term.z*r.z : xyz_term.z;                 \n"
    "   VARIABLE.z = EQUATION;                                              \n"
    "                                                                       \n"
    "   r = texRECT(epos, float2(eposX + 3,eposY)).xyz;                     \n"
    "   r = r + n_center;                                                   \n"
    "   r_sq.w = dot(r,r);                                                  \n"
    "   xyz_term.w = klm.x>0 ? xyz_term.w*r.x : xyz_term.w;                 \n"
    "   xyz_term.w = klm.x>1 ? xyz_term.w*r.x : xyz_term.w;                 \n"
    "   xyz_term.w = klm.y>0 ? xyz_term.w*r.y : xyz_term.w;                 \n"
    "   xyz_term.w = klm.y>1 ? xyz_term.w*r.y : xyz_term.w;                 \n"
    "   xyz_term.w = klm.z>0 ? xyz_term.w*r.z : xyz_term.w;                 \n"
    "   xyz_term.w = klm.z>1 ? xyz_term.w*r.z : xyz_term.w;                 \n"
    "   VARIABLE.w = EQUATION;                                              \n"
    "   xyz_term *= LARGE_MULTIPLIER;                                       \n"
    "                                                                       \n"
    "   float4 coeff = 0;                                                   \n"
    "   for(int j=0; j<ntexs; j++){                                         \n"
    "      coeff = texRECT(params,float2(bfX,2+j));                         \n"
    "      output += PREFACTOR1coeff.y*EXP_BASE(N_ONE_DIV_LN2*coeff.x*r_sq+SHIFT);    \n"
    "      if(coeff.z != 0){                                                \n"
    "      output += PREFACTOR2coeff.w*EXP_BASE(N_ONE_DIV_LN2*coeff.z*r_sq+SHIFT);    \n"
    "      }                                                                \n"
    "   }                                                                   \n"
    "   output *= xyz_term;                                                 \n"
    //*
    "   return output*EXP_SHIFT;                                            \n"
    /*/
    "   float send = 1e-35;                                                   \n"
    "   return float4(send,send,send,send);                                 \n"
    //"   return float4(5,5,5,5);                                 \n"
    //"   return texRECT(params,float2(bfX,3));\n"
    //"   return float4(texRECT(epos, float2(eposX,eposY)).xyz, eposY);\n"
    //"   return float4(texRECT(params,float2(bfX,1)).xyz, bfX);\n"
    //"   return float4(coords ,eposX, eposY);\n"
    //*/
    "}                                                                      \n";
    
  switch(which)
    {
        case psi:
        {
          findandreplace(shader,"VARIABLE", "//");
          findandreplace(shader,"PREFACTOR1", "");
          findandreplace(shader,"PREFACTOR2", "");
          break;
        }
        case grx:
        {
          findandreplace(shader,"VARIABLE", "r_extra");
          findandreplace(shader,"EQUATION", "r.x");
          findandreplace(shader,"PREFACTOR1", "(klm.x/r_extra - 2*coeff.x*r_extra)*");
          findandreplace(shader,"PREFACTOR2", "(klm.x/r_extra - 2*coeff.z*r_extra)*");
          break;
        }
        case gry:
        {
          findandreplace(shader,"VARIABLE", "r_extra");
          findandreplace(shader,"EQUATION", "r.y");
          findandreplace(shader,"PREFACTOR1", "(klm.y/r_extra - 2*coeff.x*r_extra)*");
          findandreplace(shader,"PREFACTOR2", "(klm.y/r_extra - 2*coeff.z*r_extra)*");
          break;
        }
        case grz:
        {
          findandreplace(shader,"VARIABLE", "r_extra");
          findandreplace(shader,"EQUATION", "r.z");
          findandreplace(shader,"PREFACTOR1", "(klm.z/r_extra - 2*coeff.x*r_extra)*");
          findandreplace(shader,"PREFACTOR2", "(klm.z/r_extra - 2*coeff.z*r_extra)*");
          break;
        }
        case lap:
        {
          findandreplace(shader,"VARIABLE", "r_extra");
          findandreplace(shader,"EQUATION", "dot( klm.xyz/r, (klm.xyz-1.0)/r )");
          //(a*(a-1.0)/x2 + b*(b-1.0)/y2 + c*(c-1.0)/z2 - (4.0*(a+b+c) + 6.0 - 4.0*r_sq*p0)*p0)*exp_term;
          findandreplace(shader,"PREFACTOR1",
                         "(r_extra + (-4.0*(klm.x + klm.y + klm.z - r_sq*coeff.x) - 6.0)*coeff.x)*");
          findandreplace(shader,"PREFACTOR2",
                         "(r_extra + (-4.0*(klm.x + klm.y + klm.z - r_sq*coeff.z) - 6.0)*coeff.z)*");
          break;
        }
    }
  findandreplace(shader,"WIDTH",fxo_deltaBF);
  findandreplace(shader,"HEIGHT",fxo_deltaOE);
  findandreplace(shader,"NELECTRONS",nElectrons);
  findandreplace(shader,"LARGE_MULTIPLIER","exp2((float)MULTIPLIER)");
  findandreplace(shader,"MULTIPLIER",MULTIPLIER);
  findandreplace(shader,"EXP_SHIFT","EXP_BASE((float)-SHIFT)");
  findandreplace(shader,"SHIFT",EXP_SHIFT);
  
  if(USE_BASE_2)
    {
      findandreplace(shader,"N_ONE_DIV_LN2","-1.0/log(2.0)");
      findandreplace(shader,"EXP_BASE","exp2");
    }
  else
    {
      findandreplace(shader,"N_ONE_DIV_LN2","-1.0");
      findandreplace(shader,"EXP_BASE","exp");
    }
    
  if(PRINT_SHADER)
    {
      cout << "which: " << which << endl << shader << endl;
      getchar();
    }
  return shader;
}

void GPUQMCBasisFunction::unloadData(GPUQMCFramebuffer & fb, int w, int h)
{
  fb.readFrom(0,0);
  glReadPixels(0,0,w,h,GL_RGBA,GL_FLOAT,cpuData);
  cout << "unloaded by bf\n";
  PrintRGBAPixelsBoxE(cpuData,w,h,20,20,-1,-1,true);
}

int GPUQMCBasisFunction::mapping(int i, int j, int h, int w)
{
  return 4*(i*w + j);
}

void GPUQMCBasisFunction::setUpInputs()
{
  int nGaussians, bf = 0, index;
  int count1=0, count2=0;
  for (int atom=0; atom<flags->Natoms; atom++)
    {
      for (int j=0; j<BFCoeffs(atom).getNumberBasisFunctions(); j++)
        {
        
          //Rc is in the first row. the negation is because the shader can add faster than it can subtract
          index = mapping(0,bf,basisfunctionParamsH,fxo_deltaBF);
          for (int translate=0; translate<3; translate++)
            cpuData[index + translate] = -1.0*Molecule->Atom_Positions(atom,translate);
            
          nGaussians = BFCoeffs(atom).N_Gauss(j);
          
          //k, l, m are in the second row
          index = mapping(1,bf,basisfunctionParamsH,fxo_deltaBF);
          for (int translate=0; translate<3; translate++)
            cpuData[index + translate] = BFCoeffs(atom).xyz_powers(j,translate);
          //cpuData[index + 3] = BFCoeffs(atom).xyz_powers(j,0) + BFCoeffs(atom).xyz_powers(j,1) + BFCoeffs(atom).xyz_powers(j,2);
          //the number of textures used for coefficients is also in the second row
          cpuData[index + 3] = (int)(nGaussians/2);
          if(nGaussians%2 != 0) cpuData[index + 3]++;
          
          //the gaussian coefficients
          for (int i=0; i<nGaussians; i++)
            {
              if(2 + i/2 >= basisfunctionParamsH) cout << "ERROR: cpuData not big enough for all coefficients\n";
              index = mapping(2 + i/2,bf,basisfunctionParamsH,fxo_deltaBF);
              int add1 = (i%2)?2:0;//no idea why i couldn't just stick these between the braces...
              int add2 = (i%2)?3:1;
              cpuData[index + add1] = BFCoeffs(atom).Coeffs.array()[j][i][0];//this becomes coeff.x or coeff.z
              cpuData[index + add2] = BFCoeffs(atom).Coeffs.array()[j][i][1];//this becomes coeff.y or coeff.w
            }
          bf++;
        }
    }
    
  glBindTexture(TEXTURE_TARGET, bfParametersTexID);
  glTexImage2D(TEXTURE_TARGET, 0, TEXTURE_INTERNAL_FORMAT,
               fxo_deltaBF, basisfunctionParamsH, 0, GL_RGBA, GL_FLOAT, cpuData);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  
  glBindTexture(GL_TEXTURE_RECTANGLE_NV, electronsTexID);
  glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

int GPUQMCBasisFunction::getNumIterations()
{
#ifdef PRINT_TIMINGS
  return TIMING_REPS;
#else
return 1;
#endif
}

void GPUQMCBasisFunction::operator=(GPUQMCBasisFunction & rhs)
{
  nRows = rhs.nRows;
  nCols = rhs.nCols;
  nElectrons = rhs.nElectrons;
  nBasisF = rhs.nBasisF;
  allocatedRows = rhs.allocatedRows;
  allocatedCols = rhs.allocatedCols;
  fxo_deltaOE = rhs.fxo_deltaOE;
  fxo_deltaBF = rhs.fxo_deltaBF;
  txt_deltaBF = rhs.txt_deltaBF;
  txt_deltaOE = rhs.txt_deltaOE;
  elecW = rhs.elecW;
  elecH = rhs.elecH;
  maxGaussians = rhs.maxGaussians;
  basisfunctionParamsH = rhs.basisfunctionParamsH;
  
  fragProg = rhs.fragProg;
  electronsCGP = rhs.electronsCGP;
  paramsCGP = rhs.paramsCGP;
  fxo_to_txt_CG = rhs.fxo_to_txt_CG;
  fxo_to_txt_CGP = rhs.fxo_to_txt_CGP;
  
  basisFunctionsFB = rhs.basisFunctionsFB;
  outputFB = rhs.outputFB;
  
  cpuData = (GLfloat *) calloc( nCols*fxo_deltaBF * nRows*nMats*max(fxo_deltaOE,basisfunctionParamsH) * 4 , sizeof(GLfloat) );
}

GPUQMCBasisFunction::GPUQMCBasisFunction()
{
  nRows = 0; nCols = 0;
  cpuData = (GLfloat *) calloc( 1 , sizeof(GLfloat) );
}

#endif
