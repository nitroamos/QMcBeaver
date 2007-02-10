/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "GPUGlobals.h"

#ifdef QMC_GPU

void GPUGlobals::getOpenGLError(char *msg)
{
  static char error_txt[][32] = {
                                  "GL_INVALID_ENUM",
                                  "GL_INVALID_VALUE",
                                  "GL_INVALID_OPERATION",
                                  "GL_STACK_OVERFLOW",
                                  "GL_STACK_UNDERFLOW",
                                  "GL_OUT_OF_MEMORY" };
  GLenum e = glGetError();
  if (e != GL_NO_ERROR)
    {
      cerr << "Error: OpenGL error " << error_txt[e-GL_INVALID_ENUM] << endl;
      if(msg != NULL) cerr << "Message: " << msg << endl;
    }
}

void GPUGlobals::PrintRGBAPixelsBoxF(float* pix, int w, int h)
{
  int index;
  int maxJ = 30;
  for (int i = 0; i < h; i++)
    {
      for (int j = 0; j < w && j < maxJ; j++)
        {
          index = 4*(i*w + j);
          printf("%6.3f",pix[index +0]);
          printf("%6.3f",pix[index +1]);
          printf("   ");
        }
      printf("\n");
      for (int j = 0; j < w && j < maxJ; j++)
        {
          index = 4*(i*w + j);
          printf("%6.3f",pix[index +2]);
          printf("%6.3f",pix[index +3]);
          printf("   ");
        }
      printf("\n\n");
    }
  printf("\n");
}

void GPUGlobals::PrintRGBAPixelsBoxE(float* pix, int w, int h,
                                     int wBorder, int hBorder,
                                     int maxJ, int maxI,
                                     bool printRC)
{
  int index;
  
  //cause the windows console is not "infinitely" wide like it is high.
  //meant to prevent the console from wrapping the results and making
  //the output look ugly.
  if(maxJ < 0 || maxJ > w) maxJ = w;
  if(maxI < 0 || maxI > h) maxI = h;
  
  if(wBorder > w) wBorder = -1;
  if(hBorder > h) hBorder = -1;

  int wAdjuster = 0;
  if(wBorder < 0 && maxJ%2 == 1)
    wAdjuster = 1;
  int hAdjuster = 0;
  if(hBorder < 0 && maxI%2 == 1)
    hAdjuster = 1;
    
  if(wBorder < 0 || wBorder*2 > maxJ) wBorder = (int)(maxJ/2.0);
  if(hBorder < 0 || hBorder*2 > maxI) hBorder = (int)(maxI/2.0);
  
  const char * form = "%11.3e";
  
  //this will label the rows and columns
  if(printRC)
    {
      printf("%16i",0);
      for (int j=1; j < wBorder + wAdjuster; j++)
        printf("%25i",j);
      for (int j=w-wBorder; j < w; j++)
        printf("%25i",j);
      printf("\n\n");
    }
    
  //top chunk with thickness of hBorder
  for (int i = 0; i < hBorder + hAdjuster; i++)
    {
      if(printRC)
        {
          printf("%-4i",i);
        }
      //pixels line one
      //left chunk with thickness of wBorder
      for (int j = 0; j < wBorder + wAdjuster; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +0]);
          printf(form,pix[index +1]);
          printf("   ");
        }
      //right chunk with thickness of wBorder
      for (int j = w-wBorder; j < w; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +0]);
          printf(form,pix[index +1]);
          printf("   ");
        }
        
      printf("\n");
      if(printRC)
        {
          printf("    ");
        }
        
      //pixels line two
      //left chunk with thickness of wBorder
      for (int j = 0; j < wBorder + wAdjuster; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +2]);
          printf(form,pix[index +3]);
          printf("   ");
        }
      //right chunk with thickness of wBorder
      for (int j = w-wBorder; j < w; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +2]);
          printf(form,pix[index +3]);
          printf("   ");
        }
        
      //partition between a row of pixels
      printf("\n\n");
    }
    
  //bottom chunk with thickness of hBorder
  for (int i = h-hBorder; i < h; i++)
    {
      if(printRC)
        {
          printf("%-4i",i);
        }
      //pixels line one
      //left chunk with thickness of wBorder
      for (int j = 0; j < wBorder + wAdjuster; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +0]);
          printf(form,pix[index +1]);
          printf("   ");
        }
      //right chunk with thickness of wBorder
      for (int j = w-wBorder; j < w; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +0]);
          printf(form,pix[index +1]);
          printf("   ");
        }
        
      printf("\n");
      if(printRC)
        {
          printf("    ");
        }
        
      //pixels line two
      //left chunk with thickness of wBorder
      for (int j = 0; j < wBorder + wAdjuster; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +2]);
          printf(form,pix[index +3]);
          printf("   ");
        }
      //right chunk with thickness of wBorder
      for (int j = w-wBorder; j < w; j++)
        {
          index = 4*(i*w + j);
          printf(form,pix[index +2]);
          printf(form,pix[index +3]);
          printf("   ");
        }
        
      //partition between a row of pixels
      printf("\n\n");
    }
    
  printf("\n");
  
}

void GPUGlobals::PrintRGBAPixelsColumn(float* pix, int w, int h)
{
  for (int i = 0; i < 4*h; i++)
    {
      for (int j = 0; j < w && j < 26; j++)
        {//> around 26 will wrap
          printf("%7.3f",pix[4*((h-(i/4)-1)*w + j)+i%4]);
          //the %12.8f will show pretty much the max precision
          //printf("%12.8f",pix[4*((h-(i/4)-1)*w + j)+i%4]);
        }
      printf("\n");
      if(i%4==3) printf("\n");
    }
  printf("\n");
}

void GPUGlobals::PrintRGBPixelsColumn(float* pix, int w, int h)
{
  int maxJ = 26;
  for (int i = 0; i < 3*h; i++)
    {
      for (int j = 0; j < w && j < maxJ; j++)
        {
          printf("%7.3g",pix[3*(( (i/3) )*w + j)+i%3]);
          //the %12.8f will show pretty much the max precision
          //printf("%12.8f",pix[4*((h-(i/4)-1)*w + j)+i%4]);
        }
      printf("\n");
      if(i%3==2) printf("\n");
    }
  printf("\n");
}

void GPUGlobals::PrintMatrix(Array2D<GLfloat> matrix)
{
  for(int i=0; i<matrix.dim1(); i++)
    {
      for(int j=0; j<matrix.dim2() && j < 28; j++)
        {
          printf("%7.3g", (float)matrix(i,j));
        }
      printf("\n");
    }
  printf("\n");
}

void GPUGlobals::getFactors(int num, int & fact1, int & fact2)
{
  //for hmx, cols (fact2) can be as high as 7, rows (fact1) can be as high as 21
  int max = (int) sqrt( (double)num );
  fact1 = num;
  fact2 = 1;
  
  for(int i=2; i<=max; i++)
    if( num%i == 0 )
      {
        fact1 = num/i;
        fact2 = i;// i <= num/i
      }
}

const char * GPUGlobals::gpubench_getDriverVersion(void)
{
  int infosize;
  void *info;
  int status;
  int vsinfosize;
  char *desc;
  LPVOID version=NULL;
  DWORD vLen,langD;
  static char fileVersion[256];
  static char returnValue[256];
  char *p;
  
  static char nvdriver[] = "c:\\windows\\system32\\nv4_disp.dll";
  static char atidriver[] = "c:\\windows\\system32\\atioglxx.dll";
  static char driver3dl[] = "c:\\windows\\system32\\3Dl2DD.dll";
  
  char *filename;
  
  if (gpu_make == ati_gpu) filename = atidriver;
  else if (gpu_make == nvidia_gpu) filename = nvdriver;
  else if (gpu_make == labs3d_gpu) filename = driver3dl;
  else return returnValue;
  
  //sprintf (fileVersion, "");
  
  infosize = GetFileVersionInfoSize( filename, (LPDWORD) &infosize );
  if (!infosize) return "";
  
  info = malloc (infosize);
  
  status = GetFileVersionInfo(filename, 0, infosize, info);
  if (!status)
    {
      free(info); return "";
    }
    
  status = VerQueryValue ( info, TEXT("\\VarFileInfo\\Translation"),
                           (LPVOID *) &version, (UINT *)&vLen);
  if (!status)
    {
      free(info); return "";
    }
    
  memcpy(&langD,version,4);
  sprintf(fileVersion, "\\StringFileInfo\\%02X%02X%02X%02X\\FileVersion",
          (langD & 0xff00)>>8,langD & 0xff,(langD & 0xff000000)>>24,
          (langD & 0xff0000)>>16);
          
  if (VerQueryValue(info, fileVersion, (LPVOID *)&desc, (PUINT) &vsinfosize))
    {
      p = desc;
      if (gpu_make != labs3d_gpu)
        {
          for  (   ; *p && *p != '.'; p++);
          for  (p++; *p && *p != '.'; p++);
          for  (p++; *p && *p != '.'; p++);
          p++;
        }
      else
        {
          while (*p != '\0' && *p != ' ') p++;
          if (*p == ' ') *p = '\0';
          p = desc;
        }
    }
  else
    {
      p = "";
    }
  strncpy(returnValue, p, sizeof returnValue);
  returnValue[sizeof returnValue - 1] = '\0';
  
  free (info);
  return returnValue;
}

void GPUGlobals::printVersions(ostream & strm)
{
  const char* vendor = (const char *) glGetString(GL_VENDOR);
  const char* render = (const char *) glGetString(GL_RENDERER);
  printf("%s\n", render);
  
  if (vendor == NULL)
    {
      cout << "NULL vendor string.\n";
    }
    
  if (strstr(vendor, "NVIDIA") != NULL)
    {
      cout << "nVidia";
      gpu_make = nvidia_gpu;
    }
  else if (strstr(vendor, "ATI") != NULL)
    {
      cout << "ATI";
      gpu_make = ati_gpu;
    }
  else if (strstr(vendor, "3Dlabs") != NULL)
    {
      cout << "3Dlabs";
      gpu_make = labs3d_gpu;
    }
  else
    {
      printf("Device Unknown: '%s'\n", vendor);
      gpu_make = unknown_gpu;
    }
    
  strm << " " << gpubench_getDriverVersion() << ", Cg Version: " << CG_VERSION_NUM;
  if(g_cgProfile == CG_PROFILE_FP40) strm << ", Using fp40\n";
  if(g_cgProfile == CG_PROFILE_FP30) strm << ", Using fp30\n";
  if(g_cgProfile == CG_PROFILE_FP20) strm << ", Using fp20\n";
}

void GPUGlobals::writeShader(const char * shader, const char * filename)
{
  ofstream shaderFile(filename);
  if(shaderFile.is_open())
    {
      shaderFile << shader;
      shaderFile.close();
    }
  else
    {
      cerr << "ERROR: " << filename << " failed to open\n";
    }
}

void GPUGlobals::checkExtensions()
{
  bool ok = true;
  
  if(!GL_ARB_fragment_program)
    {
      cout << "GL_FRAGMENT_PROGRAM_ARB not supported" << endl;
      ok = false;
    }
    
  if(!GL_NV_fragment_program)
    {
      cout << "GL_NV_fragment_program not supported" << endl;
      ok = false;
    }
    
  if(!GL_NV_fragment_program2)
    {
      cout << "GL_NV_fragment_program2 not supported" << endl;
      ok = false;
    }
    
  if(!GL_NV_float_buffer)
    {
      cout << "GL_NV_float_buffer not supported" << endl;
      ok = false;
    }
    
  if(!GL_EXT_framebuffer_object)
    {
      cout << "GL_EXT_framebuffer_object not supported" << endl;
      ok = false;
    }
    
  if(!GL_ARB_texture_rectangle)
    {
      cout << "ARB_texture_rectangle not supported" << endl;
      ok = false;
    }
    
  if(!ok)
    {
      cerr << "ERROR: required extensions are not supported, exiting.\n";
      exit(-1);
    }
}


void GPUGlobals::specTest()
{
  //GLenum target = GL_NV_fragment_program2;
  //GLenum target = GL_FRAGMENT_PROGRAM_NV;
  GLenum target = GL_FRAGMENT_PROGRAM_ARB;
  GLint num, nativeNum;
  glGetProgramivARB( target, GL_PROGRAM_INSTRUCTIONS_ARB, &num );
  glGetProgramivARB( target, GL_PROGRAM_NATIVE_INSTRUCTIONS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_PROGRAM_INSTRUCTIONS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_PROGRAM_TEMPORARIES_ARB, &num );
  glGetProgramivARB( target, GL_PROGRAM_NATIVE_TEMPORARIES_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_PROGRAM_TEMPORARIES_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_PROGRAM_PARAMETERS_ARB, &num );
  glGetProgramivARB( target, GL_PROGRAM_NATIVE_PARAMETERS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_PROGRAM_PARAMETERS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_PROGRAM_ATTRIBS_ARB, &num );
  glGetProgramivARB( target, GL_PROGRAM_NATIVE_ATTRIBS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_PROGRAM_ATTRIBS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_PROGRAM_ALU_INSTRUCTIONS_ARB, &num );
  glGetProgramivARB( target, GL_PROGRAM_NATIVE_ALU_INSTRUCTIONS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_PROGRAM_ALU_INSTRUCTIONS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_PROGRAM_TEX_INSTRUCTIONS_ARB, &num );
  glGetProgramivARB( target, GL_PROGRAM_NATIVE_TEX_INSTRUCTIONS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_PROGRAM_TEX_INSTRUCTIONS_ARB",num,nativeNum);
  cout << endl;
  
  glGetProgramivARB( target, GL_MAX_PROGRAM_INSTRUCTIONS_ARB, &num );
  glGetProgramivARB( target, GL_MAX_PROGRAM_NATIVE_INSTRUCTIONS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_MAX_PROGRAM_INSTRUCTIONS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_MAX_PROGRAM_TEMPORARIES_ARB, &num );
  glGetProgramivARB( target, GL_MAX_PROGRAM_NATIVE_TEMPORARIES_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_MAX_PROGRAM_TEMPORARIES_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_MAX_PROGRAM_PARAMETERS_ARB, &num );
  glGetProgramivARB( target, GL_MAX_PROGRAM_NATIVE_PARAMETERS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_MAX_PROGRAM_PARAMETERS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_MAX_PROGRAM_ATTRIBS_ARB, &num );
  glGetProgramivARB( target, GL_MAX_PROGRAM_NATIVE_ATTRIBS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_MAX_PROGRAM_ATTRIBS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_MAX_PROGRAM_ALU_INSTRUCTIONS_ARB, &num );
  glGetProgramivARB( target, GL_MAX_PROGRAM_NATIVE_ALU_INSTRUCTIONS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_MAX_PROGRAM_ALU_INSTRUCTIONS_ARB",num,nativeNum);
  glGetProgramivARB( target, GL_MAX_PROGRAM_TEX_INSTRUCTIONS_ARB, &num );
  glGetProgramivARB( target, GL_MAX_PROGRAM_NATIVE_TEX_INSTRUCTIONS_ARB, &nativeNum );
  printf("%-40s is (num, native): %2i, %2i\n", "GL_MAX_PROGRAM_TEX_INSTRUCTIONS_ARB",num,nativeNum);
  glGetProgramivARB(target, GL_MAX_PROGRAM_LOOP_COUNT_NV, &num);
  printf("%-40s is: %2i\n", "GL_MAX_PROGRAM_LOOP_COUNT_NV",num);
  glGetProgramivARB(target, GL_MAX_PROGRAM_LOOP_DEPTH_NV, &num);
  printf("%-40s is: %2i\n", "GL_MAX_PROGRAM_LOOP_DEPTH_NV",num);
  glGetProgramivARB(target, GL_MAX_PROGRAM_IF_DEPTH_NV, &num);
  printf("%-40s is: %2i\n", "GL_MAX_PROGRAM_IF_DEPTH_NV",num);
  
  getOpenGLError("Error in glGetProgramivARB");
  cout << endl;
  
  /*
  glGetIntegerv(GL_MAX_FRAGMENT_PROGRAM_LOCAL_PARAMETERS_NV, &num );
  cout << "GL_MAX_FRAGMENT_PROGRAM_LOCAL_PARAMETERS_NV is " << num << endl;
  glGetIntegerv(GL_MAX_PROGRAM_EXEC_INSTRUCTIONS_NV, &num );
  cout << "GL_MAX_PROGRAM_EXEC_INSTRUCTIONS_NV is " << num << endl;
  cout << endl;
  getOpenGLError("Error in glGetIntegerv");
  */
}

#endif

