/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "GPUQMCFramebuffer.h"
#ifdef QMC_GPU

#define TEXTURE_INTERNAL_FORMAT    GL_FLOAT_RGBA32_NV
#define TEXTURE_TARGET             GL_TEXTURE_RECTANGLE_NV

GPUQMCFramebuffer::GPUQMCFramebuffer()
{
  w = 0; h = 0;
}

GPUQMCFramebuffer::GPUQMCFramebuffer(int width, int height, int numFB, int numRT)
{
  initialize(width,height,numFB, numRT);
}

GPUQMCFramebuffer::~GPUQMCFramebuffer()
{
  for(int i=0; i<getNumFB(); i++)
    {
      for(int j=0; j<getNumRT(); j++)
        {
          glDeleteTextures(1, &attachedTextures[i][j]);
        }
      glDeleteFramebuffersEXT(1,&theFrameBuffers[i]);
    }
  theFrameBuffers.clear();
  attachedTextures.clear();
  bufferAttachments.clear();
}

void GPUQMCFramebuffer::operator=(GPUQMCFramebuffer & rhs)
{
  initialize(rhs.w, rhs.h, rhs.getNumFB(), rhs.getNumRT());
}

void GPUQMCFramebuffer::initialize(int width, int height, int numFB, int numRT)
{
  w = width;
  h = height;
  
  if(numRT < 1 || numRT > 4)
    cerr << "Error: number render targets must be between 1 and 4\n";
  if(width*height == 0)
    cerr << "Error: framebuffer must have non-zero dimensions\n";
    
  theFrameBuffers.clear();
  attachedTextures.clear();
  bufferAttachments.clear();
  theFrameBuffers.resize(numFB);
  attachedTextures.resize(numFB);
  bufferAttachments.resize(numFB);
  
  for(int i=0; i<getNumFB(); i++)
    {
      glGenFramebuffersEXT(1, &theFrameBuffers[i]);
      bufferAttachments[i] = new GLuint[4];
    }
    
  for(int i=0; i<numRT; i++)
    addTexture();
    
  checkFramebufferStatus();
}

int GPUQMCFramebuffer::getWidth()
{
  return w;
}

int GPUQMCFramebuffer::getHeight()
{
  return h;
}

void GPUQMCFramebuffer::drawTo(int whichFB)
{
  if(!testExistance(whichFB,0)) return;
  
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, theFrameBuffers[whichFB]);
  
  if(getNumFB() == 1)
    {
      glDrawBuffer(bufferAttachments[whichFB][0]);
    }
  else
    {
      glDrawBuffers(getNumRT(), bufferAttachments[whichFB]);
    }
    
  glViewport( 0, 0, getWidth(), getHeight());
  GET_GLERROR(0);
}

void GPUQMCFramebuffer::readFrom(int whichFB, int whichRT)
{
  if(!testExistance(whichFB,whichRT)) return;
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, theFrameBuffers[whichFB]);
  glReadBuffer(bufferAttachments[whichFB][whichRT]);
  GET_GLERROR(0);
}

void GPUQMCFramebuffer::cleanBuffer(int whichFB)
{
  if(!testExistance(whichFB,0)) return;
  drawTo(whichFB);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void GPUQMCFramebuffer::cleanAllBuffers()
{
  for(int i=0; i<getNumFB(); i++)
    cleanBuffer(i);
}

void GPUQMCFramebuffer::swapBuffers()
{
  cerr << "Error: Function " << __FUNCTION__ << " is not implemented yet!\n";
}

void GPUQMCFramebuffer::setupSwapBuffers(int first, int second)
{
  cerr << "Error: Function " << __FUNCTION__ << " is not implemented yet!\n";
}

GLuint GPUQMCFramebuffer::addTexture()
{
  int numFB = getNumFB();
  int numRT = getNumRT();
  
  if(numRT >= 4)
    {
      cerr << "ERROR: too many buffers requested\n";
      return 0;
    }
    
  GLuint * newTextureIDs = new GLuint[numFB];
  glGenTextures(numFB, newTextureIDs);
  
  //add 1 texture per framebuffer
  for(int i=0; i<numFB; i++)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, theFrameBuffers[i]);
      glBindTexture(TEXTURE_TARGET, newTextureIDs[i]);
      glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf(TEXTURE_TARGET, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexImage2D(TEXTURE_TARGET, 0, TEXTURE_INTERNAL_FORMAT,
                   w, h, 0, GL_RGBA, GL_FLOAT, NULL);
                   
      bufferAttachments[i][numRT] = GL_COLOR_ATTACHMENT0_EXT + numRT;
      
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                bufferAttachments[i][numRT], TEXTURE_TARGET, newTextureIDs[i], 0);
                                
      attachedTextures[i].push_back(newTextureIDs[i]);
    }
    
  return 0;
}

void GPUQMCFramebuffer::addTexture(GLuint newTexID)
{
  cerr << "Error: Function " << __FUNCTION__ << " is not implemented yet!\n";
}

void GPUQMCFramebuffer::removeTexture(int whichRT)
{
  cerr << "Error: Function " << __FUNCTION__ << " is not implemented yet!\n";
}

GLuint GPUQMCFramebuffer::getTextureID(int whichFB, int whichRT)
{
  if(!testExistance(whichFB,whichRT)) return 0;
  return attachedTextures[whichFB][whichRT];
}

int GPUQMCFramebuffer::getNumRT()
{
  return (int)attachedTextures[0].size();
}

int GPUQMCFramebuffer::getNumFB()
{
  return (int)theFrameBuffers.size();
}

bool GPUQMCFramebuffer::testExistance(int whichFB, int whichRT)
{
  if(whichRT >= getNumRT() || whichRT < 0)
    {
      cerr << "ERROR: render target " << whichRT << " does not exist.\n";
      return false;
    }
    
  if(whichFB >= getNumFB() || whichFB < 0)
    {
      cerr << "ERROR: framebuffer " << whichFB << " does not exist.\n";
      return false;
    }
  return true;
}

void GPUQMCFramebuffer::checkFramebufferStatus()
{
  GLenum status;
  status = (GLenum) glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch(status)
    {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
        {
          break;
        }
        case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
        {
          printf("Unsupported framebuffer format\n");
          break;
        }
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
        {
          printf("Framebuffer incomplete, missing attachment\n");
          break;
        }
        case GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT:
        {
          printf("Framebuffer incomplete, duplicate attachment\n");
          break;
        }
        case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
        {
          printf("Framebuffer incomplete, attached images must have same dimensions\n");
          break;
        }
        case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
        {
          printf("Framebuffer incomplete, attached images must have same format\n");
          break;
        }
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
        {
          printf("Framebuffer incomplete, missing draw buffer\n");
          break;
        }
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
        {
          printf("Framebuffer incomplete, missing read buffer\n");
          break;
        }
        default:
        {
          cerr << "Error: Framebuffer screwed up for some odd reason...\n";
          exit(-1);
        }
    }
}

#endif
