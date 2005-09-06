#ifndef GPU_FRAMEBUFFER
#define GPU_FRAMEBUFFER

#ifdef QMC_GPU

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <string>

#include "GPUGlobals.h"

/**
 This class is meant to provide an easy to understand and use wrapper for the new
 EXT_framebuffer_object OpenGL extension. To read about how the specific functions
 work, read the registry entry:
 http://oss.sgi.com/projects/ogl-sample/registry/EXT/framebuffer_object.txt
 
 To summarize the registry entry as I understand it, a framebuffer object (FBO)
 is meant to provide a way for OpenGL itself to create off screen buffers to
 graphics programs. Previously, this required using OS dependent windowing calls
 (glX, wgl, agl, etc) so this makes it easy. The advantage to using this over
 the OS dependent stuff (like the RenderTexture class) is two fold:
 1) portability
 2) ease of turning a buffer output into input
 
 using windowing commands, if you wanted to use an off screen graphics buffer 
 for input to a later function, then you had to perform a copy operation. this
 extension removes that requirement. This makes multipass algorithms (where a
 particular texture's function is switched between output and input multiple times)
 faster.
 
 this class creates two things: framebuffers and rendertextures.
 up to four rendertextures can be held by a single framebuffer, and thus is useful
 in multple render targets type shaders.
*/
class GPUQMCFramebuffer
  {
    public:
      GPUQMCFramebuffer();
      
      /**
       The basic constructor.
       @param width the dimensions must be the same for all render targets in order to be
        bound to the same shader.
       @param height
       @param numFB the number of framebuffers. e.g. 2 for multipass.
       @param numRT the number of textures to add to each framebuffer.
      */
      GPUQMCFramebuffer(int width, int height, int numFB, int numRT);
      
      /**
       I'm not sure what it would mean to create a copy of a class containing
       textures... so this simply creates an object with the same set up, just
       different textures.
      */
      void operator=(GPUQMCFramebuffer & rhs);
      
      /**
       Frees all the GPU resources  
      */
      ~GPUQMCFramebuffer();
      
      /**
       To be functional, this function must be called, either by the constructor or
       by the user.
       @param width the dimensions must be the same for all render targets in order to be
        bound to the same shader.
       @param height
       @param numFB the number of framebuffers. e.g. 2 for multipass.
       @param numRT the number of textures to add to each framebuffer.
      */
      void initialize(int width, int height, int numFB, int numRT);
      
      /**
       @return the width of the framebuffer -- is the same for all attached
        textures.
      */
      int getWidth();
      
      /**
       @return the height of the framebuffer -- is the same for all attached
        textures.
      */
      int getHeight();
      
      /**
       This function binds the framebuffer to the output, selecting
      
       It is necessary to call glViewport(0,0,w,h) before rendering
       to the FBO for the following explanation
      
       (Cliff Woolley at GPGPU forums): 
       This is one case where pbuffers differ from FBOs, in the sense that
       FBOs do not have their own GL context or window, whereas pbuffers do.
       A pbuffer therefore has its own viewport which is initialized to
       the size of the pbuffer, whereas an FBO uses the same viewport as
       the rest of the GL context in which it is used (probably that of the
       main window of the application). If the FBO is of different dimensions
       than the window, you'll have to change the viewport to match it
       when you bind the FBO and change it back when you rebind the
       default framebuffer.
      
       I figure that if glViewport needs to be called before any rendering,
       then that might as well be handled by GPUQMCFramebuffer.
      
       @param which is the index of the framebuffer to draw to. if there are
        multiple textures bound to the framebuffer, then all of them will be
        sent to glDrawBuffers.
      */
      void drawTo(int whichFB);
      
      /**
       This function sets the caller up to run glReadPixels with desired
       dimensions.
      
       @param whichFB is the framebuffer index to read from
       @param whichRT is the render target index to read from
      */
      void readFrom(int whichFB, int whichRT);
      
      /**
       This function draws {0,0,0,0} (black) into all pixels associated
       with all the textures for a framebuffer. This also calls drawTo(@param whichFB),
       so if you're going to call drawTo(@param whichFB) right after anyway,
       this makes it unnecessary.
      
       @param which is the framebuffer index to clear
      */
      void cleanBuffer(int whichFB);
      
      /**
       This will iterate through all textures associated with this GPUQMCFramebuffer
       and write {0,0,0,0} (black) into them all.
      */
      void cleanAllBuffers();
      
      /**
       This was intended to make it easier to do multipass, but it never got implemented...
       Still might do it...
      */
      void swapBuffers();
      
      /**
       This was intended to make it easier to do multipass, but it never got implemented...
      */
      void setupSwapBuffers(int first, int second);
      
      /**
       This function will check for errors in the framebuffer.
       It was taken from the simple_framebuffer_object.cpp
       sample code for the GL_EXT_framebuffer_object extension
      
       The extension file describes an OpenGL function of the same
       name... perhaps I should just get rid of this.
      */
      void checkFramebufferStatus();
      
      /**
       This function will create a new texture for each framebuffer
       and add it to the FBO.
      */
      GLuint addTexture();
      
      /**
       This seemed like it would be useful.. But I haven't found that use yet...
       so it's not implemented.
      */
      void addTexture(GLuint newTexID);
      
      /**
       This seemed like it would be useful.. But I haven't found that use yet...
       so it's not implemented.
      */
      void removeTexture(int whichRT);
      
      /**
       This function will return the ID of a specific texture.
       @param whichFB which framebuffer
       @param whichRT which render target
      */
      GLuint getTextureID(int whichFB, int whichRT);
      
      /**
       Each framebuffer in a GPUQMCFramebuffer has the same number of 
       textures.
       @return the number of rendertargets
      */
      int getNumRT();
      
      /**
       @return the number of framebuffers this GPUQMCFramebuffer has
      */
      int getNumFB();
      
      /**
       This function provides a way for any of the other functions
       to check whether their inputs are valid.
       zero is allowed in order to check for the validity of the other
       @param whichFB does this index have a framebuffer?
       @param whichRT do the framebuffers have at least this many textures?
       @return false means it's bad.
      */
      bool testExistance(int whichFB, int whichRT);
      
    private:
    
      /**
       Indexed by the number of framebuffers this GPUQMCFramebuffer has
      */
      vector< GLuint >  theFrameBuffers;
      
      /**
       Indexed first by the number of framebuffers, then by the number
       of textures each framebuffer has.
      */
      vector< vector< GLuint > >  attachedTextures;
      
      /**
       This is indexed by the number of framebuffers.
       This facilitates the use of glDrawBuffers by storing:
       bufferAttachments[i][numRT] = GL_COLOR_ATTACHMENT0_EXT + numRT;
      */
      vector< GLenum* > bufferAttachments;
      
      /**
       The dimensions of all the textures.
      */
      int w, h;
      
      /**
       There are only this many GL_COLOR_ATTACHMENT0_EXT available currently.
      */
      const static int maxBuffers = 16;
  };
#endif
#endif
  
  