/** Programmed by Amos Anderson, amosa@caltech*/
#ifndef GPU_MATRIX
#define GPU_MATRIX

#include "RenderTexture.h"
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <windows.h>
#include <gl/glut.h>
#include <cg/cgGL.h>
#include <math.h>

#include "Array2D.h"
#include "cg_programs.h"
#include "Stopwatch.h"

using namespace std;

#define NUM_REPS 100
static bool PRINT = !true;
static const char * textureMode = "rgba=32f doublebuffer texRECT rtt";

/**
justifications: 
1) i have not removed the functions from this matrix.h file into a
matrix.cpp file because it's annoying messing around with prototypes
2) a lot of this code is not highly organized and the API is not official.
this is just my code "playpen" for trying out different ideas.

i put in a lot of constructors that don't necessarily need to be there.
they are there for my personal convenience, as is the world in general.

future modifications:
1) implement more BLAS stuff. for example, special considerations
for vector operations (e.g. vector-matrix multiplies)
*/
class Matrix {
public:
	Matrix(){ rows = 0; columns = 0; textureData = 0; pixelData = 0;}

	/**the number of rows and columns must be specified in this constructor because there is not a bijection between
	data dimensions and texture dimensions*/
	Matrix(int _row, int _col, RenderTexture * data){
		columns = _col;
		rows = _row;
		textureData = data;
		pixelData = (GLfloat*)malloc( sizeof(GLfloat) * textureData->GetHeight() * textureData->GetWidth() * 4);
	}

	/**i think we're going to want to look into ways of highly optimizing this process, it will probably be the most
	used constructor*/
	Matrix(Array2D<GLfloat> & data){
		columns = data.dim2();
		rows = data.dim1();
		setup();
		loadMatrix(data);
	}

	Matrix(int _row, int _col, float param1, float param2){
		columns = _col;
		rows = _row;
		setup();
		Array2D<GLfloat> data = Array2D<GLfloat>(rows,columns);
		makeCheckerMatrix(data, param1, param2);
		//makeRandMatrix(data);
		loadMatrix(data);
	}

	Matrix(int _row, int _col){
		columns = _col;
		rows = _row;
		setup();
		Array2D<GLfloat> data = Array2D<GLfloat>(rows,columns);
		makeRandMatrix(data);
		//Stopwatch sw = Stopwatch();
		//sw.reset(); sw.start();
		loadMatrix(data);
		//sw.stop();
		//cout << "uploading took " << sw.timeMS() << " ms" << endl;
	}

	/**the initial value is the value meant to initialize the texture. if diagonal matricies are a major
	issue, then we can write a cg script for them too.*/
	Matrix(int _row, int _col, float initial, bool diag){
		columns = _col;
		rows = _row;
		setup();
		if(diag){
			Array2D<GLfloat> data = Array2D<GLfloat>(rows,columns);
			for(int i=0; i<rows; i++)
				for(int j=0; j<columns; j++){
					if(i==j) data(i,j) = initial;
					else data(i,j) = 0.0f;
				}
			loadMatrix(data);
		} else {
			//operator=(initial+3);
			Array2D<GLfloat> data = Array2D<GLfloat>(rows,columns);
			for(int i=0; i<rows; i++)
				for(int j=0; j<columns; j++){
					data(i,j) = initial;
				}
			loadMatrix(data);			
		}
	}

	Matrix(int dim, float param1, int param2){
		rows = dim;
		columns = dim;
		setup();
		Array2D<GLfloat> data = Array2D<GLfloat>(rows,columns);
		//makeCheckerMatrix(data, param1, param2);
		makeRandMatrix(data);
		loadMatrix(data);
	}

	/*i was having problems with a destructor. it might be a compiler problem?
	basically, it would run the default constructor for any object created, 
	then run the destructor on that object deleting the same memory. we'll do
	this manually.*/
	void destroy(){
		delete [] pixelData;
		//Reset both deletes the textures and returns buffer memory
		textureData->Reset(0,0);
	}

	void setup(){
		int tw = columns/2.0;
		int th = rows/2.0;
		if(columns%2 != 0) tw += 1;
		if(rows%2 != 0) th += 1;

		textureData = new RenderTexture(textureMode);
		textureData->Initialize(tw,th,true,false);
		textureData->BeginCapture();
		glDrawBuffer(GL_BACK);
		glReadBuffer(GL_BACK);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
		glDrawBuffer(GL_FRONT);
		glReadBuffer(GL_FRONT);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		textureData->BindBuffer(WGL_BACK_LEFT_ARB); 
		textureData->EndCapture();
		pixelData = new GLfloat[th*tw*4];
		//cout << "matrix created, " << pixelData << " size " << (sizeof(GLfloat) * th * tw * 4) << endl;
		//pixelData = (GLfloat*)malloc( sizeof(GLfloat) * th * tw * 4);
	}

	/** there is a question of how to store 2D data in a 1D array, which is required by OpenGl. all other functions
	in this class rely on this function (which should may be inline or #defined... i and j are the coordinates
	in the data matrix while h and w are the dimensions of the texture matrix. i don't think this mapping needs to change
	if we switch to the column-data pixel representation */
	int mapping(int i, int j, int h, int w){
		//matrix multiplication was harder using the commented mapping
		//return 4*((h-i-1)*w + j);

		//the natural mapping
		return 4*(i*w + j);
	}

	/*the purpose of this function is to load data from an Array2D object, and place it in
	texture memory. the big deal in this function is 1st) making sure that all the matrix data gets
	into pixel data, and 2nd) to assign some value to all the pixels. All the if statements attempt
	to take care of the problem when the dimensions of the data are odd. this method probably needs
	to be highly optimized because it will be used frequently.
	my idea for this is to pass the matrix data into vertex shaders, who will pass the info to fragment shaders
	in a one-to-one correspondance if (i think) there are as many verticies as pixels. this means that the 
	fragment shader won't have to do a texture lookup. this will probably save time.
	
	there are probably other C++ ways of optimizing this. valarray? inline mapping? i think that the drawPixels
	function is actually the most expensive part of this function
	
	the majority of this code is in handling fringe effects*/
	void loadMatrix(const Array2D<GLfloat> &matrix){
		int h = textureData->GetHeight();
		int w = textureData->GetWidth();
		if(PRINT) cout << "data: " << rows << "x" << columns << ", texture: " << h << "x" << w << " dimensional\n";

		GLfloat trash = 0.0;

		/*there is a choice in how to represent a matrix in a texture. either in a box:
          r   g
          b   a
        or as a column
		  r
		  g
		  b
		  a
		where each of the letters rgba represent one number in the matrix. either way, four numbers from the matrix
		are stored in each pixel. maybe i'll change this later, but for now, i choose the first representation.*/
		int i, j, index;
		for(i=0; i<h-1; i++){
            for(j=0; j<w-1; j++){
				index = mapping(i, j, h, w);
				pixelData[index   ] = matrix.get( i*2   , j*2   );
                pixelData[index +1] = matrix.get( i*2   , j*2+1 );
                pixelData[index +2] = matrix.get( i*2+1 , j*2   );
                pixelData[index +3] = matrix.get( i*2+1 , j*2+1 );
            }
        }

		//the conditional statements here are for the situations where either data dimension is not even.
		//if not even, then the pixels need to be packed in special ways.
		j = w-1;
		if(columns%2 != 0){
			for(int i=0; i<h; i++){
				index = mapping(i, j, h, w);
				pixelData[index   ] = matrix.get( i*2   , j*2   );
                pixelData[index +1] = trash;
				if(rows%2 == 0 || i<h-1)
					pixelData[index +2] = matrix.get( i*2+1 , j*2   );
				else
					pixelData[index +2] = trash;
				pixelData[index +3] = trash;
			}
		} else {
			for(int i=0; i<h; i++){
				index = mapping(i, j, h, w);
				pixelData[index   ] = matrix.get( i*2   , j*2   );
                pixelData[index +1] = matrix.get( i*2   , j*2+1 );
				if(rows%2 == 0 || i<h-1){
					pixelData[index +2] = matrix.get( i*2+1 , j*2   );
					pixelData[index +3] = matrix.get( i*2+1 , j*2+1 );
				} else {
					pixelData[index +2] = trash;
					pixelData[index +3] = trash;
				}
			}
		}

		i = h-1;
		if(rows%2 != 0){
			for(int j=0; j<w; j++){
				index = mapping(i, j, h, w);
				pixelData[index   ] = matrix.get( i*2   , j*2   );
				if(columns%2 == 0 || j<w-1)
					pixelData[index +1] = matrix.get( i*2   , j*2+1 );
				else
					pixelData[index +2] = trash;
				pixelData[index +2] = trash;
                pixelData[index +3] = trash;
			}
		} else {
			for(int j=0; j<w; j++){
				index = mapping(i, j, h, w);
				pixelData[index   ] = matrix.get( i*2   , j*2   );
                pixelData[index +2] = matrix.get( i*2+1 , j*2   );                
				if(columns%2 == 0 || j<w-1){
					pixelData[index +1] = matrix.get( i*2   , j*2+1 );
					pixelData[index +3] = matrix.get( i*2+1 , j*2+1 );
				} else {
					pixelData[index +1] = trash;
					pixelData[index +3] = trash;
				}
			}
		}
		//*
		//textureData->BindBuffer(WGL_BACK_LEFT_ARB);
		Stopwatch sw = Stopwatch();
		sw.reset(); sw.start();
		for(int j=0; j<NUM_REPS; j++){
		textureData->Bind();
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
					 w, h, 0, GL_RGBA, GL_FLOAT, pixelData);
		}
		glFinish();
		sw.stop();
		//cout << "loading took " << sw.timeMS()/(double)NUM_REPS << endl;
		moveTextureToBuffer();
		/*/
		textureData->BeginCapture();
		glDrawPixels(w, h, GL_RGBA, GL_FLOAT, pixelData);
		textureData->EndCapture();
		//*/
        getError("Error loading matrix");
	}

	/*this function is supposed to take the texture out of the GPU's memory and make it accessible
	to the rest of the c++ code by converting the texture to an Array2D object. the input print
	can allow the data to be additionally displayed in the console
	calling this function will store the results in the data class variable.*/
	void unloadMatrix(bool print){;
		int h = textureData->GetHeight();
        int w = textureData->GetWidth();			
		glFinish();
		
		//this is a *lot* faster than glGetTexImage(textureData->GetTextureTarget(),0,GL_RGBA,GL_FLOAT,result);
		Stopwatch sw = Stopwatch();
		sw.reset(); sw.start();
		for(int j=0; j<NUM_REPS; j++){
			textureData->BeginCapture();
			glReadPixels(0,0,w,h,GL_RGBA,GL_FLOAT,pixelData);
			textureData->EndCapture();
		}
		glFinish();
		sw.stop();
		//cout << "unloading took " << sw.timeMS()/(double)NUM_REPS << endl;

		if(print){
			cout << "unloading data...\n";
			//PrintRGBAPixelsColumn(result,w,h);
			PrintRGBAPixelsBox(pixelData,w,h);
		}

		getError("Error unloading matrix");
	}

	/* engine to run simple cg scripts. this is really for simple math operators who only take at most 2 parameters
	the function can be summarized as either:
	result = LHS op RHS
	or
	result = op RHS (e.g. result = scalar or (not implemented yet) result = RHS texture)
	where result and LHS are textures, and RHS can be either a scalar or a texture
	the matrix multiplication can't go through this because the C++ code needs to be different.

	RHS is (if used) a constant parameter
	idRHS is (if used) the texture id (you can't use both RHS and idRHS)
	binaryOp clues the code whether whether there is an LHS or not
	fragProg is the created program object (not bound yet)
	param1, param2 are either uniform constants or texture objects
	*/
    void runCg(	float RHS, int idRHS, bool binaryOp, CGprogram fragProg, 
								CGparameter param1, CGparameter param2){
		textureData->BeginCapture();
		cgGLBindProgram(fragProg);
		cgGLEnableProfile(g_cgProfile);
		textureData->swapBuffers();
		specTest();

		if(binaryOp){			
			cgGLSetTextureParameter(param1, textureData->GetTextureID());			
			cgGLEnableTextureParameter(param1);
			if(idRHS < 0){
				cgGLSetParameter1f(param2, RHS);
			} else {
				cgGLSetTextureParameter(param2, (unsigned int)idRHS);
				cgGLEnableTextureParameter(param2);
			}
		} else {
			if(idRHS >= 0){
				cgGLSetTextureParameter(param1, (unsigned int)idRHS);
				cgGLEnableTextureParameter(param1);
			} else {
				cgGLSetParameter1f(param1, RHS);
			}
		}		

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		int maxs = textureData->GetMaxS();
		int maxt = textureData->GetMaxT(); 
		glBegin(GL_QUADS);
			glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
			glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
			glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
			glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
		glEnd();
		textureData->EndCapture();
		glFinish();
		glFlush();
		
		if(binaryOp || idRHS >= 0){
			cgGLDisableTextureParameter(param1);
		} 
		if(binaryOp && idRHS >= 0){
			cgGLDisableTextureParameter(param2);
		}
		cgGLDisableProfile(g_cgProfile);
		
		getError("Error in runCg function");
    }    

	/*this function simply displays whatever data is currently in the textureData. It needs
	to be run through a copyTexture program to convert the texture data into color data. This is
	purely an aesthetic part of the code -- eventually, we'll just empty this function
	an idea: for display purposes, we could normalize all the data values to the largest matrix entry.
	this would make something actually viewable on the screen instead of just white*/
	void display(){
		CGprogram     _fragmentProgram;  // the fragment program used to update
		CGparameter   _textureParam;     // a parameter to the fragment program

		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           copyTexture, g_cgProfile,
                           "main", NULL);

        if(_fragmentProgram != NULL){
            cgGLLoadProgram(_fragmentProgram);
            _textureParam = cgGetNamedParameter(_fragmentProgram, "texture");
        }

		cgGLBindProgram(_fragmentProgram);
		cgGLEnableProfile(g_cgProfile);
		cgGLSetTextureParameter(_textureParam, textureData->GetTextureID());
		cgGLEnableTextureParameter(_textureParam);

		glClearColor(0.7f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glFinish();

		int maxs = textureData->GetMaxS();
		int maxt = textureData->GetMaxT(); 
		glBegin(GL_QUADS);
			glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
			glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
			glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
			glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
		glEnd();

		cgGLDisableTextureParameter(_textureParam);
		cgGLDisableProfile(g_cgProfile);
		
		getError("Error in display");
		glFinish();
	}

	void operator+=(const Matrix RHS){
		if(rows != RHS.rows || columns != RHS.columns) return;
		CGprogram     fp;             // the fragment program used to update
		CGparameter   tLHS, tRHS;     // a parameter to the fragment program

		fp = cgCreateProgram(g_cgContext, CG_SOURCE,
                           addTextures, g_cgProfile,
                           "main", NULL);

        // Create the texture parameter for the fragment program
        if(fp != NULL){
            cgGLLoadProgram(fp);
            tLHS = cgGetNamedParameter(fp, "texLHS");
			tRHS = cgGetNamedParameter(fp, "texRHS");			
        }
		runCg(0,RHS.textureData->GetTextureID(),true,fp,tLHS,tRHS);
	}

	void operator+=(const float C){
		CGprogram     _fragmentProgram;
		CGparameter   _textureParam;
		CGparameter   _valueParam;

		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           addConstant, g_cgProfile,
                           "main", NULL);

        if(_fragmentProgram != NULL){
            cgGLLoadProgram(_fragmentProgram);
            _textureParam = cgGetNamedParameter(_fragmentProgram, "texture");
			_valueParam = cgGetNamedParameter(_fragmentProgram, "value");
        }
		runCg(C,-1,true,_fragmentProgram,_textureParam,_valueParam);
		cleanFringes();
	}

	void operator*=(const float C){
		Matrix m = (*this)*C;
		this->textureData = m.textureData;
	}

	Matrix operator*(const float C){
		CGprogram     _fragmentProgram;
		CGparameter   _textureParam;
		CGparameter   _valueParam;

		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           multiplyConstant, g_cgProfile,
                           "main", NULL);

		if(_fragmentProgram != NULL){
            cgGLLoadProgram(_fragmentProgram);
            _textureParam = cgGetNamedParameter(_fragmentProgram, "texture");
			_valueParam = cgGetNamedParameter(_fragmentProgram, "value");
        }

		runCg(C,-1,true,_fragmentProgram,_textureParam,_valueParam);
		cout << "this method is broken. runCg affects this textureData, doesn't create a new one\n";
		return Matrix(0,0);
	}

	/*most of the variables are self explanatory. the startOps and stopOps direct the shader to which elements of the 
	inner product it should be accumulating because if the textures are too large, it can't do the complete summation
	in one pass due to limitations in either shader instruction count limitations or for loop limitations
	
	no attempt is made to deal with the "fringe effects" of a possible extra row or column of data. based on how
	the data is loaded, it doesn't change the answer.

	we need to look into ways of speeding this up because it really is quite slow compared to what other
	GPU algorithms are able to get. some ideas:
	1) shader in assembly
	2) instead of using uniform parameters, hardcode the values into the shader program
	
	LHS = this*RHS
	*/
	void matrixMultiplyFaster(Matrix &RHS, Matrix &LHS){
		if(columns != RHS.rows){
			cout << "inncorrect dimensions for matrix multiply" << endl;
			return;
		}
		CGprogram     fp;
		CGparameter   accum, tLHS, tRHS, startOps, stopOps;
		const int maxLoops = 255;
		//glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &maxLoops);
		const int numLoops = textureData->GetWidth();
		const int numPasses = ceil((double)numLoops/maxLoops);
		//const int numPasses = 2;
		//const char * arg[] = {"fastprecision"};
		fp = cgCreateProgram(g_cgContext, CG_SOURCE,matrixMultiply, g_cgProfile,"main", NULL);
        //fp = cgCreateProgram(g_cgContext, CG_SOURCE,testInputs, g_cgProfile,"main", NULL);

		if(fp != NULL){
            cgGLLoadProgram(fp);
			accum = cgGetNamedParameter(fp, "accum");     
            tLHS = cgGetNamedParameter(fp, "texLHS");
			tRHS = cgGetNamedParameter(fp, "texRHS");
			startOps = cgGetNamedParameter(fp, "startOps");
			stopOps = cgGetNamedParameter(fp, "stopOps");
		} else {
			cout << "error in matrixMultiply script" << endl;
			return;
		}

		RenderTexture * result = LHS.textureData;
		result->BeginCapture();
		cgGLBindProgram(fp);
		cgGLEnableProfile(g_cgProfile);
		textureData->BindBuffer(WGL_FRONT_LEFT_ARB);
		cgGLSetTextureParameter(tLHS, textureData->GetTextureID());
		cgGLEnableTextureParameter(tLHS);
		RHS.textureData->BindBuffer(WGL_FRONT_LEFT_ARB);
		cgGLSetTextureParameter(tRHS, RHS.textureData->GetTextureID());
		cgGLEnableTextureParameter(tRHS);
		cgGLSetTextureParameter(accum, result->GetTextureID());
		cgGLEnableTextureParameter(accum);

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		if(numPasses%2 == 1) result->swapBuffers();
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		result->EndCapture();		
		Stopwatch sw = Stopwatch();
		sw.reset(); sw.start();
		int maxs = result->GetMaxS();
		int maxt = result->GetMaxT();
		for(int i=0; i<numPasses; i++){
			result->BeginCapture();
			result->swapBuffers();
			
			cgGLSetParameter1f(startOps, i*maxLoops);
			cgGLSetParameter1f(stopOps, i<numPasses-1?(i+1)*maxLoops:numLoops);
			
			for(int j=0; j<NUM_REPS; j++){
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			}
			result->EndCapture();			
		}
		glFinish();
		glFlush();
		sw.stop();
		//cout << "multiplying took " << sw.timeMS()/(double)NUM_REPS << endl;

		
		cgGLDisableTextureParameter(accum);
		cgGLDisableTextureParameter(tLHS);
		cgGLDisableTextureParameter(tRHS);
		cgGLDisableProfile(g_cgProfile);
		//LHS.unloadMatrix(true);
		//unloadMatrix(true);
		//RHS.unloadMatrix(true);
		getError("Error in matrix multiply");
	}

		//LHS = this*RHS
	void matrixMultiplyDoNothing(Matrix &RHS, Matrix &LHS){
		if(columns != RHS.rows){
			cout << "inncorrect dimensions for matrix multiply" << endl;
			return;
		}
		CGprogram     fp;
		CGparameter   accum, tLHS, tRHS, startOps, stopOps;
		const int maxLoops = 255;
		//glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &maxLoops);
		const int numLoops = textureData->GetWidth();
		const int numPasses = ceil((double)numLoops/maxLoops);
		//const char * arg[] = {"fastprecision"};
        fp = cgCreateProgram(g_cgContext, CG_SOURCE,testInputs, g_cgProfile,"main", NULL);

		if(fp != NULL){
            cgGLLoadProgram(fp);
			accum = cgGetNamedParameter(fp, "accum");     
            tLHS = cgGetNamedParameter(fp, "texLHS");
			tRHS = cgGetNamedParameter(fp, "texRHS");
			startOps = cgGetNamedParameter(fp, "startOps");
			stopOps = cgGetNamedParameter(fp, "stopOps");
		} else {
			cout << "error in matrixMultiply script" << endl;
			return;
		}

		RenderTexture * result = LHS.textureData;
		result->BeginCapture();
		cgGLBindProgram(fp);
		cgGLEnableProfile(g_cgProfile);
		textureData->BindBuffer(WGL_FRONT_LEFT_ARB);
		cgGLSetTextureParameter(tLHS, textureData->GetTextureID());
		cgGLEnableTextureParameter(tLHS);
		RHS.textureData->BindBuffer(WGL_FRONT_LEFT_ARB);
		cgGLSetTextureParameter(tRHS, RHS.textureData->GetTextureID());
		cgGLEnableTextureParameter(tRHS);
		cgGLSetTextureParameter(accum, result->GetTextureID());
		cgGLEnableTextureParameter(accum);

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		if(numPasses%2 == 1) result->swapBuffers();
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		result->EndCapture();		

		int maxs = result->GetMaxS();
		int maxt = result->GetMaxT();
		for(int i=0; i<numPasses; i++){
			result->BeginCapture();
			result->swapBuffers();
			
			cgGLSetParameter1f(startOps, i*maxLoops);
			cgGLSetParameter1f(stopOps, i<numPasses-1?(i+1)*maxLoops:numLoops);
			for(int j=0; j<NUM_REPS; j++){
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			}
			
			result->EndCapture();			
		}

		glFinish();
		glFlush();
		cgGLDisableTextureParameter(accum);
		cgGLDisableTextureParameter(tLHS);
		cgGLDisableTextureParameter(tRHS);
		cgGLDisableProfile(g_cgProfile);
		//LHS.unloadMatrix(true);
		//unloadMatrix(true);
		//RHS.unloadMatrix(true);
		getError("Error in matrix multiply");
	}

	void operator=(const float C){
		CGprogram     _fragmentProgram;
		CGparameter   _valueParam;

		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           assignConst, g_cgProfile,
                           "main", NULL);
		/*
		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           assignConstFixRight, g_cgProfile,
                           "main", NULL);
		*/
        if(_fragmentProgram != NULL){
            cgGLLoadProgram(_fragmentProgram);
			_valueParam = cgGetNamedParameter(_fragmentProgram, "v");
        }
		runCg(C,-1,false,_fragmentProgram,_valueParam,_valueParam);
	}

	/*this function makes some arbitrary data to view. creates a checker type board.
	the param input is to allow easy creation of different matricies*/
	void makeCheckerMatrix(Array2D<GLfloat> &matrix, GLfloat param1, GLint param2){
		GLfloat c;
		for (int i = 0; i < rows; i+=2) {
			for (int j = 0; j < columns; j+=2) {
				//alternates every 8th matrix (not pixel) element
				c = (((i&param2)==0)^((j&param2)==0))*param1;				
				matrix(i,j) = c/3.0;
				if(j<columns-1 || columns%2 == 0)
				matrix(i,j+1) = c/6.0;
				if(i<rows-1 || rows%2 == 0)
				matrix(i+1,j) = c/2.0;
				if((j<columns-1 && i<rows-1) || (columns%2 == 0 && rows%2 == 0))
				matrix(i+1,j+1) = 1.0;				
			}
		}
	} 

	void makeRandMatrix(Array2D<GLfloat> &matrix){
        for(int i=0; i<matrix.dim1(); i++)
			for(int j=0; j<matrix.dim2(); j++){
				matrix(i,j) = (GLfloat)(rand()/33000.0);
				//mat[i][j] = (GLfloat)( (int)(mat[i][j]*5) );
			}
    }

    void PrintMatrix(Array2D<GLfloat> matrix){
        //if(!PRINT) return;
		for(int i=0; i<matrix.dim1(); i++){
            for(int j=0; j<matrix.dim2() && j < 28; j++){
                printf("%7.3g", (float)matrix(i,j));
            }
            printf("\n");
        }
        printf("\n");
    }

    void PrintArray(float* pix, int dim){
        if(!PRINT) return;
        for(int i=0; i<dim; i++)
            printf("%10.10g",pix[i]);
        cout << endl;
    }

	/*provides a text representation of the pixels:
		r	g
		b	a  */
	void PrintRGBAPixelsBox(float* pix, int w, int h){
		int index;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w && j < 12; j++) {
				index = mapping(i, j, h, w);
				printf("%6.3f",pix[index +0]);
				printf("%6.3f",pix[index +1]);
				printf("   ");
			}
			printf("\n");
			for (int j = 0; j < w && j < 12; j++) {
				index = mapping(i, j, h, w);
				printf("%6.3f",pix[index +2]);
				printf("%6.3f",pix[index +3]);
				printf("   ");
			}
			printf("\n\n");  
		}
		printf("\n");
	}

	/*This function will print pixels to a terminal screen in a column format. this representation
	isn't as wide in a terminal screen*/
	void PrintRGBAPixelsColumn(float* pix, int w, int h){
		//see loadMatrix for disc of matrix-pixel mapping
		for (int i = 0; i < 4*h; i++) {
			for (int j = 0; j < w && j < 26; j++) {//> around 26 will wrap 
				printf("%7.3f",pix[4*((h-(i/4)-1)*w + j)+i%4]);
				//the %12.8f will show pretty much the max precision 
				//printf("%12.8f",pix[4*((h-(i/4)-1)*w + j)+i%4]);
			}
			printf("\n");
			if(i%4==3) printf("\n");
		}
		printf("\n");
	}

	int getRows(){ return rows; }
	int getColumns(){ return columns; }
	bool isNull(){ return (rows == 0 || columns == 0); }
	
	GLfloat& operator()(int i,int j){
		return pixelData[mapping(i,j,textureData->GetHeight(),textureData->GetWidth())];
	}

	void moveTextureToBuffer(){
		CGprogram     _fragmentProgram;
		CGparameter   _textureParam;

		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           copyTexture, g_cgProfile,
                           "main", NULL);

        if(_fragmentProgram != NULL){
            cgGLLoadProgram(_fragmentProgram);
            _textureParam = cgGetNamedParameter(_fragmentProgram, "texture");
        }

		textureData->BeginCapture();
		cgGLBindProgram(_fragmentProgram);
		cgGLEnableProfile(g_cgProfile);
		cgGLSetTextureParameter(_textureParam, textureData->GetTextureID());
		cgGLEnableTextureParameter(_textureParam);
		
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		int maxs = textureData->GetMaxS();
		int maxt = textureData->GetMaxT(); 
		glBegin(GL_QUADS);
			glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
			glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
			glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
			glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
		glEnd();
		glFinish();
		textureData->EndCapture();

		cgGLDisableTextureParameter(_textureParam);
		cgGLDisableProfile(g_cgProfile);		
		getError("Error in display");
	}

	void cleanFringes(){
		CGprogram     _fragmentProgram;
		CGparameter   _textureParam, _valueParam;

		if(rows%2 == 1){
			_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
							fixRightFringe, g_cgProfile,
							"main", NULL);

			if(_fragmentProgram != NULL){
				cgGLLoadProgram(_fragmentProgram);
				_textureParam = cgGetNamedParameter(_fragmentProgram, "texture");
				_valueParam = cgGetNamedParameter(_fragmentProgram, "value");
			}

			textureData->BeginCapture();
			textureData->swapBuffers();
			cgGLBindProgram(_fragmentProgram);
			cgGLEnableProfile(g_cgProfile);
			cgGLSetTextureParameter(_textureParam, textureData->GetTextureID());
			cgGLEnableTextureParameter(_textureParam);

			glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			int maxs = textureData->GetMaxS();
			int maxt = textureData->GetMaxT();
			cgGLSetParameter1f(_valueParam,maxs-0.51);
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			glFinish();
			textureData->EndCapture();

			cgGLDisableTextureParameter(_textureParam);
			cgGLDisableProfile(g_cgProfile);
		}

		if(columns%2 == 1){
			_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
							fixBottomFringe, g_cgProfile,
							"main", NULL);

			if(_fragmentProgram != NULL){
				cgGLLoadProgram(_fragmentProgram);
				_textureParam = cgGetNamedParameter(_fragmentProgram, "texture");
				_valueParam = cgGetNamedParameter(_fragmentProgram, "value");
			}

			textureData->BeginCapture();
			textureData->swapBuffers();
			cgGLBindProgram(_fragmentProgram);
			cgGLEnableProfile(g_cgProfile);
			cgGLSetTextureParameter(_textureParam, textureData->GetTextureID());
			cgGLEnableTextureParameter(_textureParam);

			glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			int maxs = textureData->GetMaxS();
			int maxt = textureData->GetMaxT();
			cgGLSetParameter1f(_valueParam,maxs-0.51);
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			glFinish();
			textureData->EndCapture();

			cgGLDisableTextureParameter(_textureParam);
			cgGLDisableProfile(g_cgProfile);
		}

		getError("Error in display");
	}
	/*to do: implement method of testing whether unloadMatrix actually needs to be called*/
	Array2D<GLfloat> getData(){
		unloadMatrix(false);
		int h = textureData->GetHeight();
        int w = textureData->GetWidth();
		int index;
		Array2D<GLfloat> data = Array2D<GLfloat>(rows,columns);

		//  r   g    or    x   y
        //  b   a          z   w
		//this section could be optimized
        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++){
				index = mapping(i, j, h, w);
                data(i*2,j*2) = pixelData[index   ];
				if(j*2+1 < columns)
                data(i*2,j*2+1) = pixelData[index +1];
                if(i*2+1 < rows)
				data(i*2+1,j*2) = pixelData[index +2];
				if(i*2+1 < rows && j*2+1 < columns)
                data(i*2+1,j*2+1) = pixelData[index +3];
            }
        }
		//cout << "recieved...\n";
		//cout << endl << data << endl;
		return data;
	}

	void setTexture(int _rows, int _columns, RenderTexture * newTex){
		rows = _rows; columns = _columns;
		textureData = newTex;
	}

protected:
	//implement this!
	//bool sync;

	//there is no bijection between these numbers and the texture dimensions
	int rows, columns;
	
	//eventually to be replaced with superbuffers!!
	RenderTexture *textureData;

	//making a buffer of GLfloat*'s available aleviates a huge memory leak problem
	GLfloat* pixelData;
};

#endif