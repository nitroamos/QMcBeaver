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

static bool PRINT = false;
static const char * textureMode = "rgba=32f texRECT rtt";

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
	Matrix(){ rows = 0; columns = 0; textureData = 0;}

	/**the number of rows and columns must be specified in this constructor because there is not a bijection between
	data dimensions and texture dimensions*/
	Matrix(int _row, int _col, RenderTexture * data){
		columns = _col;
		rows = _row;
		textureData = data;
	}

	/**i think we're going to want to look into ways of highly optimizing this process, it will probably be the most
	used constructor*/
	Matrix(Array2D<GLfloat> data){
		columns = data.dim2();
		rows = data.dim1();
		int tw = columns/2.0;
		int th = rows/2.0;
		if(columns%2 != 0) tw += 1;
		if(rows%2 != 0) th += 1;

		textureData = new RenderTexture(textureMode);
		textureData->Initialize(tw,th,true,false);	
		loadMatrix(data);
	}

	Matrix(int _row, int _col, float param1, float param2){
		columns = _col;
		rows = _row;
		int tw = columns/2.0;
		int th = rows/2.0;
		if(columns%2 != 0) tw += 1;
		if(rows%2 != 0) th += 1;

		textureData = new RenderTexture(textureMode);
		textureData->Initialize(tw,th,true,false);
		
		data = Array2D<GLfloat>(rows,columns);
		makeCheckerMatrix(data, param1, param2);
		//makeRandMatrix(data);
		loadMatrix(data);
	}

	Matrix(int _row, int _col){
		columns = _col;
		rows = _row;
		int tw = columns/2.0;
		int th = rows/2.0;
		if(columns%2 != 0) tw += 1;
		if(rows%2 != 0) th += 1;

		textureData = new RenderTexture(textureMode);
		textureData->Initialize(tw,th,true,false);
		
		data = Array2D<GLfloat>(rows,columns);
		makeRandMatrix(data);
		loadMatrix(data);
	}

	/**the initial value is the value meant to initialize the texture. if diagonal matricies are a major
	issue, then we can write a cg script for them too.*/
	Matrix(int _row, int _col, float initial, bool diag){
		columns = _col;
		rows = _row;
		int tw = columns/2.0;
		int th = rows/2.0;
		if(columns%2 != 0) tw += 1;
		if(rows%2 != 0) th += 1;

		textureData = new RenderTexture(textureMode);
		textureData->Initialize(tw,th,true,false);
		
		if(diag){
			data = Array2D<GLfloat>(rows,columns);
			GLfloat** mat = data.array();
			for(int i=0; i<rows; i++)
				for(int j=0; j<columns; j++){
					if(i==j) mat[i][j] = initial;
					else mat[i][j] = 0.0f;
				}
			loadMatrix(data);
		} else {
			operator=(initial);
		}
	}

	Matrix(int dim, float param1, int param2){
		rows = dim;
		columns = dim;
		int td = dim/2.0;
		if(dim%2 != 0) td += 1;

		textureData = new RenderTexture(textureMode);
		textureData->Initialize(td,td,true,false);

		data = Array2D<GLfloat>(rows,columns);
		//makeCheckerMatrix(data, param1, param2);
		makeRandMatrix(data);
		loadMatrix(data);
	}

	/**we really need a deconstructor. but i haven't quite figured out how to deconstruct a RenderTexture yet*/
	void release(){
		data.deallocate();
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
	void loadMatrix(Array2D<GLfloat> &matrix){
		int h = textureData->GetHeight();
		int w = textureData->GetWidth();
		if(PRINT) cout << "data: " << rows << "x" << columns << ", texture: " << h << "x" << w << " dimensional\n";
        GLfloat* texels = new GLfloat[h*w*4];
		GLfloat** mat = matrix.array();

		//trash is the arbitrary value assigned to pixel channels not used by actual data
		//changing this shouldn't affect the validity of the actual answer
		GLfloat trash = 0.5;

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
				texels[index   ] = mat[ i*2 ][ j*2 ];
                texels[index +1] = mat[ i*2 ][j*2+1];
                texels[index +2] = mat[i*2+1][ j*2 ];
                texels[index +3] = mat[i*2+1][j*2+1];
            }
        }

		//the conditional statements here are for the situations where either data dimension is not even.
		//if not even, then the pixels need to be packed in special ways.
		j = w-1;
		if(columns%2 != 0){
			for(int i=0; i<h; i++){
				index = mapping(i, j, h, w);
				texels[index   ] = mat[ i*2 ][ j*2 ];
                texels[index +1] = trash;
				if(rows%2 == 0 || i<h-1)
					texels[index +2] = mat[i*2+1][ j*2 ];
				else
					texels[index +2] = trash;
				texels[index +3] = trash;
			}
		} else {
			for(int i=0; i<h; i++){
				index = mapping(i, j, h, w);
				texels[index   ] = mat[ i*2 ][ j*2 ];
                texels[index +1] = mat[ i*2 ][j*2+1];
				if(rows%2 == 0 || i<h-1){
					texels[index +2] = mat[i*2+1][ j*2 ];
					texels[index +3] = mat[i*2+1][j*2+1];
				} else {
					texels[index +2] = trash;
					texels[index +3] = trash;
				}
			}
		}

		i = h-1;
		if(rows%2 != 0){
			for(int j=0; j<w; j++){
				index = mapping(i, j, h, w);
				texels[index   ] = mat[ i*2 ][ j*2 ];
				if(columns%2 == 0 || j<w-1)
					texels[index +1] = mat[ i*2 ][j*2+1];
				else
					texels[index +2] = trash;
				texels[index +2] = trash;
                texels[index +3] = trash;
			}
		} else {
			for(int j=0; j<w; j++){
				index = mapping(i, j, h, w);
				texels[index   ] = mat[ i*2 ][ j*2 ];
                texels[index +2] = mat[i*2+1][ j*2 ];                
				if(columns%2 == 0 || j<w-1){
					texels[index +1] = mat[ i*2 ][j*2+1];
					texels[index +3] = mat[i*2+1][j*2+1];
				} else {
					texels[index +1] = trash;
					texels[index +3] = trash;
				}
			}
		}

		textureData->BeginCapture();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawPixels(w, h, GL_RGBA, GL_FLOAT, texels);
		textureData->EndCapture();

        getError("Error loading matrix");
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
    RenderTexture * runCg(	float RHS, int idRHS, bool binaryOp, CGprogram fragProg, 
								CGparameter param1, CGparameter param2){
		RenderTexture * result = new RenderTexture(textureMode);
		result->Initialize(textureData->GetWidth(),textureData->GetHeight(),true,false);
		
		//how many hours did it take me to discover that this line needs to be BEFORE the cg stuff??? :-(
		//someone mentioned to me on the opengl forum that with list sharing, it shouldn't matter where
		//this line is... but either i don't understand list sharing (likely) or they were wrong
		//but BeginCapture captures more than the output data, it must "capture" the state information as well
		result->BeginCapture();

		cgGLBindProgram(fragProg);
		cgGLEnableProfile(g_cgProfile);
		
		//specTest();

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
			
		result->EndCapture();
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
		return result;
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

	/*this function is supposed to take the texture out of the GPU's memory and make it accessible
	to the rest of the c++ code by converting the texture to an Array2D object. the input print
	can allow the data to be additionally displayed in the console*/
	Array2D<GLfloat> unloadMatrix(bool print){
		int h = textureData->GetHeight();
        int w = textureData->GetWidth();
		int index;
		GLfloat* result = (GLfloat*)malloc( sizeof(GLfloat) * h * w * 4);
		Array2D<GLfloat> matrix(rows, columns);
		GLfloat** mat = matrix.array();
		
		glFinish();

		textureData->Bind();
		glGetTexImage(textureData->GetTextureTarget(),0,GL_RGBA,GL_FLOAT,result);
		if(print){
			//PrintRGBAPixelsColumn(result,w,h);
			PrintRGBAPixelsBox(result,w,h);
		}

		//  r   g    or    x   y
        //  b   a          z   w
        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++){
				index = mapping(i, j, h, w);
                mat[ i*2 ][ j*2 ] = result[index   ];
				if(j*2+1 < columns)
                mat[ i*2 ][j*2+1] = result[index +1];
                if(i*2+1 < rows)
				mat[i*2+1][ j*2 ] = result[index +2];
				if(i*2+1 < rows && j*2+1 < columns)
                mat[i*2+1][j*2+1] = result[index +3];
            }
        }
		getError("Error unloading matrix");
		return matrix;
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
		RenderTexture * temp = runCg(0,RHS.textureData->GetTextureID(),true,fp,tLHS,tRHS);
		delete textureData;
		textureData = temp;
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
		RenderTexture * temp = runCg(C,-1,true,_fragmentProgram,_textureParam,_valueParam);
		delete textureData;
		textureData = temp;
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

		RenderTexture * temp = runCg(C,-1,true,_fragmentProgram,_textureParam,_valueParam);
		return Matrix(rows,columns,temp);
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
	*/
	Matrix operator*(const Matrix RHS){
		if(columns != RHS.rows){
			cout << "inncorrect dimensions for matrix multiply" << endl;
			return Matrix();
		}
		CGprogram     fp;
		CGparameter   accum, tLHS, tRHS, startOps, stopOps;
		int resultR = textureData->GetHeight(), resultC = RHS.textureData->GetWidth();
		int maxLoops;
		glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &maxLoops);
		int numLoops = textureData->GetWidth();
		int numPasses = ceil((double)numLoops/maxLoops);
		if(PRINT) cout << "numPasses is " << numPasses << endl;

		fp = cgCreateProgram(g_cgContext, CG_SOURCE,
                           matrixMultiply, g_cgProfile,
                           "main", NULL);
        if(fp != NULL){
            cgGLLoadProgram(fp);
			accum = cgGetNamedParameter(fp, "accum");     
            tLHS = cgGetNamedParameter(fp, "texLHS");
			tRHS = cgGetNamedParameter(fp, "texRHS");
			startOps = cgGetNamedParameter(fp, "startOps");
			stopOps = cgGetNamedParameter(fp, "stopOps");
		} else {
			cout << "error in matrixMultiply script" << endl;
			return *this;
		}


		Matrix initial = Matrix(rows,RHS.columns,0.0,false);
		RenderTexture * result[2] = {new RenderTexture(textureMode),initial.textureData};
		result[0]->Initialize(resultC,resultR,true,false);
		bool writingTexture = true;
		//result[writingTexture]->BeginCapture();

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		int maxs = result[writingTexture]->GetMaxS();
		int maxt = result[writingTexture]->GetMaxT();
		for(int i=0; i<numPasses; i++){
			writingTexture = !writingTexture;
			//result[writingTexture]->BeginCapture(result[!writingTexture]);
			result[writingTexture]->BeginCapture();

			//i really think there should be a way of NOT binding the shader and NOT assigning the uniforms for every pass
			//but it works reasonably and i can't get the other, more intuitive way to work
			cgGLBindProgram(fp);
			cgGLEnableProfile(g_cgProfile);
			cgGLSetTextureParameter(tLHS, textureData->GetTextureID());
			cgGLEnableTextureParameter(tLHS);
			cgGLSetTextureParameter(tRHS, RHS.textureData->GetTextureID());
			cgGLEnableTextureParameter(tRHS);
			cgGLSetTextureParameter(accum, result[!writingTexture]->GetTextureID());
			cgGLEnableTextureParameter(accum);

			cgGLSetParameter1f(startOps, i*maxLoops);
			cgGLSetParameter1f(stopOps, i<numPasses-1?(i+1)*maxLoops:numLoops);
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			if(PRINT) cout << "pass number " << (i+1) << ": " << (i*maxLoops) << "-" << (i<numPasses-1?(i+1)*maxLoops:resultC) << " complete" << endl;
			result[writingTexture]->EndCapture();
		}
		
		glFinish();//~600 ms for 1000x1000. the card saves the actual calculations until this line.
		glFlush();

		cgGLDisableTextureParameter(accum);
		cgGLDisableTextureParameter(tLHS);
		cgGLDisableTextureParameter(tRHS);
		cgGLDisableProfile(g_cgProfile);
		
		if(PRINT) cout << "multiplication complete\n";
		getError("Error in matrix multiply");
		return Matrix(rows,RHS.columns,result[writingTexture]);
	}

	void operator=(const float C){
		CGprogram     _fragmentProgram;
		CGparameter   _valueParam;

		_fragmentProgram = cgCreateProgram(g_cgContext, CG_SOURCE,
                           assignConst, g_cgProfile,
                           "main", NULL);

        if(_fragmentProgram != NULL){
            cgGLLoadProgram(_fragmentProgram);
			_valueParam = cgGetNamedParameter(_fragmentProgram, "v");
        }
		RenderTexture * temp = runCg(C,-1,false,_fragmentProgram,_valueParam,_valueParam);
		delete textureData;
		textureData = temp;
	}

	/*this function makes some arbitrary data to view. creates a checker type board.
	the param input is to allow easy creation of different matricies*/
	void makeCheckerMatrix(Array2D<GLfloat> &matrix, GLfloat param1, GLint param2){
		GLfloat c;
		GLfloat** mat = matrix.array();
		for (int i = 0; i < rows; i+=2) {
			for (int j = 0; j < columns; j+=2) {
				//alternates every 8th matrix (not pixel) element
				c = (((i&param2)==0)^((j&param2)==0))*param1;				
				mat[ i ][ j ] = c/3.0;
				if(j<columns-1 || columns%2 == 0)
				mat[ i ][j+1] = c/6.0;
				if(i<rows-1 || rows%2 == 0)
				mat[i+1][ j ] = c/2.0;
				if((j<columns-1 && i<rows-1) || (columns%2 == 0 && rows%2 == 0))
				mat[i+1][j+1] = 1.0;				
			}
		}
	} 

	void makeRandMatrix(Array2D<GLfloat> &matrix){
		GLfloat** mat = matrix.array();
        for(int i=0; i<rows; i++)
			for(int j=0; j<columns; j++){
                mat[i][j] = (GLfloat)(rand()/33000.0);
				//mat[i][j] = (GLfloat)( (int)(mat[i][j]*5) );
			}
    }

    void PrintMatrix(Array2D<GLfloat> matrix){
        //if(!PRINT) return;
		GLfloat** mat = matrix.array();
		for(int i=0; i<matrix.dim1(); i++){
            for(int j=0; j<matrix.dim2() && j < 28; j++){
                printf("%7.3g", (float)mat[i][j]);
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

protected:
	Array2D<GLfloat> data;//at the end, eliminate this var so all mem is in texture?
    int rows, columns;
	RenderTexture *textureData;
};

#endif