#ifndef GPUMATRIX
#define GPUMATRIX

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include "GPUGlobals.h"
#include "RenderTexture.h"
#include <math.h>

#include "Array1D.h"
#include "Array2D.h"
#include "Stopwatch.h"

/*
deltaOE = 38
deltaBF = 258
1/deltaOE = 0.026315789
*/

static const char * qmcMatrixProgram =
"float4 main(in float2 position : TEX0,                              \n"
"            uniform int  deltaOE,                                   \n"
"            uniform int  deltaBF,                                   \n"
"            uniform int  startOps,                                  \n"
"            uniform int  stopOps,                                   \n"
"            uniform samplerRECT  accum,                             \n"
"            uniform samplerRECT  elbf,                              \n"
"            uniform samplerRECT  oebf) : COLOR                      \n"
"{                                                                   \n"
"    int c = position.x/deltaOE;                                     \n"
"    int coeffRow = fmod(position.x,deltaOE);                        \n"
"    int shift = c*deltaBF;                                          \n"
"    float4 output = texRECT(accum, position);                       \n"
"    float4 a, b;                                                    \n"
"    for (int i = startOps; i < stopOps; i++) {                      \n"
"       a  = texRECT(elbf, float2(i + shift, position.y));           \n"
"       b  = float4(1,2,3,4);                                        \n"
"       b  = texRECT(oebf, float2(i, coeffRow));                     \n"
"       output += a.xxzz*b.xzxz;                                     \n"
"       output += a.yyww*b.ywyw;                                     \n"
"    }                                                               \n"
"    return output;                                                  \n"
"}                                                                   \n";

/**The Coeff and basisFunctions submatricies have the same shape because the
Coeff matrix is passed in transposed.*/
class GPUQMCMatrix {
public:
	GPUQMCMatrix(){
		nOrbitals = 0;
		nBasisF = 0;
		nCalcs = 0;
	}

	GPUQMCMatrix(Array2D<QMCfloat> & Coeffs, int numCalcs){
		nCalcs = numCalcs;
		nOrbitals = Coeffs.dim1();
		nBasisF = Coeffs.dim2();

		deltaBF = nBasisF/2.0;
		deltaOE = nOrbitals/2.0;
		if(nBasisF%2 != 0)   deltaBF += 1;
		if(nOrbitals%2 != 0) deltaOE += 1;

		gpuData = new RenderTexture("rgba=32f doublebuffer texRECT rtt");
		gpuData->Initialize(nCalcs*deltaOE,5*deltaOE,true,false);
		cpuData = new GLfloat[ nCalcs*deltaOE * 5*deltaOE * 4 ];
		//pixelData = new GLfloat[ 5*deltaOE * 2*deltaBF * 4 ];
		pixelData = (GLfloat *) calloc( 5*deltaOE * 2*deltaBF * 4 , sizeof(GLfloat) );
		glGenTextures(2, operands);
		fp = cgCreateProgram(g_cgContext, CG_SOURCE,qmcMatrixProgram, g_cgProfile,"main", NULL);

		if(fp != NULL){
            cgGLLoadProgram(fp);
			tELxOR     = cgGetNamedParameter(fp, "accum");     
            tELxBF     = cgGetNamedParameter(fp, "elbf");
			tORxBF     = cgGetNamedParameter(fp, "oebf");
			startOps   = cgGetNamedParameter(fp, "startOps");
			stopOps    = cgGetNamedParameter(fp, "stopOps");
			cg_deltaOE = cgGetNamedParameter(fp, "deltaOE");
			cg_deltaBF = cgGetNamedParameter(fp, "deltaBF");
		} else {
			cout << "error in matrixMultiply script" << endl;
			return;
		}

		numLoops = deltaBF;

		//*
		//maxLoops = 255;
		maxLoops = 50;
		/*/
		glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &maxLoops);
		//*/
		
		//*
		numPasses = ceil((double)numLoops/maxLoops);
		/*/
		numPasses = 2;
		//*/

		loadCoeffs(Coeffs);
	}

	void destroy(){
		delete [] cpuData;
		gpuData->Reset(0,0);
	}

	void runCalculation(ArrayGPU & basisFunctions, ArrayGPU & results){
		if(nCalcs == 0) cout << "ERROR: initialize GPUQMCMatrix first!\n";
		if(basisFunctions.dim1() != nCalcs*5) cout << "ERROR: wrong num calcs!\n";
		Stopwatch sw = Stopwatch();
		if(showTimings){ sw.reset(); sw.start(); }
		loadBasisF(basisFunctions);
		if(showTimings){ sw.stop(); cout << " loading: " << sw.timeMS(); }
		cleanRenderTexture();

		if(showTimings){ sw.reset(); sw.start(); }

		gpuData->BeginCapture();
		cgGLBindProgram(fp);
		cgGLEnableProfile(g_cgProfile);

		if(numPasses%2 == 1) gpuData->swapBuffers();

		cgGLSetParameter1f(cg_deltaOE, deltaOE);
		cgGLSetParameter1f(cg_deltaBF, deltaBF);
		cgGLSetTextureParameter(tELxBF, operands[basis]);
		cgGLSetTextureParameter(tORxBF, operands[coeff]);
		cgGLSetTextureParameter(tELxOR, gpuData->GetTextureID());
		cgGLEnableTextureParameter(tELxBF);
		cgGLEnableTextureParameter(tORxBF);
		cgGLEnableTextureParameter(tELxOR);
		gpuData->EndCapture();		

		int maxs = gpuData->GetMaxS();
		int maxt = gpuData->GetMaxT();
		for(int i=0; i<numPasses; i++){
			gpuData->BeginCapture();
			gpuData->swapBuffers();
			
			cgGLSetParameter1f(startOps, i*maxLoops);
			cgGLSetParameter1f(stopOps, i<numPasses-1?(i+1)*maxLoops:numLoops);
			
			//for(int i=0; i<50; i++){
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			//}
			gpuData->EndCapture();			
		}
		glFinish();
		glFlush();
		
		cgGLDisableTextureParameter(tELxBF);
		cgGLDisableTextureParameter(tORxBF);
		cgGLDisableTextureParameter(tELxOR);
		cgGLDisableProfile(g_cgProfile);

		if(showTimings){ sw.stop(); cout << " cg: " << sw.timeMS(); }

		if(showTimings){ sw.reset(); sw.start(); }
		unloadData(results);
		if(showTimings){ sw.stop(); cout << " unloading: " << sw.timeMS(); }
		
		getOpenGLError("Error in QMC matrix multiply");
	}

private:
	void cleanRenderTexture(){
		gpuData->BeginCapture();
		glDrawBuffer(GL_BACK);
		glReadBuffer(GL_BACK);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
		glDrawBuffer(GL_FRONT);
		glReadBuffer(GL_FRONT);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		gpuData->BindBuffer(WGL_BACK_LEFT_ARB); 
		gpuData->EndCapture();
	}

	void loadCoeffs(Array2D<QMCfloat> & Coeffs){
		ArrayGPU coeffA = ArrayGPU(1);
		coeffA(0) = Coeffs;
		glActiveTextureARB(GL_TEXTURE1_ARB);
		glBindTexture(GL_TEXTURE_RECTANGLE_NV, operands[coeff]);
		mapData(coeffA,true);
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
					 deltaBF, deltaOE, 0, GL_RGBA, GL_FLOAT, 
					 pixelData);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        getOpenGLError("Error loading Coeffs");
	}

	void loadBasisF(ArrayGPU & basisFunctions){
		glActiveTextureARB(GL_TEXTURE0_ARB);
		glBindTexture(GL_TEXTURE_RECTANGLE_NV, operands[basis]);
		mapData(basisFunctions,false);
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
					 nCalcs*deltaBF, 5*deltaOE, 0, GL_RGBA, GL_FLOAT, 
					 pixelData);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        getOpenGLError("Error loading basis function data");
	}

	/**r and c indicate which submatrix we're in, i and j indicate which value in that submatrix*/
	inline int operandMapping(int r, int c, int i, int j, int numCols){
		//4*( (row selection)*(total width) + (column selection) );
		return 4*( (r*deltaOE + i)*numCols*deltaBF + (c*deltaBF + j) );
	}

	void mapData(ArrayGPU &data, bool isCoeff){
		Stopwatch sw = Stopwatch();
		if(showTimings){ sw.reset(); sw.start(); }
		int numRows = isCoeff? 1 : 5;
		int numCols = isCoeff? 1 : nCalcs;
		int h = deltaOE;
		int w = deltaBF;

		GLfloat fringe = 0.0;
		Array2D<QMCfloat> * matrix;

		//we iterate over all the matricies		
		for(int c = 0; c < numCols; c++){
			for(int r = 0; r < numRows; r++){
				matrix = & data( c*numRows + r);
				
				int i, j, index;
				for(i=0; i<h-1; i++){
					for(j=0; j<w-1; j++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix->get( i*2   , j*2   );
						pixelData[index +1] = matrix->get( i*2   , j*2+1 );
						pixelData[index +2] = matrix->get( i*2+1 , j*2   );
						pixelData[index +3] = matrix->get( i*2+1 , j*2+1 );
					}
				}

				j = w-1;
				if(nBasisF%2 != 0){
					for(int i=0; i<h; i++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix->get( i*2   , j*2   );
						pixelData[index +1] = fringe;
						if(nOrbitals%2 == 0 || i<h-1)
							pixelData[index +2] = matrix->get( i*2+1 , j*2   );
						else
							pixelData[index +2] = fringe;
						pixelData[index +3] = fringe;
					}
				} else {
					for(int i=0; i<h; i++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix->get( i*2   , j*2   );
						pixelData[index +1] = matrix->get( i*2   , j*2+1 );
						if(nOrbitals%2 == 0 || i<h-1){
							pixelData[index +2] = matrix->get( i*2+1 , j*2   );
							pixelData[index +3] = matrix->get( i*2+1 , j*2+1 );
						} else {
							pixelData[index +2] = fringe;
							pixelData[index +3] = fringe;
						}
					}
				}

				i = h-1;
				if(nOrbitals%2 != 0){
					for(int j=0; j<w; j++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix->get( i*2   , j*2   );
						if(nBasisF%2 == 0 || j<w-1)
							pixelData[index +1] = matrix->get( i*2   , j*2+1 );
						else
							pixelData[index +2] = fringe;
						pixelData[index +2] = fringe;
						pixelData[index +3] = fringe;
					}
				} else {
					for(int j=0; j<w; j++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix->get( i*2   , j*2   );
						pixelData[index +2] = matrix->get( i*2+1 , j*2   );                
						if(nBasisF%2 == 0 || j<w-1){
							pixelData[index +1] = matrix->get( i*2   , j*2+1 );
							pixelData[index +3] = matrix->get( i*2+1 , j*2+1 );
						} else {
							pixelData[index +1] = fringe;
							pixelData[index +3] = fringe;
						}
					}
				}
			}
		}
		if(showTimings){ sw.stop(); cout << " mapping: " << sw.timeMS(); }
		//PrintRGBAPixelsBox(pixelData,numCols*deltaBF,numRows*deltaOE);
		//cout << endl << data << endl;
	}

	void unloadData(ArrayGPU & results){
		gpuData->BeginCapture();
		glReadPixels(0,0,nCalcs*deltaOE,5*deltaOE,GL_RGBA,GL_FLOAT,cpuData);
		gpuData->EndCapture();
		//Array2D<double> temp;
		int rstart, cstart, index;
		for(int r=0; r<5; r++){
			rstart = r*deltaOE;
			for(int c=0; c<nCalcs; c++){
				cstart = c*deltaOE;
				//temp = (results(c*5 + r));
				for(int i=0; i<deltaOE; i++){
					for(int j=0; j<deltaOE; j++){
						index = 4*( (rstart + i)*nCalcs*deltaOE + (cstart + j) );
						(results(c*5 + r))(i*2,j*2)     = cpuData[index   ];
						if(j*2+1 < nOrbitals)
						(results(c*5 + r))(i*2,j*2+1)   = cpuData[index +1];
						if(i*2+1 < nOrbitals)
						(results(c*5 + r))(i*2+1,j*2)   = cpuData[index +2];
						if(i*2+1 < nOrbitals && j*2+1 < nOrbitals)
						(results(c*5 + r))(i*2+1,j*2+1) = cpuData[index +3];
					}
				}
			}
		}
		//PrintRGBAPixelsBox(cpuData,nCalcs*deltaOE,5*deltaOE);
		//cout << endl << results << endl;
	}

	/**The texture id's of the operands*/
	GLuint operands[2];

	/**operand indicies*/
	const static int basis = 0;
	const static int coeff = 1;

	/**indicies for psi, grad psi, and lap psi*/
	const static int psi = 0;
	const static int grx = 1;
	const static int gry = 2;
	const static int grz = 3;
	const static int lap = 4;

	/**Number of simultaneous configurations processed*/
	int nCalcs;

	/**There is no bijection between the dimensions of CPU data and texture dimensions*/
	int nOrbitals, nBasisF;

	/**Internal dimensions of gpuData*/
	int deltaOE, deltaBF;

	/**Compiled version of qmcMatrixProgram*/
	CGprogram fp;

	/**Input parameters of qmcMatrixProgram*/
	CGparameter tELxOR, tELxBF, tORxBF, startOps, stopOps, cg_deltaOE, cg_deltaBF;

	/**Multiplication parameters*/
	int maxLoops, numLoops, numPasses;

	/**This is the GPU data structure meant to hold the result of the calculation*/
	RenderTexture * gpuData;

	/**This is the CPU data structure meant to hold the result of the calculation*/
	GLfloat * cpuData;
	GLfloat * pixelData;
};
#endif