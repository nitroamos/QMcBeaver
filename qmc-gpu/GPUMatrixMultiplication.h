#include "RenderTexture.h"
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <windows.h>
#include <gl/glut.h>
#include <cg/cgGL.h>
#include <math.h>

#include "Array1D.h"
#include "Array2D.h"
#include "Stopwatch.h"

static const char * qmcMatrixProgram =
"float4 main(in float2 position : TEX0,                              \n"
"            uniform int  delta,                                     \n"
"            uniform int  startOps,                                  \n"
"            uniform int  stopOps,                                   \n"
"            uniform samplerRECT  accum,                             \n"
"            uniform samplerRECT  texLHS,                            \n"
"            uniform samplerRECT  texRHS) : COLOR                    \n"
"{                                                                   \n"
"    int c = position.x/delta;                                       \n"
"    int shift = c*delta                                             \n"
"    float4 output = texRECT(accum, position);                       \n"
"    for (int i = startOps; i < stopOps; i++) {                      \n"
"       float4 a  = texRECT(texLHS, float2(i + shift, position.y));  \n"
"       float4 b  = texRECT(texRHS, float2(i, position.y));          \n"
"       output += a.xxzz*b.xzxz;                                     \n"
"       output += a.yyww*b.ywyw;                                     \n"
"    }                                                               \n"
"    return output;                                                  \n"
"}                                                                   \n";

typedef Array1D< Array2D<double> > ArrayGPU;

/**The Coeff and basisFunctions submatricies have the same shape because the
Coeff matrix is passed in transposed.*/
class GPUQMCMatrix {
public:
	GPUQMCMatrix(){
		nOrbitals = 0;
		nBasisF = 0;
		nCalcs = 0;
	}

	GPUQMCMatrix(Array2D<double> & Coeffs, int numCalcs){
		nCalcs = numCalcs;
		nOrbitals = Coeffs.dim1();
		nBasisF = Coeffs.dim2();
		glGenTextures(2, operands);
		loadCoeffs(Coeffs);
		initialize();
	}

	void destroy(){
		delete [] cpuData;
		gpuData->Reset(0,0);
	}

	ArrayGPU runCalculation(ArrayGPU basisFunctions){
		if(basisFunctions.dim1() != nCalcs) cout << "wrong num calcs";
		loadBasisF(basisFunctions);
		cleanRenderTexture();

		gpuData->BeginCapture();
		cgGLBindProgram(fp);
		cgGLEnableProfile(g_cgProfile);

		if(numPasses%2 == 1) gpuData->swapBuffers();

		cgGLSetParameter1f(delta, deltaOE);
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
			
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
			gpuData->EndCapture();			
		}
		glFinish();
		glFlush();
		
		cgGLDisableTextureParameter(tELxBF);
		cgGLDisableTextureParameter(tORxBF);
		cgGLDisableTextureParameter(tELxOR);
		cgGLDisableProfile(g_cgProfile);
		getError("Error in QMC matrix multiply");
	}

	double getPsi(int whichCalc, int whichElec, int whichOrb){
		return resultMapping(psi,whichCalc,whichElec,whichOrb);		
	}

	double getGrad(int whichCalc, int xyz, int whichElec, int whichOrb){
		return resultMapping(grx+xyz,whichCalc,whichElec,whichOrb);		
	}

	double getLaplacian(int whichCalc, int whichElec, int whichOrb){
		return resultMapping(lap,whichCalc,whichElec,whichOrb);
	}

private:
	void initialize(){
		deltaBF = nBasisF/2.0;
		deltaOE = nOrbitals/2.0;
		if(nBasisF%2 != 0)   deltaBF += 1;
		if(nOrbitals%2 != 0) deltaOE += 1;

		gpuData = new RenderTexture("rgba=32f doublebuffer texRECT rtt");
		gpuData->Initialize(nCalcs*deltaOE,5*deltaOE,true,false);
		cpuData = new GLfloat[ nCalcs*deltaOE * 5*deltaOE * 4 ];

		fp = cgCreateProgram(g_cgContext, CG_SOURCE,qmcMatrixProgram, g_cgProfile,"main", NULL);
        //fp = cgCreateProgram(g_cgContext, CG_SOURCE,testInputs, g_cgProfile,"main", NULL);

		if(fp != NULL){
            cgGLLoadProgram(fp);
			tELxOR   = cgGetNamedParameter(fp, "accum");     
            tELxBF   = cgGetNamedParameter(fp, "texLHS");
			tORxBF   = cgGetNamedParameter(fp, "texRHS");
			startOps = cgGetNamedParameter(fp, "startOps");
			stopOps  = cgGetNamedParameter(fp, "stopOps");
			delta    = cgGetNamedParameter(fp, "delta");
		} else {
			cout << "error in matrixMultiply script" << endl;
			return;
		}

		numLoops = deltaOE;

		//*
		maxLoops = 255;
		/*/
		glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &maxLoops);
		//*/
		
		//*
		numPasses = ceil((double)numLoops/maxLoops);
		/*/
		numPasses = 2;
		//*/
	}

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

	void loadCoeffs(Array2D<double> & Coeffs){
		ArrayGPU coeffA = ArrayGPU(1);
		coeffA(0) = Coeffs;
		glActiveTextureARB(GL_TEXTURE1_ARB);
		glBindTexture(GL_TEXTURE_RECTANGLE_NV, operands[coeff]);
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
					 deltaBF, deltaOE, 0, GL_RGBA, GL_FLOAT, 
					 mapData(coeffA,true));
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        getError("Error loading Coeffs");
	}

	void loadBasisF(ArrayGPU basisFunctions){
		glActiveTextureARB(GL_TEXTURE0_ARB);
		glBindTexture(GL_TEXTURE_RECTANGLE_NV, operands[basis]);
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
					 nCalcs*deltaBF, 5*deltaOE, 0, GL_RGBA, GL_FLOAT, 
					 mapData(basisFunctions,false));
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        getError("Error loading basis function data");
	}

	/**r and c indicate which submatrix we're in, i and j indicate which value in that submatrix*/
	int operandMapping(int r, int c, int i, int j, int numCols){
		//4*( (row selection)*(total width) + (column selection) );
		return 4*( (r*deltaOE + i)*numCols*deltaBF + (c*deltaBF + j) );
	}

	/**r and c indicate which submatrix we're in, i and j indicate which value in that submatrix*/
	int resultMapping(int r, int c, int i, int j){
		//4*( (row selection)*(total width) + (column selection) );
		return 4*( (r*deltaOE + i)*nCalcs*deltaOE + (c*deltaOE + j) );
	}

	GLfloat * mapData(ArrayGPU &data, bool isCoeff){
		int numRows = isCoeff? 1 : 5;
		int numCols = isCoeff? 1 : nCalcs;
		int h = deltaOE;
		int w = deltaBF;
		GLfloat * pixelData = new GLfloat[ numRows*deltaOE * numCols*deltaBF * 4 ];
		GLfloat fringe = 0.0;
		Array2D<double> matrix;

		//we iterate over all the matricies		
		for(int c = 0; c < numCols; c++){
			for(int r = 0; r < numRows; r++){
				matrix = data( c*numRows + r);
				
				int i, j, index;
				for(i=0; i<h-1; i++){
					for(j=0; j<w-1; j++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix.get( i*2   , j*2   );
						pixelData[index +1] = matrix.get( i*2   , j*2+1 );
						pixelData[index +2] = matrix.get( i*2+1 , j*2   );
						pixelData[index +3] = matrix.get( i*2+1 , j*2+1 );
					}
				}

				j = w-1;
				if(nBasisF%2 != 0){
					for(int i=0; i<h; i++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix.get( i*2   , j*2   );
						pixelData[index +1] = fringe;
						if(nOrbitals%2 == 0 || i<h-1)
							pixelData[index +2] = matrix.get( i*2+1 , j*2   );
						else
							pixelData[index +2] = fringe;
						pixelData[index +3] = fringe;
					}
				} else {
					for(int i=0; i<h; i++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix.get( i*2   , j*2   );
						pixelData[index +1] = matrix.get( i*2   , j*2+1 );
						if(nOrbitals%2 == 0 || i<h-1){
							pixelData[index +2] = matrix.get( i*2+1 , j*2   );
							pixelData[index +3] = matrix.get( i*2+1 , j*2+1 );
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
						pixelData[index   ] = matrix.get( i*2   , j*2   );
						if(nBasisF%2 == 0 || j<w-1)
							pixelData[index +1] = matrix.get( i*2   , j*2+1 );
						else
							pixelData[index +2] = fringe;
						pixelData[index +2] = fringe;
						pixelData[index +3] = fringe;
					}
				} else {
					for(int j=0; j<w; j++){
						index = operandMapping(r,c,i,j,numCols);
						pixelData[index   ] = matrix.get( i*2   , j*2   );
						pixelData[index +2] = matrix.get( i*2+1 , j*2   );                
						if(nBasisF%2 == 0 || j<w-1){
							pixelData[index +1] = matrix.get( i*2   , j*2+1 );
							pixelData[index +3] = matrix.get( i*2+1 , j*2+1 );
						} else {
							pixelData[index +1] = fringe;
							pixelData[index +3] = fringe;
						}
					}
				}
			}
		}
		return pixelData;
	}

	Array2D<double> getSlaterMatrix(int whichCalc){
		Array2D<double> matrix = Array2D<double>(nOrbitals,nOrbitals);
		for(int i=0; i<matrix.dim1(); i++)
			for(int j=0; j<matrix.dim2(); j++)
				matrix(i,j) = getPsi(whichCalc,i,j);
		return matrix;
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
	CGparameter tELxOR, tELxBF, tORxBF, startOps, stopOps, delta;

	/**Multiplication parameters*/
	int maxLoops, numLoops, numPasses;

	/**This is the GPU data structure meant to hold the result of the calculation*/
	RenderTexture * gpuData;

	/**This is the CPU data structure meant to hold the result of the calculation*/
	GLfloat * cpuData;
};