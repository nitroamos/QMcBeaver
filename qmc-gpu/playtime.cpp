#include <iostream>
#include <stdlib.h>

#include <windows.h>
#include <gl/glew.h>
#include <gl/glut.h>
#include <cg/cgGL.h>
#include <time.h>

#include "Array2D.h"
#include "Stopwatch.h"
#include "playtime.h"
#include "matrix.h"
#include "dgemm.h"
#include "cppblas.h"
#include "basisfunction.h"

using namespace std;

typedef double mytype;
Matrix todisplay, other;
int rows, columns;

void testBasisFunction(int dim){
	Stopwatch sw = Stopwatch();
	long timeGPU=0, timeCPU=0;

	//there has got to be a "more elegant" way of initializing an Array2D than this...
	//problem is probably just that i'm dumb and forgot pointer stuff
	Array2D<double> coeffs = Array2D<double>(8,2);
	double temp[8][2] = {
		{6.66500000E+03,	0.363584299905},
		{1.00000000E+03,	0.674985792448},
		{2.28000000E+02,	1.1316199392},
		{6.47100000E+01,	1.65300920982},
		{2.10600000E+01,	1.92382064638},
		{7.49500000E+00,	1.44727831645},
		{2.79700000E+00,	0.439163129646},
		{5.21500000E-01,	0.00664580366553}
	};
	for(int i=0; i<coeffs.dim1(); i++)
		for(int j=0; j<coeffs.dim2(); j++)
			coeffs(i,j) = temp[i][j];

	int numBF = dim*dim;
	Array2D<float> R = Array2D<float>(numBF,4);
    for(int i=0; i<R.dim1(); i++)
		for(int j=0; j<R.dim2(); j++)
            R(i,j) = (double)(rand()/33000.0);
	

	BasisFunction bfCPU = BasisFunction(numBF,coeffs,0,0,0,coeffs.dim1(),0.0,0.0,0.0,"purple monkey");
	BasisFunction bfGPU = BasisFunction(numBF,coeffs,0,0,0,coeffs.dim1(),0.0,0.0,0.0,"purple monkey");	

	sw.reset();
	sw.start();
	bfCPU.calculateBasisFunctions(R,false);
	sw.stop();
	timeCPU = sw.timeMS();
	
	int numIters = 10;
	bfGPU.calculateBasisFunctions(R,true);
	for(int i=0; i<numIters; i++){
		sw.reset(); sw.start();
		bfGPU.calculateBasisFunctions(R,true);	
		sw.stop();
		timeGPU += sw.timeMS();
		//cout << i << ") took " << sw.timeMS() << endl;
	}
	if(numIters != 0) timeGPU /= numIters;

	if(false){
		cout << "calculating errors...\n";
		mytype gpuResult;
		mytype cpuResult;
		double avgError = 0;
		double lrgError = 0;
		double error = 0;
		for(int i=0; i<numBF; i++){
			bfGPU.getPsi(i,&gpuResult);
			bfCPU.getPsi(i,&cpuResult);
			if(i < 20 && !true){
				cout << "gpupsi " << gpuResult << " cpupsi " << cpuResult << endl;
			}
			if(abs(gpuResult - cpuResult) == 0) break;
			error = (gpuResult - cpuResult)/cpuResult;
			avgError += error;
			if(error > lrgError) lrgError = error;
		}
		avgError /= numBF;
		printf("average error was %8.4e, and the largest error was %8.4e\n",avgError,lrgError);
	}
	printf("%10i%10i%10i%10i\n",dim,dim*dim,timeCPU,timeGPU);
	coeffs.deallocate();
	R.deallocate();
	bfGPU.destroy();
	bfCPU.destroy();
}

void makeRandMatrix(Array2D<float> &matrix){
    for(int i=0; i<matrix.dim1(); i++)
		for(int j=0; j<matrix.dim2(); j++){
            matrix(i,j) = (float)(rand()/33000.0);
			//matrix(i,j) = (float)( (int)(matrix(i,j)*10000) );
		}
}
void makeRandMatrix(Array2D<double> &matrix){
    for(int i=0; i<matrix.dim1(); i++)
		for(int j=0; j<matrix.dim2(); j++){
            matrix(i,j) = (double)(rand()/33000.0);
			//matrix(i,j) = (double)( (int)(matrix(i,j)*10000) );
		}
}

void PrintArray(Array2D<GLfloat> matrix) {
	for(int i=0; i<matrix.dim1(); i++){
	    for(int j=0; j<matrix.dim2() && j < 28; j++){
		    printf("%12.3g", matrix(i,j));
        }
        printf("\n");
    }
	printf("\n");
}

/**this just does matrix multiply on the GPU. used to make sure the matrix multiply
is doing things correctly*/
void checkMatrixMultiply(){
	Array2D<GLfloat> C,D;
	Stopwatch sw = Stopwatch();
	long timeGPU=0, timeCPU=0;
	int numIters = 100, M, N, K;
	bool verbose = false;
	bool errorCheck = true;
	bool dispTimings = false;
	
	/*
	K = 1000;
	M = K;
	N = K;
	/*/
	K = 570;
	M = 76;
	N = 76;
	//*/
	
	cout << "beginning matrix multiply test: C(" << M << ", " << N << ") = A(" << M << ", " 
		 << K << ") * B(" << K << ", " << N << ")...\n\n";
	
	Array2D<GLfloat> A = Array2D<GLfloat>(M,K);
	makeRandMatrix(A);
	Array2D<GLfloat> B = Array2D<GLfloat>(K,N);
	makeRandMatrix(B);

	Matrix Amat = Matrix(A);
	Matrix Bmat = Matrix(B);
	Matrix Cmat = Matrix(M,N);

	Amat = Matrix(A);
	sw.reset(); sw.start();
	for(int i=0; i<numIters; i++)
		Amat = Matrix(A);
	sw.stop();
	cout << "average time loading " << M << "x" << K << " matrix was " << sw.timeMS()/(double)numIters << endl << endl;

	Amat.matrixMultiplyDoNothing(Bmat, Cmat);
	sw.reset(); sw.start();
	for(int i=0; i<1; i++)
		Amat.matrixMultiplyDoNothing(Bmat, Cmat);	
	sw.stop();
	cout << "average time doing nothing " << M << "x" << N << " matrix was " << sw.timeMS()/(double)numIters << endl << endl;

	Amat.matrixMultiplyFaster(Bmat, Cmat);	
	sw.reset(); sw.start();
	for(int i=0; i<1; i++)
		Amat.matrixMultiplyFaster(Bmat, Cmat);	
	sw.stop();
	timeGPU = sw.timeMS();
	timeGPU /= (float)(numIters+0.5);
	cout << "average time multiplying " << M << "x" << N << " matrix was " << timeGPU << endl << endl;

	C = Cmat.getData();	
	sw.reset(); sw.start();
	for(int i=0; i<numIters; i++)
		C = Cmat.getData();	
	sw.stop();
	cout << "average time unloading " << M << "x" << N << " matrix was " << sw.timeMS()/(double)numIters << endl << endl;
	
	if(verbose){
		cout << "matrix A (" << A.dim1() << ", " << A.dim2() << ") is: " << endl;
		PrintArray(A);
		cout << "matrix B (" << B.dim1() << ", " << B.dim2() << ") is: " << endl;
		PrintArray(B);
		cout << "matrix C (" << C.dim1() << ", " << C.dim2() << ") is: " << endl;
		PrintArray(C);
	}

	if(errorCheck){
		sw.reset(); sw.start();	
		for(int i=0; i<numIters; i++)			
			D = A*B;
		sw.stop();
		timeCPU = sw.timeMS();
		cout << "average time CPU " << M << "x" << N << " matrix was " << sw.timeMS()/(double)numIters << endl << endl;
		timeCPU /= (float)(numIters+0.5);

		bool same = true;
		int largestI = 0, largestJ = 0;
		double currentError, largestError = 0, error = 0;
		for(int i=0; i<C.dim1(); i++)
			for(int j=0; j<C.dim2(); j++){
				currentError = abs( (C(i,j)-D(i,j))/D(i,j) );
				if(currentError > largestError){
					largestError = currentError;
					largestI = i; largestJ = j;
				}
				error += currentError;
				if( currentError > 0.0001){
					same = false;
				}
			}
		cout << "largest error at (" << largestI << ", " << largestJ << ") with GPU,CPU "
			 << C(largestI,largestJ) << ", " << D(largestI,largestJ) << " with error " << largestError << endl;
		error = error/(C.dim1()*C.dim2());
		cout << "average rel error is " << error << endl;
		printf("GPU time = %6i\nCPU time = %6i\nfactor = %8.4g\n",timeGPU,timeCPU,(double)timeCPU/(double)timeGPU);
		if(same){
			cout << "(" << D.dim1() << ", " << D.dim2() << ") dimensional matrix multiplication a success" << endl;
		} else {
			cout << "error in the calculation" << endl;
			if(verbose && false){
				cout << "matrix supposed to be: " << endl;
				PrintArray(D);
			}
		}
	} 

	rows = Cmat.getRows();
	columns = Cmat.getColumns();
	perspective(false);
	todisplay = Cmat;
	todisplay.display();
}

void basic_dgemmF (const int lda, const int M, const int N, const int K,
			 const mytype *A, const mytype *B, mytype *C)
{
	int i, j, k;
	for (i = 0; i < M; ++i) {
		const register mytype *Ai_ = A + i*lda;
		for (j = 0; j < N; ++j) {
			const register mytype *B_j = B + j*lda;
			register mytype cij = C [j*lda + i];
			for (k = 0; k < K; ++k) {
				cij += Ai_[k] * B_j[k];
			}
			C[j*lda + i] = cij;
		}
	}
}

/**this prints timings for all matrix multiply methods available*/
void matrixMultiplyTimings(int d){
	Stopwatch sw = Stopwatch();
	long time, basic;
	printf("%5i",d);

	double * E = new double[d * d];
	double * F = new double[d * d];
	double * G = new double[d * d];
	for(int i=0; i<d*d; i++){
		E[i] = (double)(rand()/33000.0);
		F[i] = (double)(rand()/33000.0);
		G[i] = (double)0;
	}
	sw.reset();
	sw.start();
	basic_dgemm(d,d,d,d,E,F,G);
	sw.stop();
	basic = sw.timeMS();
	printf("%20i",basic);

	Array2D<mytype> A = Array2D<mytype>(d,d),B = Array2D<mytype>(d,d),C;
	makeRandMatrix(A);
	makeRandMatrix(B);
	sw.reset();
	sw.start();
	C = A*B;
	sw.stop();
	printf("%20i",sw.timeMS());
	A.deallocate(); B.deallocate(); C.deallocate();

	if(!true){
		for(int l=32; l<=1024; l*=2)
			for(int s=32; s<=l; s*=2){
				set_l_block_size(l);
				set_s_block_size(s);
				sw.reset();
				sw.start();
				square_dgemm(d,E,F,G);
				sw.stop();
			}
	} else {
		int l = 256, s = 256;
		set_l_block_size(l);
		set_s_block_size(s);
		sw.reset();
		sw.start();
		square_dgemm(d,E,F,G);
		sw.stop();
		printf("%20i",sw.timeMS());
	}

	sw.reset();
	sw.start();
	//for some reason, even though this line works, if i uncomment it, the debugger can't work! :-(
	cblas_dgemm(CBLAS_ORDER(CblasRowMajor),CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasNoTrans),d,d,d,1.0,E,d,F,d,0.0,G,d);
	sw.stop();
	printf("%20i",sw.timeMS());
	delete[] E; delete[] F;	delete[] G;

	Matrix Amat = Matrix(d,d);
	Matrix Bmat = Matrix(d,d);
	Matrix Cmat = Matrix(d,d);
	int numIters = 5;
	long timeGPU = 0;
	Amat.matrixMultiplyFaster(Bmat, Cmat);
	for(int i=0; i<numIters; i++){
		sw.reset(); sw.start();
		Amat.matrixMultiplyFaster(Bmat, Cmat);	
		sw.stop();
		timeGPU += sw.timeMS();
		//cout << i << ") took " << sw.timeMS() << endl;
	}
	timeGPU /= numIters;
	printf("%20i\n",timeGPU);
	Amat.destroy(); Bmat.destroy(); Cmat.destroy();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
		case 27://esc
			exit(0);
			break;
		case '1':
			printf("%10s%10s%10s%10s\n","Dim","#BF","CPU","GPU");
			if(!true){
				for(int i=100; i<=1500; i+=100)
					testBasisFunction(i);
			} else {
				testBasisFunction(1000);
			}
			break;
		case '2':
			checkMatrixMultiply();			
			break;
		case '3':
			printf("%5s%20s%20s%20s%20s%20s\n","dim","simple (pointers)","Array2D","DGEMM","ATLAS","GPU");
			if(true){			
				for(int i=300; i<=1000; i+=100)
					matrixMultiplyTimings(i);			
			} else {
				matrixMultiplyTimings(1000);
			}	
			cout << "done with matrixMultiplyTimings" << endl;
			break;
		case 'a':
			todisplay += 0.1;
			todisplay.display();
			break;
		case 'm':
			todisplay *= 0.9;
			todisplay.display();
			break;
		case 'n':
			todisplay = todisplay * 0.9;
			todisplay.display();
			break;
		case 'r':
			todisplay = Matrix(rows,columns,1.0,8.0f);
			todisplay.display();
			break;
		case 's':
			todisplay = 0.5;
			todisplay.display();
			break;
		case 't':
			todisplay += other;
			todisplay.display();
			break;
		case 'u':
		case 'U':
			todisplay.unloadMatrix(true);
			break;
		case '='://+
			rows *= 2;
			columns *= 2;
			perspective(false);
			todisplay = Matrix(rows,columns,1.0,8.0f);
			other = Matrix(rows,columns,0.3,16.0f);
			todisplay.display();
			break;
		case '-':
			rows /= 2;
			columns /= 2;
			perspective(false);
			todisplay = Matrix(rows,columns,1.0,8.0f);
			other = Matrix(rows,columns,0.3,16.0f);
			todisplay.display();
			break;
	}
}

	// Called when Cg detects an error
void cgErrorCallback(){
	CGerror lastError = cgGetError();
 
	if(lastError){
		printf("%s\n\n", cgGetErrorString(lastError));
		printf("%s\n", cgGetLastListing(g_cgContext));
		printf("Cg error!\n");
	}
} 

// Called at startup
void initializeCg(){
	cgSetErrorCallback(cgErrorCallback);

	g_cgContext = cgCreateContext();
    
	// get the best profile for this hardware
	g_cgProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
	if(g_cgProfile == CG_PROFILE_FP40) cout << "Using fp40\n";
	if(g_cgProfile == CG_PROFILE_FP30) cout << "Using fp30\n";
	if(g_cgProfile == CG_PROFILE_FP20) cout << "Using fp20\n";
	assert(g_cgProfile != CG_PROFILE_UNKNOWN);
	cgGLSetOptimalOptions(g_cgProfile);	
}

/**this method will print out the limitations placed on shader programs by the profile currently loaded*/
void specTest(){
	
	if(!GL_NV_fragment_program2) cout << "no GL_NV_fragment_program2" << endl;
	if(!GL_FRAGMENT_PROGRAM_NV) cout << "no GL_FRAGMENT_PROGRAM_NV" << endl;
	if(!GL_FRAGMENT_PROGRAM_ARB) cout << "no GL_FRAGMENT_PROGRAM_ARB" << endl;
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
	getError("Error in glGetProgramivARB");
	cout << endl;


	/*
	glGetIntegerv(GL_MAX_FRAGMENT_PROGRAM_LOCAL_PARAMETERS_NV, &num );
	cout << "GL_MAX_FRAGMENT_PROGRAM_LOCAL_PARAMETERS_NV is " << num << endl;
	glGetIntegerv(GL_MAX_PROGRAM_EXEC_INSTRUCTIONS_NV, &num );
	cout << "GL_MAX_PROGRAM_EXEC_INSTRUCTIONS_NV is " << num << endl;
	cout << endl;
	getError("Error in glGetIntegerv");
	*/
}

void init(){
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
    glutInitWindowSize(250, 250);
    glutInitWindowPosition(500, 500);
	glutCreateWindow("Hello, GPGPU!");
	
	glDrawBuffer(GL_BACK);
	glReadBuffer(GL_BACK);
	glDisable (GL_DEPTH_TEST);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glDrawBuffer(GL_FRONT);
	glReadBuffer(GL_FRONT);
	glDisable(GL_DEPTH_TEST);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	int err = glewInit();
    if (GLEW_OK != err)
    {
        // problem: glewInit failed, something is seriously wrong
        fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
        getchar();
		exit(-1);
    }

    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc (keyboard);
	getError(0);
}

/*if there is an error, msg will be printed*/
void getError(char *msg){
    static char	error_txt[][32]	=	{
        "GL_INVALID_ENUM",
        "GL_INVALID_VALUE",
        "GL_INVALID_OPERATION",
        "GL_STACK_OVERFLOW",
        "GL_STACK_UNDERFLOW",
        "GL_OUT_OF_MEMORY" };
    GLenum e = glGetError();
	if (e != GL_NO_ERROR){
        cout << "OpenGL error " << error_txt[e-GL_INVALID_ENUM] << endl;
		if(msg != NULL) cout << "Message: " << msg << endl;
	}

}

int main(int argc, char* argv[])
{
	init();
	initializeCg();
	srand( (unsigned)time( NULL ) );
	rows = 10;
	columns = 20;

	todisplay = Matrix(rows,columns,1.0,8.0f);
	other = Matrix(rows,columns,0.3,16.0f);

	//there's no way of avoiding glutMainLoop is there?
	//e.g. if i could manually show the window and manually proceed to manupulate it
	//i guess for the purposes of programming, this needs to rely on the callback functions
	//to perform any calculations.
	glutMainLoop();
	getchar();
	return 0;
}

void idle(){
	//glutPostRedisplay();
}

// GLUT display function
void display(){
	todisplay.display();
}

// GLUT reshape function
void reshape(int w, int h){	
	window_w = w;
	window_h = h;
	perspective(false);
}

//must call this function if the size of the matrix has been changed
//if transformed = true, then the texture will appear inset in the window (prettier?)
void perspective(bool transformed){
	if(transformed){
		if (window_h == 0) window_h = 1;
		glViewport(0, 0, window_w, window_h);
		glMatrixMode(GL_PROJECTION);    
		glLoadIdentity();               
		gluOrtho2D(-1, 1, -1, 1); 
		gluPerspective(60.0, (GLfloat) window_w/(GLfloat) window_h, 1.0, 30.0);
		//glTranslatef(1.0, -1.0, 0.0);
		glMatrixMode(GL_MODELVIEW);     
		glLoadIdentity();
		glTranslatef(0.0, 0.0, -2.0);;
	} else {
		glViewport(0, 0, (GLsizei)(columns/2), (GLsizei)(rows/2));
		glMatrixMode(GL_PROJECTION);    
		glLoadIdentity();               
		gluOrtho2D(-1, 1, -1, 1);       
		glMatrixMode(GL_MODELVIEW);     
		glLoadIdentity(); 
	}
}