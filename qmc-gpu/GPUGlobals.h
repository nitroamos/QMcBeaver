#ifndef GPUGLOBALS
#define GPUGLOBALS

#include <windows.h>
#include <gl/glew.h>
#include <gl/glut.h>
#include <cg/cgGL.h>
#include "Array2D.h"

typedef Array1D< Array2D<QMCfloat> > ArrayGPU;

#define QMC_GPU
static const bool showTimings = true;

extern CGcontext g_cgContext;
extern CGprofile g_cgProfile;

/*if there is an error, msg will be printed*/
static void getOpenGLError(char *msg){
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

/*provides a text representation of the pixels:
r	g
b	a  */
static void PrintRGBAPixelsBox(float* pix, int w, int h){
	int index;
	int maxJ = 30;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w && j < maxJ; j++) {
			index = 4*(i*w + j);
			printf("%6.3f",pix[index +0]);
			printf("%6.3f",pix[index +1]);
			printf("   ");
		}
		printf("\n");
		for (int j = 0; j < w && j < maxJ; j++) {
			index = 4*(i*w + j);
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
static void PrintRGBAPixelsColumn(float* pix, int w, int h){
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

static void PrintMatrix(Array2D<GLfloat> matrix){
	for(int i=0; i<matrix.dim1(); i++){
		for(int j=0; j<matrix.dim2() && j < 28; j++){
			printf("%7.3g", (float)matrix(i,j));
		}
		printf("\n");
	}
	printf("\n");
}

/**this method will print out the limitations placed on shader programs by the profile currently loaded*/
static void specTest(){
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