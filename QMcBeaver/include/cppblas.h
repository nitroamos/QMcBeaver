/*
  BLAS might be 0.1 % faster than CBLAS, so choice
  here depends on whether it's easier to compile with
  BLAS or CBLAS.

  configure.py has been set up to assume BLAS.

  If you uncomment this flag, then your library line will probably
  need to look something like:
  LINK = -lm -L$(HOME)/lib -llapack -lf77blas -lcblas -latlas -lg2c 
*/
//#define USE_CBLAS

//this will concatenate an underscore at the end of the function name
//used on gcc/linux, osx, 
#define FORTRAN_FUNC(x) x ## _
//keep the function name as is
//used on AIX
//#define FORTRAN_FUNC(x) (x)

extern "C"
{
#ifdef USE_CBLAS
  /*
    If you want to use CBLAS instead of BLAS, this
    is the header file.
  */
#include "cblas.h"

#else
  /*
    These are the function prototypes for fortran BLAS libraries.
    The calling convention might be different depending on compiler.
    That is, you may need to remove (or add) "_" from the name.
  */
  int FORTRAN_FUNC(dgemm)(const char* transa, const char* transb,
	     const int* m, const int* n, const int* k,
	     const double* alpha, const double* a, const int* lda,
	     const double* b, const int* ldb, const double* beta,
	     double* c, const int* ldc);
  int FORTRAN_FUNC(sgemm)(const char* transa, const char* transb,
	     const int* m, const int* n, const int* k,
	     const float* alpha, const float* a, const int* lda,
	     const float* b, const int* ldb, const float* beta,
	     float* c, const int* ldc);
  double FORTRAN_FUNC(ddot)(const int* N, const double *X, const int* incX,
	       const double *Y, const int* incY);
  float FORTRAN_FUNC(sdot)(const int* N, const float *X, const int* incX,                                                                                                          
	      const float *Y, const int* incY);

  void FORTRAN_FUNC(daxpy)(const int* N, const double *a,
	      const double * x, const int* incx,
	      double *y, const int* incy);
#endif
}
