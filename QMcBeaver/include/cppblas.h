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
  int dgemm_(const char* transa, const char* transb,
	     const int* m, const int* n, const int* k,
	     const double* alpha, const double* a, const int* lda,
	     const double* b, const int* ldb, const double* beta,
	     double* c, const int* ldc);
  int sgemm_(const char* transa, const char* transb,
	     const int* m, const int* n, const int* k,
	     const float* alpha, const float* a, const int* lda,
	     const float* b, const int* ldb, const float* beta,
	     float* c, const int* ldc);
  double ddot_(const int* N, const double *X, const int* incX,
	       const double *Y, const int* incY);
  float sdot_(const int* N, const float *X, const int* incX,                                                                                                          
	      const float *Y, const int* incY);

  void daxpy_(const int* N, const double *a,
	      const double * x, const int* incx,
	      double *y, const int* incy);
#endif
}
