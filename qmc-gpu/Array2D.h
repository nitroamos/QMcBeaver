//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   

Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/

#ifndef Array2D_H
#define Array2D_H

#include "Stopwatch.h"
#include <iostream>
#include "cppblas.h"
#include <assert.h>


using namespace std;
//#define USEATLAS

/**
A 2-dimensional template for making arrays.  All of the memory allocation
and deallocation details are dealt with by the class.
*/


template <class T> class Array2D
{
private:
	/**
	Number of elements in the array's first dimension.
	*/

	int n_1;


	/**
	Number of elements in the array's second dimension.
	*/

	int n_2;

	/**
	Array containing the data.
	*/
	//T** pArray;
	T * pArray;

public:
	/**
	Gets the number of elements in the array's first dimension.

	@return number of elements in the array's first dimension.
	*/
	int dim1(){return n_1;}

	/**
	Gets the number of elements in the array's second dimension.

	@return number of elements in the array's second dimension.
	*/
	int dim2(){return n_2;}

	/**
	Gets the total number of elements in the array.

	@return total number of elements in the array.
	*/
	int size(){return n_1*n_2;}

	/**
	Gets a pointer to an array containing the array elements.  
	The ordering of this array is NOT specified.  
	*/
	T* array(){ return pArray; }
	/**
	Allocates memory for the array.

	@param i size of the array's first dimension.
	@param j size of the array's second dimension.
	*/

	void allocate(int i, int j){
		if( n_1 != i || n_2 != j ){
			deallocate();

			n_1 = i;
			n_2 = j;

			if(n_1 >= 1 && n_2 >= 1){
				pArray = new T[ n_1*n_2 ];
			} else {
				pArray = 0;
			}
		}
	}


	/**
	Deallocates memory for the array.
	*/

	void deallocate()
	{
		delete [] pArray;
		pArray = 0;

		n_1 = 0;
		n_2 = 0;
	}


	/**
	Sets two arrays equal.
	*/

	void operator=(const Array2D & rhs)
	{
		if(n_1 != rhs.n_1 || n_2 != rhs.n_2) allocate(rhs.n_1,rhs.n_2);
		memcpy(pArray, rhs.pArray, sizeof(T)*n_1*n_2);
	}


	/**
	Sets all of the elements in an array equal to the same value.
	*/

	void operator=(const T C)
	{
		if(C == 0) {
			memset(pArray,0,sizeof(T)*n_1*n_2);
			return;
		}
		for(int i=0; i<n_1*n_2; i++)
			pArray[i] = C;
	}


	/**
	Returns the matrix product of two arrays.
	*/

#ifdef USEATLAS
	Array2D<double> operator*(const Array2D<double> & rhs)
	{
		if(n_2 != rhs.n_1)
		{
			cerr << "ERROR: Matrices of incorrect dimensions are being"
				<< " multiplied!" << endl;
			exit(1);
		}
		Array2D<T> TEMP(n_1,rhs.n_2);
		TEMP = 0;
/*
	Array2D<int> A = Array2D<int>(d-1,d);
	Array2D<int> B = Array2D<int>(d,d+1);
	for(int i=0; i<A.dim1(); i++)
		for(int j=0; j<A.dim2(); j++)
			A(i,j) = i + 2*j;
	for(int i=0; i<B.dim1(); i++)
		for(int j=0; j<B.dim2(); j++)
			B(i,j) = i + 2*j;
	Array2D<int> C = A*B;

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,const enum CBLAS_TRANSPOSE TransB,
				 const int M, const int N, const int K,
				 const double alpha, const double *A, const int lda,
				 const double *B, const int ldb,
				 const double beta, double *C, const int ldc);
*/
		cblas_dgemm(CBLAS_ORDER(CblasRowMajor),CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasNoTrans),
			n_1, rhs.n_2, n_2,
			1.0, pArray, n_2,
			rhs.pArray, rhs.n_2,
			0.0, TEMP.pArray, TEMP.n_2);
		return TEMP;
	}

	Array2D<float> operator*(const Array2D<float> & rhs)
	{
		if(n_2 != rhs.n_1)
		{
			cerr << "ERROR: Matrices of incorrect dimensions are being"
				<< " multiplied!" << endl;
			exit(1);
		}
		Array2D<T> TEMP(n_1,rhs.n_2);
		TEMP = 0;
		cblas_sgemm(CBLAS_ORDER(CblasRowMajor),CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasNoTrans),
			n_1, rhs.n_2, n_2,
			1.0, pArray, n_2,
			rhs.pArray, rhs.n_2,
			0.0, TEMP.pArray, TEMP.n_2);
		return TEMP;
	}
#else
	Array2D operator*(const Array2D & rhs)
	{
		if(n_2 != rhs.n_1)
		{
			cerr << "ERROR: Matrices of incorrect dimensions are being"
				<< " multiplied!" << endl;
			exit(1);
		}
		
		Array2D<T> TEMP(n_1,rhs.n_2);
		TEMP = 0;
		T * A = pArray;
		T * B = rhs.pArray;
		T * C = TEMP.pArray;
		int M = n_1;
		int N = rhs.n_2;
		int K = n_2;
		for (int i = 0; i < M; ++i) {
			register T *Ai_ = A + i*K;
			for (int j = 0; j < N; ++j) {
				register T cij = 0;
				for (int k = 0; k < K; ++k) {
					cij += Ai_[k] * B[k*N+j];
				}
				C[i*N + j] = cij;
			}
		}
		return TEMP;
	}
#endif


	/**
	Returns the product of an array and a scalar.
	*/

	Array2D operator*(const T C)
	{
		Array2D<T> TEMP(n_1,n_2);
		for(int i=0; i<n_1*n_2; i++)
			TEMP.pArray[i] = C*pArray[i];
		return TEMP;
	}


	/**
	Sets this array equal to itself times a scalar value.
	*/

	void operator*=(const T C)
	{
		for(int i=0; i<n_1*n_2; i++)
			pArray[i] *= C;
	}


	/**
	Sets this array equal to itself divided by a scalar value.
	*/

	void operator/=(const T C)
	{
		T inv = 1.0/C;
		operator*=(inv);
	}


	/**
	Creates an array.
	*/

	Array2D(){pArray = 0; n_1 = 0; n_2 = 0;}


	/**
	Creates an array and allocates memory.

	@param i size of the array's first dimension.
	@param j size of the array's second dimension.
	*/

	Array2D(int i, int j) {
		pArray = 0;
		n_1 = 0; n_2 = 0;
		allocate(i,j);
	}


	/**
	Creates an array and sets it equal to another array.

	@param rhs array to set this array equal to.
	*/

	Array2D( const Array2D<T> & rhs)
	{
		n_1 = 0;
		n_2 = 0;
		pArray = 0;
		allocate(rhs.n_1,rhs.n_2);
		operator=(rhs);
	}


	/**
	Destroy's the array and cleans up the memory.
	*/

	~Array2D(){deallocate();}

	inline int map(int i, int j) const{
		return n_2*i + j;
		//return n_1*j + i;
	}

	/**
	Accesses element <code>(i,j)</code> of the array.
	*/
	T& operator()(int i,int j){
		assert(pArray != 0);
		return pArray[map(i,j)];
	}

	T get(int i, int j) const{
		assert(pArray != 0);
		return pArray[map(i,j)];
	}

	/**
	Prints the array to a stream.
	*/

	friend ostream& operator<<(ostream & strm, const Array2D<T> & rhs)
	{
		for(int i=0; i<rhs.n_1;i++)
		{
			for(int j=0; j<rhs.n_2; j++)
			{
				strm << rhs.pArray[rhs.map(i,j)] << "\t";
			}
			strm << endl;
		}
		return strm;
	}
};

#endif
