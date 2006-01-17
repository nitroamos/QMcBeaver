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

#include <iostream>
#include "cppblas.h"
#include <assert.h>
#include "Array1D.h"
#include "math.h"

#ifdef USELAPACK
#include "cpplapack.h"
#endif

//Array2D is included in all the classes where this change is relevant
#if defined SINGLEPRECISION
typedef float  qmcfloat;
#else
typedef double qmcfloat;
#endif

#if defined SINGLEPRECISION || defined QMC_GPU
#define REALLYTINY 1e-35f
#else
#define REALLYTINY 1e-250
#endif

using namespace std;

/**
   A 2-dimensional template for making arrays.  All of the memory allocation
   and deallocation details are dealt with by the class. Some of the operators
   require the Array2D to contain doubles or floats.
   
   NOTICE: All Array2D a operator * Array2D b methods now require b to be transposed.
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
    T * pArray;

  public:
    /**
    Gets a pointer to an array containing the array elements.  
    The ordering of this array is NOT specified.  
    */
    T* array();

    /**
    Gets the number of elements in the array's first dimension.

    @return number of elements in the array's first dimension.
    */
    int dim1()
    {
      return n_1;
    }

    /**
    Gets the number of elements in the array's second dimension.

    @return number of elements in the array's second dimension.
    */
    int dim2()
    {
      return n_2;
    }

    /**
    Gets the total number of elements in the array.

    @return total number of elements in the array.
    */
    int size()
    {
      return n_1*n_2;
    }

    /**
    Allocates memory for the array.

    @param i size of the array's first dimension.
    @param j size of the array's second dimension.
    */

    void allocate(int i, int j)
    {
      if( n_1 != i || n_2 != j )
        {
          deallocate();

          n_1 = i;
          n_2 = j;

          if(n_1 >= 1 && n_2 >= 1)
            {
              pArray = new T[ n_1*n_2 ];
            }
          else
            {
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
       The flag USEATLAS dictates whether Array2D will use the library to speed it's math or not.
       Actually, if the dimensions of the data are small, using ATLAS may actually be slower.
       FYI: DGEMM stands for Double-precision GEneral Matrix-Matrix multiplication
       DGEMM convention: MxN = MxK * KxN; lda, ldb, and ldc are the n_2 of their respective Array2Ds
    */

#ifdef USEATLAS
    /**This matrix multiplication requires rhs to be transposed.*/
    Array2D<double> operator*(const Array2D<double> & rhs);

    /**This matrix multiplication requires rhs to be transposed.*/
    Array2D<float> operator*(const Array2D<float> & rhs);

    /**rhs is NOT transposed*/
    Array2D<double> multiply(const Array2D<double> & rhs);

    /**rhs is NOT transposed*/
    Array2D<float> multiply(const Array2D<float> & rhs);

    /**
    This uses ATLAS to calculate a dot product between all the elements in one
    Array2D and all the elements in another Array2D
    */
    double dotAllElectrons(const Array2D<double> & rhs);

    /**
    This uses ATLAS to calculate a dot product between all the elements in one
    row from one Array2D and one row in another Array2D
    */
    double dotOneElectron(const Array2D<double> & rhs, int whichElectron);

    float dotAllElectrons(const Array2D<float> & rhs);

    float dotOneElectron(const Array2D<float> & rhs, int whichElectron);

#else

    /**This requires rhs to be transposed*/
    Array2D operator*(const Array2D & rhs);

    /**rhs is NOT transposed*/
    Array2D multiply(const Array2D & rhs);

    /**
    Calculates a dot product between all the elements in one
    Array2D and all the elements in another Array2D
    */
    T dotAllElectrons(const Array2D<T> & rhs);

    /**
    Calculates a dot product between all the elements in one
    row from one Array2D and one row in another Array2D
    */
    T dotOneElectron(const Array2D<T> & rhs, int whichElectron);
#endif

    /**
    Calculates the inverse and determinant of a matrix using a 
    dense LU solver.  This 
    method scales as \f$O(1 N^3)\f$.

    @param a a \f$N \times N\f$ matrix 
    @param inv inverse of a is returned here
    @param det determinant of a is returned here
    @param calcOK returns false if the calculation is singular and true otherwise
    */

    void determinant_and_inverse(Array2D<T> &inv, double& det, bool *calcOK);

    /**
    LU decomposition using the algorithm in numerical recipes for a dense matrix.


    @param a a \f$N \times N\f$ matrix which is destroyed during the operation.
    The resulting LU decompositon is placed here.
    @param indx a \f$N\f$ dimensional array which records the row permutation 
    from partial pivoting.
    @param d used to give det(a) the correct sign
    @param calcOK returns false if the calculation is singular and true otherwise
    */

    void ludcmp(int *indx, double *d, bool *calcOK);


    /**
    LU backsubstitution using the algorithm in numerical
    recipes for a dense matrix. 

    @param a the LU decomposition of a matrix produced by ludcmp
    @param indx a \f$N\f$ dimensional array which records the row permutation 
    from partial pivoting generated by ludcmp
    @param b the \f$N\f$ dimensional array right hand side of the system of 
    equations to solve
    */

    void lubksb(int *indx, Array1D<T>& b);

    /**
    Sets two arrays equal. memcpy would probably break if T is a class. is this a problem?
    if you want to set to Array2D's containing a class equal, then use the (i,j) operator.
    */
    void operator=(const Array2D & rhs)
    {
      if(n_1 != rhs.n_1 || n_2 != rhs.n_2) allocate(rhs.n_1,rhs.n_2);
      memcpy(pArray, rhs.pArray, sizeof(T)*n_1*n_2);
    }

    /**
    Sets all of the elements in an array equal to the same value. There's an obvious shortcut if
    the constant happens to be 0. in retrospect, while memset works well for doubles, floats, etc,
    using this operator would probably crash the program if T was a class.
    */
    void operator=(const T C)
    {
      if(C == 0)
        {
          memset(pArray,0,sizeof(T)*n_1*n_2);
          return;
        }
      for(int i=0; i<n_1*n_2; i++)
        pArray[i] = C;
    }

    /**
    Returns the product of an array and a scalar.
    */
    Array2D operator*(const T C);


    /**
    Sets this array equal to itself times a scalar value.
    */
    void operator*=(const T C);


    /**
    Sets this array equal to itself divided by a scalar value.
    */

    void operator/=(const T C);


    /**
    Creates an array.
    */
    Array2D()
    {
      pArray = 0; n_1 = 0; n_2 = 0;
    }


    /**
    Creates an array and allocates memory.

    @param i size of the array's first dimension.
    @param j size of the array's second dimension.
    */
    Array2D(int i, int j)
    {
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
    This function makes it easy to replace a row with an Array1D
    */
    void setRow(int whichRow, Array1D<T> & rhs);

    /**
    This function takes advantage of the row-major format of the matricies
    to copy numRows worth of data from one Array2D to another
    */
    void setRows(int to, int from, int numRows, const Array2D<T> & rhs);

    /**
    Copies a column from one Array2D to another
    */
    void setColumn(int to, int from, const Array2D<T> & rhs);

    /**
    Destroy's the array and cleans up the memory.
    */
    ~Array2D()
    {
      deallocate();
    }

    /**
    This particular choice indicates row-major format.
    */
    inline int map(int i, int j) const
      {
        return n_2*i + j;
      }

    /**
    Accesses element <code>(i,j)</code> of the array.
    */
    T& operator()(int i,int j)
    {
      assert(pArray != 0);
      return pArray[map(i,j)];
    }

    /**
    Accesses element <code>(i,j)</code> of the array, except this prevents
    modification of that element.
    */
    T get(int i, int j) const
        {
          assert(pArray != 0);
          return pArray[map(i,j)];
        }

    /**
    Prints the array to a stream.
    printf("%11.3e",pix[index +1]);
    printf("   ");
    */
    template <class TT>
    friend ostream& operator<<(ostream & strm, const Array2D<TT> & rhs);
  };

#endif
