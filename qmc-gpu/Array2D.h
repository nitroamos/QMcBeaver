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

using namespace std;


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

	T** pArray;

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
	T** array(){return pArray;}

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
				pArray = new T*[n_1];
				for(int ii=0; ii< n_1; ii++)
				{
					pArray[ii] = new T[n_2];
				}
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
		for(int i=0; i<n_1; i++) delete [] pArray[i];
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

		for(int i=0; i<n_1;i++)
			for(int j=0; j<n_2;j++)
				pArray[i][j] = rhs.pArray[i][j];
	}


	/**
	Sets all of the elements in an array equal to the same value.
	*/

	void operator=(const T C)
	{
		for(int i=0; i<n_1;i++)
			for(int j=0; j<n_2;j++)
				pArray[i][j] = C;
	}


	/**
	Returns the matrix product of two arrays.
	*/

	Array2D operator*(const Array2D & rhs)
	{
		if(n_2 != rhs.n_1)
		{
			cerr << "ERROR: Matrices of incorrect dimensions are being"
				<< " multiplied!" << endl;
			exit(1);
		}
		Array2D<T> TEMP(n_1,rhs.n_2);
		for(int i=0;i<n_1;i++)
		{
			for(int j=0;j<rhs.n_2;j++)
			{
				TEMP.pArray[i][j] = 0;
				for(int k=0;k<n_2;k++)
				{
					TEMP.pArray[i][j] += pArray[i][k] * rhs.pArray[k][j];
				}
			}
		}
		return TEMP;
	}


	/**
	Returns the product of an array and a scalar.
	*/

	Array2D operator*(const T C)
	{
		Array2D<T> TEMP(n_1,n_2);
		for(int i=0;i<n_1;i++)
		{
			for(int j=0;j<n_2;j++)
			{
				TEMP.pArray[i][j] = C*pArray[i][j];
			}
		}
		return TEMP;
	}


	/**
	Sets this array equal to itself times a scalar value.
	*/

	void operator*=(const T C)
	{
		for(int i=0;i<n_1;i++)
		{
			for(int j=0;j<n_2;j++)
			{
				pArray[i][j] *= C;
			}
		}
	}


	/**
	Sets this array equal to itself divided by a scalar value.
	*/

	void operator/=(const T C)
	{
		for(int i=0;i<n_1;i++)
		{
			for(int j=0;j<n_2;j++)
			{
				pArray[i][j] /= C;
			}
		}
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

	Array2D(int i, int j) {pArray = 0; n_1 = 0; n_2 = 0; allocate(i,j);}


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
		for(int i = 0; i<n_1; i++)
			for(int j = 0; j<n_2; j++) pArray[i][j] = rhs.pArray[i][j];
	}


	/**
	Destroy's the array and cleans up the memory.
	*/

	~Array2D(){deallocate();}


	/**
	Accesses element <code>(i,j)</code> of the array.
	*/
	T& operator()(int i,int j){return pArray[i][j];}


	/**
	Prints the array to a stream.
	*/

	friend ostream& operator<<(ostream & strm, const Array2D<T> & rhs)
	{
		for(int i=0; i<rhs.n_1;i++)
		{
			for(int j=0; j<rhs.n_2; j++)
			{
				strm << rhs.pArray[i][j] << "\t";
			}
			strm << endl;
		}
		return strm;
	}

};

#endif
