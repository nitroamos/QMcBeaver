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

#ifndef Array1D_H
#define Array1D_H

#include <iostream>

using namespace std;


/**
  A 1-dimensional template for making arrays.  All of the memory allocation
  and deallocation details are dealt with by the class.
*/

template <class T> class Array1D
{
 private:
  /**
     Number of elements in the array's first dimension.
  */

  int n_1;


  /**
     Array containing the data.
  */

  T* pArray;

public:
  /**
    Gets the number of elements in the array's first dimension.

    @return number of elements in the array's first dimension.
    */
  int dim1(){return n_1;}

  /**
    Gets the total number of elements in the array.

    @return total number of elements in the array.
    */
  int size(){return n_1;}

  /**
    Gets a pointer to an array containing the array elements.  The ordering of this
    array is NOT specified.  
    */
  T* array(){return pArray;}

  /**
     Allocates memory for the array.

     @param i size of the array's first dimension.
  */

  void allocate(int i)
    {
      if( n_1 != i )
	{
	  deallocate();
	  
	  n_1 = i;
	  
	  if(n_1 >= 1) 
	    {
	      pArray = new T[n_1];
	    }
	  else 
	    {
	      n_1 = 0;
	      pArray = 0;
	    }
	}
    }


  /**
     Deallocates memory for the array.
  */

  void deallocate()
    {
      if( n_1 > 0 )
	{
	  delete [] pArray;
	  pArray = 0;

	  n_1 = 0;
	}
    }


  /**
     Sets two arrays equal.
  */

  void operator=(const Array1D & rhs)
  {
   if(n_1 != rhs.n_1) allocate(rhs.n_1);

   for(int i=0; i<n_1;i++) pArray[i] = rhs.pArray[i];
  }


  /**
     Sets all of the elements in an array equal to the same value.
  */

  void operator=(const T C)
  {
    for(int i=0; i<n_1; i++)
      {
	pArray[i] = C;
      }
  }


  /**
     Returns the dot product of two arrays.
  */

  T operator*( const Array1D & rhs)
    {
      if( n_1 != rhs.n_1 )
	{
	  cerr << "ERROR: incorrect sizes for two Array1D's in Array1D*"
	       << "Array1D" << endl;
	  exit(1);
	}
      if( n_1 <= 0 )
	{
	  cerr << "ERROR: Array1D of size 0 used in Array1D*Array1D" << endl;
	  exit(1);
	}

      double temp = 0.0;
      for( int i=0; i<n_1; i++ )
	{
	  temp += pArray[i] * rhs.pArray[i];
	}
      return temp;
    }


  /**
     Returns the product of an array and a double.
  */

  Array1D operator*( const double rhs)
    {
      Array1D <T> A(n_1);
      for( int i=0; i<n_1; i++ )
	{
	  A.pArray[i] = rhs*pArray[i];
	}
      return A;
    }


  /**
     Returns the sum of two arrays.
  */

  Array1D operator+( const Array1D & rhs)
    {
      if( n_1 != rhs.n_1 )
	{
	  cerr << "ERROR: incorrect sizes for two Array1D's in Array1D+"
	       << "Array1D" << endl;
	  exit(1);
	}
      if( n_1 <= 0 )
	{
	  cerr << "ERROR: Array1D of size 0 used in Array1D+Array1D" << endl;
	  exit(1);
	}
      Array1D <T> A(n_1);

      for( int i=0; i<n_1; i++ )
	{
	  A.pArray[i] = pArray[i] + rhs.pArray[i];
	}
      return A;
    }


  /**
     Returns the difference of two arrays.
  */

  Array1D operator-( const Array1D & rhs)
    {
      if( n_1 != rhs.n_1 )
	{
	  cerr << "ERROR: incorrect sizes for two Array1D's in Array1D-"
	       << "Array1D" << endl;
	  exit(1);
	}
      if( n_1 <= 0 )
	{
	  cerr << "ERROR: Array1D of size 0 used in Array1D-Array1D" << endl;
	  exit(1);
	}
      Array1D <T> A(n_1);

      for( int i=0; i<n_1; i++ )
	{
	  A.pArray[i] = pArray[i] - rhs.pArray[i];
	}
      return A;
    }


  /**
     Sets this array equal to itself times a scalar value.
  */

  void operator*=(const T C)
    {
      for(int i=0;i<n_1;i++)
        {
	  pArray[i] *= C;
        }
    }


  /**
     Sets this array equal to itself divided by a scalar value.
  */

  void operator/=(const T C)
    {
      for(int i=0;i<n_1;i++)
        {
	  pArray[i] /= C;
        }
    }


  /**
     Creates an array.
  */

  Array1D(){pArray = 0; n_1 = 0;}


  /**
     Creates an array and allocates memory.
     
     @param i size of the array's first dimension.
  */

  Array1D(int i){pArray = 0; n_1 = 0; allocate(i);}


  /**
     Creates an array and sets it equal to another array.

     @param rhs array to set this array equal to.
  */

  Array1D( const Array1D & rhs)
    {
      n_1 = 0;
      pArray = 0;
      allocate(rhs.n_1);
      for(int i=0; i<n_1; i++) pArray[i] = rhs.pArray[i];
    }


  /**
     Destroy's the array and cleans up the memory.
  */

  ~Array1D(){deallocate();}

  /**
    Accesses element <code>(i)</code> of the array.
    */
  T& operator()(int i){return pArray[i];}

  /**
     Prints the array to a stream.
  */

  friend ostream& operator<<(ostream & strm, const Array1D<T> & rhs)
    {
      for( int i=0; i<rhs.n_1; i++ )  
	{
	  strm << rhs.pArray[i] << "\t";
	}
      strm << endl;
      return strm;
    }
};

#endif






