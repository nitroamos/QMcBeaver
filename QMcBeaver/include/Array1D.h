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
#include <iomanip>
#include <math.h>
#include "cppblas.h"
#include <assert.h>

using namespace std;

/**
A 1-dimensional template for making arrays.  All of the memory allocation
and deallocation details are dealt with by the class.

Some of the methods can now be run with ATLAS. Some of the other methods (e.g. rotation)
can also be sent to ATLAS, but these aren't bottlenecks...
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

    int dim1() const
      {
        return n_1;
      }

    /**
    Gets the total number of elements in the array.
    @return total number of elements in the array.
    */

    int size() const
      {
        return n_1;
      }

    /**
    Gets a pointer to an array containing the array elements.  The ordering of
    this array is NOT specified.
    */

    T* array()
    {
      return pArray;
    }

    /**
    Allocates memory for the array.
    @param i size of the array's first dimension.
    */

    void allocate(int i)
    {
#ifdef QMC_DEBUG
      /*
	//There are a number of places in the code that allocate
	//negative numbers, e.g. Polynomial::initialize,
	//i don't want to change them, so i'm
	//commenting out this warning.
      if ( i < 0 )
        {
          cerr << "Warning: invalid dimension in Array1D::allocate\n"
	       << " i = " << i << endl;
        }
      */
#endif
      if ( n_1 != i )
        {
          deallocate();
          n_1 = i;
          if (n_1 >= 1)
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
       This is for compatability with std::vector. It just
       calls allocate(int i).
    */
    void resize(int i)
    {
      allocate(i);
    }

    /**
     This is for compatability with std::vector. It just
     deallocates.
    */
    void clear()
    {
      deallocate();
    }

    /**
       Deallocates memory for the array.
    */

    void deallocate()
    {
      if ( n_1 > 0 )
        {
          delete [] pArray;
          pArray = 0;
          n_1 = 0;
        }
    }

#ifdef USEBLAS
    T operator*( const Array1D<double> & rhs)
    {
#ifdef QMC_DEBUG
      assert(n_1 > 0);
      assert(pArray);
      assert(rhs.pArray);
#endif
#ifdef USE_CBLAS
      return cblas_ddot(n_1, pArray, 1,rhs.pArray, 1);
#else
      int inc = 1;
      return ddot_(&n_1, pArray, &inc, rhs.pArray, &inc);
#endif
    }

    T operator*( const Array1D<float> & rhs)
    {
#ifdef QMC_DEBUG
      assert(n_1 > 0);
      assert(pArray);
      assert(rhs.pArray);
#endif
#ifdef USE_CBLAS
      return cblas_sdot(n_1, pArray, 1,rhs.pArray, 1);
#else
      int inc = 1;
      return sdot_(&n_1, pArray, &inc, rhs.pArray, &inc);
#endif
    }

    //this method has not been tested with ATLAS
    Array1D operator+( const Array1D & rhs)
    {
      Array1D <T> A(rhs);
#ifdef USE_CBLAS
      cblas_daxpy(n_1, 1.0, pArray,1, A.pArray, 1);
#else
      int inc = 1;
      double a = 1.0;
      daxpy_(&n_1, &a, pArray, &inc, A.pArray, &inc);
#endif
      return A;
    }

    //this method has not been tested with ATLAS
    void operator+=( const Array1D & rhs)
    {
#ifdef USE_CBLAS
      cblas_daxpy(n_1, 1.0, rhs.pArray,1, pArray, 1);
#else
      int inc = 1;
      double a = 1.0;
      daxpy_(&n_1, &a, rhs.pArray, &inc, pArray, &inc);
#endif
    }

#else
    /**
    Returns the dot product of two arrays.
    */

    T operator*( const Array1D & rhs)
    {
#ifdef QMC_DEBUG
      if ( n_1 != rhs.n_1 )
        {
          cerr << "ERROR: incorrect sizes for two Array1D's in Array1D*"
          << "Array1D" << endl;
          assert(0);
        }
      if ( n_1 <= 0 )
        {
          cerr << "ERROR: Array1D of size 0 used in Array1D*Array1D" << endl;
          assert(0);
        }
      assert(pArray);
      assert(rhs.pArray);
#endif
      T temp = 0.0;
      for ( int i=0; i<n_1; i++ )
        {
          temp += pArray[i] * rhs.pArray[i];
        }
      return temp;
    }

    /**
    Returns the sum of two arrays.
    */

    Array1D operator+( const Array1D & rhs)
    {
#ifdef QMC_DEBUG
      if ( n_1 != rhs.n_1 )
        {
          cerr << "ERROR: incorrect sizes for two Array1D's in Array1D+"
          << "Array1D" << endl;
          exit(1);
        }
      if ( n_1 <= 0 )
        {
          cerr << "ERROR: Array1D of size 0 used in Array1D+Array1D" << endl;
          exit(1);
        }
      assert(pArray);
      assert(rhs.pArray);
#endif
      Array1D <T> A(n_1);

      for ( int i=0; i<n_1; i++ )
        {
          A.pArray[i] = pArray[i] + rhs.pArray[i];
        }
      return A;
    }

    void operator+=( const Array1D<T> & rhs)
      {
	for(int i=0; i<n_1; i++)
	  pArray[i] += rhs.pArray[i];
      }
    
#endif
    /**
    Sets two arrays equal.
    */

    void operator=(const Array1D & rhs)
    {
      if (n_1 != rhs.n_1) allocate(rhs.n_1);

      for (int i=0; i<n_1;i++) pArray[i] = rhs.pArray[i];
      //memcpy(pArray, rhs.pArray, sizeof(T)*n_1);
    }

    /**
    Sets all of the elements in an array equal to the same value.
    */

    void operator=(const T C)
    {
      if (C == 0)
        {
          memset(pArray,0,sizeof(T)*n_1);
          return;
        }
      for (int i=0; i<n_1; i++)
        pArray[i] = C;
    }

    /**
    Returns the product of an array and a double.
    */

    Array1D operator*( const double rhs)
    {
      Array1D <T> A(n_1);
      for ( int i=0; i<n_1; i++ )
        {
          A.pArray[i] = rhs*pArray[i];
        }
      return A;
    }

    /**
    Sets this array equal to itself times a scalar value.
    */

    void operator*=(const T C)
    {
      for (int i=0;i<n_1;i++)
        {
          pArray[i] *= C;
        }
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
    Returns the difference of two arrays.
    */

    Array1D operator-( const Array1D & rhs)
    {
      if ( n_1 != rhs.n_1 )
        {
          cerr << "ERROR: incorrect sizes for two Array1D's in Array1D-"
          << "Array1D" << endl;
          exit(1);
        }
      if ( n_1 <= 0 )
        {
          cerr << "ERROR: Array1D of size 0 used in Array1D-Array1D" << endl;
          exit(1);
        }
      Array1D <T> A(n_1);

      for ( int i=0; i<n_1; i++ )
        {
          A.pArray[i] = pArray[i] - rhs.pArray[i];
        }
      return A;
    }

    /**
    Conjugates a quaternion.
    */

    void quaternion_conjugate()
    {
      if ( n_1 != 4 )
        {
          cerr << "ERROR: incorrect size for quaternion conjugate." << endl;
          exit(1);
        }
      for (int i=1; i<4; i++)
        {
          pArray[i] *= -1;
        }
    }

    /**
    Returns the product of two quaternions.
    @param rhs quaternion to multiply by this quaternion.
    @return product of these two quaternions.
    */

    Array1D quaternion_product( const Array1D & rhs )
    {
      if ( n_1 != 4 || rhs.n_1 != 4 )
        {
          cerr << "ERROR: incorrect size for quaternion product." << endl;
          exit(1);
        }
      Array1D<T> product(4);
      product.pArray[0] = pArray[0]*rhs.pArray[0] - pArray[1]*rhs.pArray[1] -
                          pArray[2]*rhs.pArray[2] - pArray[3]*rhs.pArray[3];
      product.pArray[1] = pArray[0]*rhs.pArray[1] + pArray[1]*rhs.pArray[0] +
                          pArray[2]*rhs.pArray[3] - pArray[3]*rhs.pArray[2];
      product.pArray[2] = pArray[0]*rhs.pArray[2] - pArray[1]*rhs.pArray[3] +
                          pArray[2]*rhs.pArray[0] + pArray[3]*rhs.pArray[1];
      product.pArray[3] = pArray[0]*rhs.pArray[3] + pArray[1]*rhs.pArray[2] -
                          pArray[2]*rhs.pArray[1] + pArray[3]*rhs.pArray[0];
      return product;
    }

    /**
    Rotates a point by an angle about an axis.
    @param angle
    @param axis 3D coords of axis of unit length about which we are rotating.
    */

    void rotate(Array1D axis, double angle)
    {
      if (n_1 != 3)
        {
          cerr << "ERROR: incorrect size for rotating point." << endl;
          exit(1);
        }
      Array1D<T> q_point(4), q_axis(4), q1(4), q2(4);
      q_point.pArray[0] = 0.0;
      q_axis.pArray[0] = cos(angle/2);
      for (int i=1; i<4; i++)
        {
          q_point.pArray[i] = pArray[i-1];
          q_axis.pArray[i] = axis.pArray[i-1]*sin(angle/2);
        }
      q1 = q_axis.quaternion_product(q_point);
      q_axis.quaternion_conjugate();
      q2 = q1.quaternion_product(q_axis);
      for (int i=1; i<4; i++)
        {
          pArray[i-1] = q2.pArray[i];
        }
    }

    /**
    Creates an array.
    */

    Array1D()
    {
      pArray = 0; n_1 = 0;
    }

    /**
    Creates an array and allocates memory.
    @param i size of the array's first dimension.
    */

    Array1D(int i)
    {
      pArray = 0; n_1 = 0; allocate(i);
    }

  /**
     Create a random array.
     It uses the UNIX rand function, which will
     produce a uniform random number, [0,1)
     
     @param shift will be added to the random value
     @param range after adding shift, we'll multiply by range
     @param iseed will be used to initialize srand, if iseed > 0
  */
  Array1D(int i, double shift, double range, long & iseed)
    {
      pArray = 0;
      n_1 = 0;
      allocate(i);
      
      if(iseed > 0) srand(iseed);
      for(int idx=0; idx<n_1; idx++)
        {
          // between 0 and 1.0
          double val = (double)(rand())/RAND_MAX;
          val = range * (val + shift);
          pArray[idx] = (T)val;
        }
    }
  
    /**
    Creates an array and sets it equal to another array.
    @param rhs array to set this array equal to.
    */

    Array1D(const Array1D & rhs)
    {
      n_1 = 0;
      pArray = 0;
      allocate(rhs.n_1);
      for (int i=0; i<n_1; i++) pArray[i] = rhs.pArray[i];
    }

    /**
    Destroys the array and cleans up the memory.
    */

    ~Array1D()
    {
      deallocate();
    }

    /**
    Accesses element <code>(i)</code> of the array.
    */

    T& operator()(int i)
    {
#ifdef QMC_DEBUG
      if( i < 0 || i >= n_1 )
	{
	  cerr << "Error: Array1D::operator() has i = " << i
	       << " where n_1 = " << n_1 << endl;
	  assert( i >= 0 && i < n_1 );
	}
      assert(pArray);
#endif
      return pArray[i];
    }

    /**
     It's nice to have a const accessor
    */
    T get(int i) const
      {
#ifdef QMC_DEBUG
	assert( i >= 0 && i < n_1 );
	assert(pArray);
#endif
        return pArray[i];
      }

    T& operator[](int i)
    {
#ifdef QMC_DEBUG
      assert( i >= 0 && i < n_1 );
      assert(pArray);
#endif
      return pArray[i];
    }

    /**
    Prints the array to a stream.
    */

    friend ostream& operator<<(ostream & strm, const Array1D<T> & rhs)
    {
      rhs.printArray1D(strm, 15, 5, 8, ',', false);
      return strm;
    }

  /**
     This makes it easier to print out the Array1D
     for input into Mathematica. In Mathematica, you will
     still have to replace all the scientific notational "e" with "*^".
  */
  void printArray1D(ostream & strm, int numSigFig, int columnSpace, int maxI, char sep, bool sci) const
    {
      strm.precision(4);
      if(maxI < 0) maxI = n_1;
      strm << "{";
      for(int i=0; i<n_1;i++)
        {
	  if(i % maxI == 0 && i != 0)
	    strm << endl << " ";

	  if(sci) strm.setf(ios::scientific);
	  else strm.unsetf(ios::scientific);

	  strm << setprecision(numSigFig)
	       << setw(numSigFig+columnSpace)
	       << pArray[i];
	  if(i<n_1-1) strm << sep;
        }
      strm << " }";
    }

    /**
       This function assumes that the Array1D was initialized to the expected
       number of parameters.

       This function will read into Array1D a set of parameters, separated
       by whitespace. It will read in as many numbers as it can until
       it runs into a letter [a-zA-Z]. At that point, if too many parameters
       were entered, it will print a warning, using 'label' to help the user
       find the problem.

       If too few parameters were entered, then it will initialize the rest to
       the provided initialization value.

       @param strm the input stream where the parameters are ready to be read
       @param initialization the value to use to fill in missing entries
       @param label printed with the warning for loading too many variables
       @return the number of parameters that were available
    */
    int read(istream & strm, T initialization, string label)
      {
	  int pi = 0;
	  strm >> ws;
	  string temp;
	  while(!isalpha(strm.peek()))
	    {
	      strm >> temp;
	      if(pi < n_1)
		pArray[pi] = (T)(atof(temp.c_str()));
	      pi++;
	      strm >> ws;
	    }

	  if(pi > n_1)
	    clog << "Warning: you entered " << pi << " parameters, when only "
		 << n_1 << " were expected for " << label << endl;  

	  for (int i=pi; i < n_1; i++)
	    pArray[i] = initialization;

	  return pi;
      }
  };

#endif
