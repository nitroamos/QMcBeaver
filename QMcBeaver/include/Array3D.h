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


#ifndef Array3D_H
#define Array3D_H


/**
  A 3-dimensional template for making arrays.  All of the memory allocation
  and deallocation details are dealt with by the class.
*/


template <class T> class Array3D
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
     Number of elements in the array's third dimension.
  */

  int n_3;


  /**
     Array containing the data.
  */

  T*** pArray;

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
    Gets the number of elements in the array's third dimension.

    @return number of elements in the array's third dimension.
    */
  int dim3(){return n_3;}

  /**
    Gets the total number of elements in the array.

    @return total number of elements in the array.
    */
  int size(){return n_1*n_2*n_3;}

  /**
    Gets a pointer to an array containing the array elements.  
    The ordering of this array is NOT specified.  
    */
  T*** array(){return pArray;}

  /**
     Allocates memory for the array.

     @param i size of the array's first dimension.
     @param j size of the array's second dimension.
     @param k size of the array's third dimension.
  */

  void allocate(int i, int j, int k)
    {
      if( n_1 != i || n_2 != j || n_3 != k )
	{
	  deallocate();

	  n_1 = i;
	  n_2 = j;
	  n_3 = k;

	  if(n_1 >= 1 && n_2 >= 1 && n_3 >= 1)
	    {
	      pArray = new T**[n_1];
	      for(int ii=0; ii< n_1; ii++)
		{
		  pArray[ii] = new T*[n_2];
		  for(int jj=0; jj<n_2; jj++) pArray[ii][jj] = new T[n_3];
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
     for(int i=0; i<n_1; i++)
	  {
	   for(int j=0; j<n_2; j++) delete [] pArray[i][j];
	   delete [] pArray[i];
     }
     delete [] pArray;
     pArray = 0;

     n_1 = 0;
     n_2 = 0;
     n_3 = 0;
    }


  /**
     Sets two arrays equal.
  */

  void operator=(const Array3D & rhs)
    {
     if(n_1 != rhs.n_1 || n_2 != rhs.n_2 || n_3 != rhs.n_3)
            allocate(rhs.n_1,rhs.n_2,rhs.n_3);

      for(int i=0; i<n_1;i++)
        for(int j=0; j<n_2;j++)
          for(int k=0; k<n_3;k++)
            pArray[i][j][k] = rhs.pArray[i][j][k];
    }


  /**
     Creates an array.
  */

  Array3D(){pArray = 0; n_1 = 0; n_2 = 0; n_3 = 0;}


  /**
     Creates an array and allocates memory.
     
     @param i size of the array's first dimension.
     @param j size of the array's second dimension.
     @param k size of the array's third dimension.
  */

  Array3D(int i, int j, int k)
    {pArray = 0; n_1 = 0; n_2 = 0; n_3 = 0; allocate(i,j, k);}


  /**
     Creates an array and sets it equal to another array.

     @param rhs array to set this array equal to.
  */

  Array3D( const Array3D<T> & rhs)
    {
      n_1 = 0;
      n_2 = 0;
      n_3 = 0;
      pArray = 0;
      allocate(rhs.n_1, rhs.n_2, rhs.n_3);

      for(int i=0; i<n_1; i++)
        for(int j=0; j<n_2; j++)
	  for(int k=0; k<n_3; k++)
	    pArray[i][j][k] = rhs.pArray[i][j][k];
    }


  /**
     Destroy's the array and cleans up the memory.
  */

  ~Array3D(){deallocate();}

 /**
    Accesses element <code>(i,j,k)</code> of the array.
    */
  T& operator()(int i,int j, int k){return pArray[i][j][k];}

};

#endif
