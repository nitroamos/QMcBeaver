/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#ifndef QMCConfigIO_H
#define QMCConfigIO_H

#include <iostream>
#include <fstream>
#include <math.h>

#include "Array1D.h"
#include "Array2D.h"
#include <string>

using namespace std;

#ifdef USEHDF5
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

typedef struct {
  double SCF_Lap;
  double weight;
  double PE;
  double lnJ;
} ct_typ;

#endif

/**
  This class is meant to handle all of the file I/O for configuration files,
  which can easily get up to several gigabytes for longer runs on larger
  molecules.
  
  It has two modes -- text and binary. I was going to add HDF5, but it looks
  like they're going to be adding a C++ PacketTable class in an upcoming
  release, so I'll wait until then since that seems like the way to go.
  
  There are a couple advantages to binary output. First, a quick measurement
  with HMX reveals that binary files are about 50% the size of an equivalent
  text output. Second, binary output means that the code does not have to convert
  doubles etc into chars, possibly preserving a little more precision as well as
  effort. Third, the text output option uses endl frequently. endl flushes the buffer
  every time, meaning the code has to access the disk quite often. Binary does
  automatically flush -- we can control when.
  
  some notes:
  1) obviously, the binary output cfgs files are only meant to be used by the program
  within the same run. e.g. they are not portable at all.
  
  2) i played around with gzip on the binary output, and the file size difference was
  very small. i conclude from this that compression doesn't make any difference.

  3) Since there can only ever be one file (per node), all the data members are static.
*/
class QMCConfigIO
{
public:
  
  /**
    Constructor. Assigns zero to the fstream pointer.
  */
  QMCConfigIO();

  /**
    Destructor. Calls close().
  */
  ~QMCConfigIO();
  
  /**
    Constructor. Assigns zero to the fstream pointer.
    
    @param since two of the arrays in the file depend on this,
      we might as well keep this number in the class.
  */
  QMCConfigIO(int numElectrons);

  /**
    Returns a string with the name of the file
    this class is writing/reading. If the file
    is opened for writing, it is opened in
    ios_base::trunc mode.
    
    @param filename the name of the file to open
    @param forWriting whether this instance of QMCConfig is
      meant to be reading from or writing to the file.
	*/
  void open(string filename, bool forWriting);

  /**
    Closes the file and then deletes the reference and
    sets the pointer to zero.
  */
  void close();

  /**
    Are there more records in the file?
  */
  bool eof();
    
  /**
     Returns a string with the name of the file
     this class is writing/reading.
     
     @return The name of the file.
  */
  string getFilename();
  
  /**
     @return The number of configurations we wrote
     into the file.
  */
  int getNumberConfigsWritten();

 
  /**
    This will set all the precision requirements we want
    for our output stream.
  */
  void setPrecision();

  /**
     This will write the contents of all the arguments
     sent to this function into the configuration file.
     It will only print numElectrons x 3 elements from
     the arrays.
     
     @param R position vectors for all the electrons
     @param SCF_Laplacian_PsiRatio the calculated "D1"
     @param SCF_Grad_PsiRatio gradient vectors for all the electrons
     @param lnJastrow the result of the Jastrow functions
     @param PE the potential energy for the configuration
     @param weight a weight associated with the configuration
  */
  void writeCorrelatedSamplingConfiguration(
    Array2D<double> & R,
    double SCF_Laplacian_PsiRatio,
    Array2D<double> & SCF_Grad_PsiRatio,
    double lnJastrow,
    double PE,
    double weight);

  /**
     This will read one set of parameters from the configuration
     file and store them in the arguments to this function.
     
     @param R position vectors for all the electrons
     @param SCF_Laplacian_PsiRatio the calculated "D1"
     @param SCF_Grad_PsiRatio gradient vectors for all the electrons
     @param lnJastrow the result of the Jastrow functions
     @param PE the potential energy for the configuration
     @param weight a weight associated with the configuration
  */  
  void readCorrelatedSamplingConfiguration(
    Array2D<double> & R,
    double & SCF_Laplacian_PsiRatio,
    Array2D<double> & SCF_Grad_PsiRatio,
    double & lnJastrow,
    double & PE,
    double & weight);
    
  /**
    Just in case we want to augment our output file...
    Anything piped through this function will be printed
    in text, regardless of whether we're in binary mode.
  */
  template<class T> 
  fstream &  operator<<(T rhs)
  {
    *config_strm << rhs;
    return *config_strm;
  }

private:  
  /**
    This is a pointer we can use to store our the information
    for the stream.
  */
#ifdef USEHDF5
  static H5File * h5_f;

  static DataSet * dst_r;
  static DataSet * dst_g;
  static DataSet * dst_o;
  static CompType * ct;
  ct_typ packet;

#else
  static fstream *config_strm;
#endif

  /**
    The filename we opened the stream with.
  */
  static string filename;
  
  /**
    Are we reading or writing?
  */
  static bool areWriting;
  
  /**
    Number of electrons.
  */  
  static int numElectrons;

  /**
     Number of configurations written
  */
  static int numWritten;
  static int numRead;
};

#endif

