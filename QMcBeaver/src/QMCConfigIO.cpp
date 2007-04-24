/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "QMCConfigIO.h"

/**
  If you want the configuration file to be in binary,
  set this to true. Binary files are much smaller,
  meaning that they take less time to write/read.

  Perhaps eventually I'll add an HDF5 option.
*/
static const bool inBinary = true;

QMCConfigIO::QMCConfigIO()
{
  if(!inBinary)
    {
      clog << "Warning: CFGS file is not being written in binary!!" << endl;
    }
  config_strm = 0;
  areWriting = true;
  numWritten = 0;
  filename = "";
}

QMCConfigIO::~QMCConfigIO()
{
  close();
}

QMCConfigIO::QMCConfigIO(int numE)
{
  if(!inBinary)
    {
      clog << "Warning: CFGS file is not being written in binary!!" << endl;
    }
  config_strm = 0;
  areWriting = true;
  numWritten = 0;
  filename = "";
  numElectrons = numE;
}

void QMCConfigIO::open(string name, bool forWriting)
{
  if( name == filename && areWriting == forWriting && config_strm != 0 )
    return;

  areWriting = forWriting;
  filename = name;
  close();
    
  if(filename == "")
  {
    cerr << "Warning: empty filename" << endl;
    return;
  }
  
  ios_base::openmode mode;
  
  if(forWriting)
  {
    numWritten = 0;
    mode = ios_base::trunc;
    mode |= ios_base::out;
  }
  else 
  {
    mode = ios_base::in;
  }
  
  if(inBinary) mode |= ios_base::binary;
  
  config_strm = new fstream(filename.c_str(), mode);
  
  if( !(*config_strm) || config_strm->bad() )
  {
    cerr << "Error: Can't open " << filename;
    if(config_strm->rdstate() & ios::badbit)
      cerr << " badbit ";
    if(config_strm->rdstate() & ios::failbit)
      cerr << " failbit ";
    if(config_strm->rdstate() & ios::eofbit)
      cerr << " eofbit ";
    cerr << endl;
    exit(1);
  }

  //this is just to make sure that the EOF flag is cleared
  config_strm->clear();
}

void QMCConfigIO::writeCorrelatedSamplingConfiguration(Array2D<double> & R,
                                                       double SCF_Laplacian_PsiRatio,
                                                       Array2D<double> & SCF_Grad_PsiRatio,
                                                       double lnJastrow,
                                                       double PE,
						       double weight)
{

  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (write sample)." << endl;
    return;    
  }
  
  if(!areWriting)
  {
    cerr << "Warning: the file is open for reading, not writing." << endl;
    return;
  }
  
  if(inBinary)
  {
    config_strm->write( (char*) &numElectrons, sizeof(int) );
    config_strm->write( (char*) R.array(), sizeof(*R.array())*R.dim1()*3 );
    config_strm->write( (char*) &SCF_Laplacian_PsiRatio, sizeof(double) );
    config_strm->write( (char*) SCF_Grad_PsiRatio.array(),
                        sizeof(*SCF_Grad_PsiRatio.array())*SCF_Grad_PsiRatio.dim1()*3 );
    config_strm->write( (char*) &lnJastrow, sizeof(double) );
    config_strm->write( (char*) &PE, sizeof(double) ); 
    config_strm->write( (char*) &weight, sizeof(double) ); 
  }
  else
  {
    setPrecision();
    fstream & strm = *config_strm;
    
    strm << "&" << endl;
    strm << "R\t" << numElectrons << endl;
    for(int i=0;i<numElectrons;i++)
    {
      //now we are printing out all the electrons
      
      strm << "\t" << i << "\t";
      strm << R(i,0) << "\t";
      strm << R(i,1) << "\t";
      strm << R(i,2) << endl;
    }
    
    // Print the sum of the laplacians from the alpha and beta electrons
    strm << "D1" << endl;
    strm << "\t" << SCF_Laplacian_PsiRatio << endl;
    
    // print the GradPsiRatio excluding the Jastrow
    strm << "D2" << endl;
    
    for(int i=0; i<numElectrons; i++)
    {
      for(int j=0; j<3; j++)
      {
        strm << "\t" << SCF_Grad_PsiRatio(i,j);
      }
      strm << endl;
    }
    
    strm << "lnJ\t" << endl;
    strm << "\t" << lnJastrow << endl;
    strm << "PE\t" << endl;
    strm << "\t" << PE << endl;
    strm << "W\t" << endl;
    strm << "\t" << weight << endl;
  }

  numWritten++;
}

void QMCConfigIO::readCorrelatedSamplingConfiguration(Array2D<double> & R,
                                                      double & SCF_Laplacian_PsiRatio,
                                                      Array2D<double> & SCF_Grad_PsiRatio,
                                                      double & lnJastrow,
                                                      double & PE,
						      double & weight)
{
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (read sample)." << endl;
    return;    
  }
  
  if(areWriting)
  {
    cerr << "Warning: the file is open for writing, not reading." << endl;
    return;
  }
  
  if(inBinary)
  {
    config_strm->read( (char*) &numElectrons, sizeof(int) );
    config_strm->read( (char*) R.array(), sizeof(*R.array())*R.dim1()*3 );
    config_strm->read( (char*) &SCF_Laplacian_PsiRatio, sizeof(double) );
    config_strm->read( (char*) SCF_Grad_PsiRatio.array(),
                        sizeof(*SCF_Grad_PsiRatio.array())*SCF_Grad_PsiRatio.dim1()*3 );
    config_strm->read( (char*) &lnJastrow, sizeof(double) );
    config_strm->read( (char*) &PE, sizeof(double) ); 
    config_strm->read( (char*) &weight, sizeof(double) ); 
  }
  else
  {
    string stemp;
    int itemp;
    fstream & strm = *config_strm;
    
    strm >> stemp;
    strm >> stemp;
    strm >> stemp;
    
    // Read in the electronic positions
    for(int i=0; i<numElectrons; i++)
    {
      strm >> itemp;

      strm >> R(i,0);
      strm >> R(i,1);
      strm >> R(i,2);
    }
    
    // read D1
    strm >> stemp;
    strm >> SCF_Laplacian_PsiRatio;
    
    // read D2
    strm >> stemp;
    for(int i=0; i<numElectrons; i++)
    {
      strm >> SCF_Grad_PsiRatio(i,0);
      strm >> SCF_Grad_PsiRatio(i,1);
      strm >> SCF_Grad_PsiRatio(i,2);
    }
    
    // read lnJ
    strm >> stemp;
    strm >> lnJastrow;

    // read PE
    strm >> stemp;
    strm >> PE;

    // read the weight
    strm >> stemp;
    strm >> weight;
  }
  config_strm->flush();
}

void QMCConfigIO::close()
{
  if(config_strm != 0){
    config_strm->close();
    delete config_strm;
    config_strm = 0;
  }
}

string QMCConfigIO::getFilename()
{
  return filename;
}

bool QMCConfigIO::eof()
{
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (eof)." << endl;
    return true;    
  }
    
  return config_strm->eof();
}

void QMCConfigIO::setPrecision()
{
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (setPrecision)." << endl;
    return;    
  }
  
  config_strm->setf(ios::fixed,ios::floatfield);
  config_strm->precision(10);
}

