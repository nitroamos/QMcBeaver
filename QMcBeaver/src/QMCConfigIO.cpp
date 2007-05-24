/*
  Copyright (c) Amos G. Anderson 2005
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "QMCConfigIO.h"

/**
   You want to be writting in binary, unless you're debugging.
   The inBinary parameter is ignored if HDF5 is enabled.

   Testing on matrix at least, HDF5 is slower than writing in binary,
   when the speeds maybe should be comparable.
   So, more tuning is probably appropriate on a system by system
   basis.

   Each configuration writes/reads 4 doubles and 2 arrays.
   The big question is whether these should go into one
   dataset, or two datasets (the arrays are the same dimensions),
   or 3 datasets.

   This does not use the high level HDF5 API.
*/
static const bool inBinary = true;

int QMCConfigIO::numRead      = 0;
int QMCConfigIO::numWritten   = 0;
int QMCConfigIO::numElectrons = 0;
bool QMCConfigIO::areWriting  = false;
string QMCConfigIO::filename  = "";

#ifdef USEHDF5
H5File * QMCConfigIO::h5_f      = 0;
DataSet * QMCConfigIO::dst_r    = 0;
DataSet * QMCConfigIO::dst_g    = 0;
DataSet * QMCConfigIO::dst_o    = 0;
CompType * QMCConfigIO::ct      = 0;
#else
fstream * QMCConfigIO::config_strm = 0;
#endif

QMCConfigIO::QMCConfigIO()
{
#ifdef USEHDF5
  h5_f  = 0;
  dst_r = 0;
  dst_o = 0;
  dst_g = 0;
  ct = 0;
#else
  config_strm = 0;
#endif
  areWriting = true;
  numWritten = 0;
  numRead    = 0;
  filename = "";
}

QMCConfigIO::~QMCConfigIO()
{
  close();
}

QMCConfigIO::QMCConfigIO(int numE)
{
#ifdef USEHDF5
  h5_f  = 0;
  dst_r = 0;
  dst_o = 0;  
  dst_g = 0;  
  ct = 0;
#else
  if(!inBinary)
    {
      clog << "Warning: CFGS file is not binary, severe performance penalty!!" << endl;
    }

  config_strm = 0;
#endif
  areWriting = true;
  numWritten = 0;
  numRead    = 0;
  filename = "";
  numElectrons = numE;
}

void QMCConfigIO::open(string name, bool forWriting)
{
#ifdef USEHDF5
  //I think it's a good idea to use the common HDF5 file extension
  name += ".h5";
  if( name == filename && areWriting == forWriting && h5_f != 0 )
    return;
#else
  if( name == filename && areWriting == forWriting && config_strm != 0 )
    return;
#endif

  areWriting = forWriting;
  filename = name;
  close();

  if(forWriting)
    numWritten = 0;

  if(filename == "")
  {
    clog << "Warning: empty filename" << endl;
    return;
  }
  
#ifdef USEHDF5

  if( !forWriting && !H5File::isHdf5(filename) )
    {
      clog << "Error: the file " << filename << " is not an HDF5 file!" << endl;
    }

  bool printError = false;
  try {
    //Exception::dontPrint();

    ct = new CompType(sizeof(ct_typ));
    ct->insertMember("SCF_Lap",  HOFFSET(ct_typ, SCF_Lap ), PredType::NATIVE_DOUBLE);
    ct->insertMember("weight",   HOFFSET(ct_typ, weight  ), PredType::NATIVE_DOUBLE);
    ct->insertMember("PE",       HOFFSET(ct_typ, PE      ), PredType::NATIVE_DOUBLE);
    ct->insertMember("lnJ",      HOFFSET(ct_typ, lnJ     ), PredType::NATIVE_DOUBLE);

    const hsize_t dim = 3*numElectrons;

    /*
      See User Guide, Chapter 2, section 8 for descriptions
      the "core" driver keeps it all in RAM, paramter is how much memory to increment
      the "sec2" driver is default, minimal buffering      
    */
    FileAccPropList   access;
    //access.setCore(1024, (hbool_t)false);
    //access.setStdio();

    FileCreatPropList create;

    if(forWriting)
      {
	h5_f = new H5File(filename, H5F_ACC_TRUNC,  create.getId(), access.getId());
	DSetCreatPropList dcpl_r, dcpl_o, dcpl_g;

	/*
	  Higher number corresponds to more compression.
	  0 through 9 allowed
	*/
	const int compressionLevel = 0;
	
	/*
	  The larger this parameter is, the more HDF5 will buffer
	  before writting. Testing on matrix suggests this only
	  makes a couple % of difference.
	*/
	const hsize_t chunkFactor = 100;


	//what's the best way to order the dimensions?
	//which one should be unlimited?
	hsize_t max_rgr_size[2] = {H5S_UNLIMITED, dim};
	hsize_t     rgr_size[2] = {chunkFactor,   dim};
	DataSpace dsp1(2,rgr_size,max_rgr_size);
	dcpl_r.setChunk(2,rgr_size);
	dcpl_r.setDeflate(compressionLevel);
	dcpl_r.setFillTime(H5D_FILL_TIME_NEVER);
	dcpl_r.setAllocTime(H5D_ALLOC_TIME_LATE);
	dst_r = new DataSet(h5_f->createDataSet("R",PredType::NATIVE_DOUBLE,dsp1,dcpl_r));
	dsp1.close();

	DataSpace dsp3(2,rgr_size,max_rgr_size);
	dcpl_g.setChunk(2,rgr_size);
	dcpl_g.setDeflate(compressionLevel);
	dcpl_g.setFillTime(H5D_FILL_TIME_NEVER);
	dcpl_g.setAllocTime(H5D_ALLOC_TIME_LATE);
	dst_g = new DataSet(h5_f->createDataSet("Grad",PredType::NATIVE_DOUBLE,dsp3,dcpl_g));
	dsp3.close();
	
	hsize_t max_oth_size[1] = {H5S_UNLIMITED};
	hsize_t     oth_size[1] = {chunkFactor};
	DataSpace dsp2(1,oth_size,max_oth_size);
	dcpl_o.setChunk(1,oth_size);
	dcpl_o.setDeflate(compressionLevel);
	dcpl_o.setFillTime(H5D_FILL_TIME_NEVER);
	dcpl_o.setAllocTime(H5D_ALLOC_TIME_LATE);
	dst_o = new DataSet(h5_f->createDataSet("Others",*ct,dsp2,dcpl_o)); 
	dsp2.close();

	numWritten = 0;

      } else {
	h5_f = new H5File(filename, H5F_ACC_RDONLY, create.getId(), access.getId());
	dst_r = new DataSet(h5_f->openDataSet("R")); 
	dst_g = new DataSet(h5_f->openDataSet("Grad")); 
	dst_o = new DataSet(h5_f->openDataSet("Others")); 
	numRead = 0;
      }
  }

  // catch failure caused by the H5File operations
  catch( FileIException error )
    {
      error.printError();
      printError = true;
    }
  
  // catch failure caused by the DataSet operations
  catch( DataSetIException error )
    {
      error.printError();
      printError = true;
    }
  
  // catch failure caused by the DataSpace operations
  catch( DataSpaceIException error )
    {
      error.printError();
      printError = true;
    }
  
  // catch failure caused by the DataSpace operations
  catch( DataTypeIException error )
    {
      error.printError();
      printError = true;
    }

  if(printError)
    {
      cout << "Error: HDF5 can't open file " << filename << endl;
    }

#else  
  ios_base::openmode mode;
  
  if(forWriting)
    {
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
#endif
}

void QMCConfigIO::writeCorrelatedSamplingConfiguration(Array2D<double> & R,
                                                       double SCF_Laplacian_PsiRatio,
                                                       Array2D<double> & SCF_Grad_PsiRatio,
                                                       double lnJastrow,
                                                       double PE,
						       double weight)
{
  if(!areWriting)
    {
      cerr << "Warning: the file is open for reading, not writing." << endl;
      return;
    }

#ifdef USEHDF5

  packet.SCF_Lap = SCF_Laplacian_PsiRatio;
  packet.weight  = weight;
  packet.PE      = PE;
  packet.lnJ     = lnJastrow;
  
  hsize_t      dim       = 3*numElectrons;
  hsize_t      dims1[2]  = {1,              dim};  
  hsize_t      size1[2]  = {numWritten+100, dim};
  hsize_t      offs1[2]  = {numWritten,       0};

  dst_r->extend( size1 );
  DataSpace ms1(2, dims1);
  DataSpace fs1,fs3;
  fs1 = dst_r->getSpace();
  fs1.selectHyperslab(H5S_SELECT_SET, dims1, offs1);
  dst_r->write(R.array(), PredType::NATIVE_DOUBLE, ms1, fs1);

  dst_g->extend( size1 );
  DataSpace ms3(2, dims1);
  fs3 = dst_g->getSpace();
  fs3.selectHyperslab(H5S_SELECT_SET, dims1, offs1);
  dst_g->write(SCF_Grad_PsiRatio.array(), PredType::NATIVE_DOUBLE, ms3, fs3);
  fs3.close();
  ms3.close();

  hsize_t      dims2[1]  = {1};  
  hsize_t      size2[1]  = {numWritten+100};
  hsize_t      offs2[1]  = {numWritten};
  dst_o->extend( size2 );

  DataSpace ms2(1, dims2);
  DataSpace fs2;
  fs2 = dst_o->getSpace();
  fs2.selectHyperslab(H5S_SELECT_SET, dims2, offs2);
  dst_o->write(&packet, (*ct), ms2, fs2);
  fs2.close();
  ms2.close();
  
#else
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (write sample)." << endl;
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
#endif

  numWritten++;
}

void QMCConfigIO::readCorrelatedSamplingConfiguration(Array2D<double> & R,
                                                      double & SCF_Laplacian_PsiRatio,
                                                      Array2D<double> & SCF_Grad_PsiRatio,
                                                      double & lnJastrow,
                                                      double & PE,
						      double & weight)
{
  if(areWriting)
    {
      cerr << "Warning: the file is open for writing, not reading." << endl;
      return;
    }

#ifdef USEHDF5
  int err = 0;

  hsize_t      dim       = 3*numElectrons;
  hsize_t      dims1[2]  = {1,          dim};  
  hsize_t      offs1[2]  = {numRead,      0};
  offs1[0] = numRead;
  
  DataSpace ms1(2, dims1);
  DataSpace fs1;
  fs1 = dst_r->getSpace();
  fs1.selectHyperslab(H5S_SELECT_SET, dims1, offs1);
  dst_r->read(R.array(), PredType::NATIVE_DOUBLE, ms1, fs1);
  fs1.close();

  DataSpace ms3(2, dims1);
  DataSpace fs3;
  fs3 = dst_g->getSpace();
  fs3.selectHyperslab(H5S_SELECT_SET, dims1, offs1);
  dst_g->read(SCF_Grad_PsiRatio.array(), PredType::NATIVE_DOUBLE, ms3, fs3);
  fs3.close();
  
  hsize_t      dims2[1]  = {1};  
  hsize_t      offs2[1]  = {numRead};
  offs2[0] = numRead;
  
  DataSpace ms2(1, dims2);
  DataSpace fs2;
  fs2 = dst_o->getSpace();
  fs2.selectHyperslab(H5S_SELECT_SET, dims2, offs2);  
  dst_o->read(&packet, (*ct), ms2, fs2);
  fs2.close();
  
  if(err < 0) exit(0);
  SCF_Laplacian_PsiRatio = packet.SCF_Lap;
  weight                 = packet.weight;
  PE                     = packet.PE;
  lnJastrow              = packet.lnJ;
  
#else
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (read sample)." << endl;
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
#endif

  numRead++;
}

void QMCConfigIO::close()
{
#ifdef USEHDF5
  if(h5_f != 0){
    h5_f->close();
    delete h5_f;     h5_f  = 0;
    delete dst_r;     dst_r  = 0;
    delete dst_o;     dst_o  = 0;
    delete dst_g;     dst_g  = 0;
    delete ct;       ct    = 0;
  }
#else
  if(config_strm != 0){
    config_strm->close();
    delete config_strm;
    config_strm = 0;
  }
#endif
}

string QMCConfigIO::getFilename()
{
  return filename;
}

bool QMCConfigIO::eof()
{
#ifdef USEHDF5
  int err;

  if(numRead < numWritten)
    return false;
  return true;

#else
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (eof)." << endl;
    return true;    
  }
    
  return config_strm->eof();
#endif
}

void QMCConfigIO::setPrecision()
{
#ifdef USEHDF5
#else
  if(config_strm == 0)
  {
    cerr << "Warning: config_strm is zero (setPrecision)." << endl;
    return;    
  }
  
  config_strm->setf(ios::fixed,ios::floatfield);
  config_strm->precision(10);
#endif
}

