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

#include "QMCMolecule.h"

QMCMolecule::QMCMolecule()
{
  Natoms = 0;
}

void QMCMolecule::initialize(int nAtoms)
{
  Natoms = nAtoms;
}

int QMCMolecule::getNumberAtoms()
{
  return Natoms;
}


QMCMolecule QMCMolecule::operator=( const QMCMolecule & rhs )
{
  Natoms         = rhs.Natoms;
  Atom_Labels    = rhs.Atom_Labels;
  Atom_Positions = rhs.Atom_Positions;
  Z              = rhs.Z;
  return *this;
}

int QMCMolecule::findClosestNucleusIndex(Array2D<double> & x, int index)
{
  // Find the closest nucleus to x

  int closest_nucleus_index               = 0;
  double closest_nucleus_distance_squared = 1e100;
  
  for(int nucleus=0; nucleus<getNumberAtoms(); nucleus++)
    {
      double r2 = 0.0;
      
      for(int xyz=0; xyz<3; xyz++)
	{
	  double temp = x(index,xyz) - Atom_Positions(nucleus,xyz);
	  
	  r2 += temp * temp;
	}
      
      if( r2 < closest_nucleus_distance_squared )
	{
	  closest_nucleus_distance_squared = r2;
	  closest_nucleus_index            = nucleus;
	}
    }

  return closest_nucleus_index;
}

istream& operator >>(istream &strm, QMCMolecule &rhs)
{
  rhs.Atom_Labels.allocate(rhs.Natoms);
  rhs.Atom_Positions.allocate(rhs.Natoms,3);
  rhs.Z.allocate(rhs.Natoms);

  for(int i=0; i<rhs.Natoms; i++)
    {
      strm >> rhs.Atom_Labels(i);
      rhs.Atom_Labels(i) = 
	StringManipulation::toFirstUpperRestLower(rhs.Atom_Labels(i));

      if(rhs.Atom_Labels(i) == "&"){
	cerr << "ERROR: End of geometry section reached prematurely. We were expecting " << rhs.Natoms
	     << " atoms, but we only got " << i << ".\n";
	exit(1);
      }

      strm >> rhs.Z(i);
      strm >> rhs.Atom_Positions(i,0);
      strm >> rhs.Atom_Positions(i,1);
      strm >> rhs.Atom_Positions(i,2);
    }
  return strm;
}

void QMCMolecule::read(string runfile)
{
  if(Natoms == 0) return;

  ifstream input_file(runfile.c_str());

  if(!input_file)
    {
      cerr << "ERROR: Could not open " << runfile << "!" << endl;
      exit(1);
    }

  string temp_string;
  input_file >> temp_string;
  while((temp_string != "&geometry") && (input_file.eof() != 1))
    {
      input_file >> temp_string;
    }

  if(temp_string == "&geometry") input_file >> *this;

  // Generate the list of nuclei types
  list<string> nuclei;
  nuclei.push_back( Atom_Labels(0) );
  for( int i=0; i<Atom_Labels.dim1(); i++ )
  {
    bool found = false;
    for( list<string>::iterator it=nuclei.begin(); it!=nuclei.end(); ++it )
    {
      if( *it == Atom_Labels(i) )
      {
        found = true;
      }
    }

    if( !found )
    {
      nuclei.push_back( Atom_Labels(i) );
    }
  }

  NucleiTypes.allocate( (int)nuclei.size() );
  int Position = 0;
  for( list<string>::iterator it=nuclei.begin(); it!=nuclei.end(); ++it )
  {
    NucleiTypes(Position) = *it;
    Position++;
  }
}

ostream& operator <<(ostream& strm, QMCMolecule& rhs)
{
 strm << "&geometry" << endl;
 for(int i=0; i<rhs.Natoms; i++)
 {
  strm << rhs.Atom_Labels(i) << "\t";
  strm << rhs.Z(i) << "\t";
  for(int j=0; j<3; j++) strm << rhs.Atom_Positions(i,j) << "\t";
  strm << endl;
 }
 strm << "&" << endl;
 return strm;
}



