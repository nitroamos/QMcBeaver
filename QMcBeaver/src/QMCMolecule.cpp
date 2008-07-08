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
#include "Lebedev-Laikov.h"
#include "MathFunctions.h"
#include <iomanip>

QMCMolecule::QMCMolecule()
{
  Natoms = 0;
}

void QMCMolecule::initialize(int nAtoms, int GridLevel)
{
  Natoms    = nAtoms;
  gridLevel = GridLevel;
}

int QMCMolecule::getNumberAtoms()
{
  return Natoms;
}

int QMCMolecule::getNuclearCharge()
{
  int sum = 0;
  for(int i=0; i<Natoms; i++)
    {
      sum += Zeff(i);
    }
  return sum;
}

QMCMolecule QMCMolecule::operator=( const QMCMolecule & rhs )
{
  Natoms         = rhs.Natoms;
  Atom_Labels    = rhs.Atom_Labels;
  Atom_Positions = rhs.Atom_Positions;
  Z              = rhs.Z;

  gridLevel      = rhs.gridLevel;
  Zeff           = rhs.Zeff;
  usesPseudo     = rhs.usesPseudo;
  Vlocal         = rhs.Vlocal;
  Vnonlocal      = rhs.Vnonlocal;
  grid           = rhs.grid;
  gridWeights    = rhs.gridWeights;
  gridLegendre   = rhs.gridLegendre;

  return *this;
}

int QMCMolecule::findClosestNucleusIndex(Array2D<double> & riA, int e)
{
  // Find the closest nucleus to x
  int closest_nucleus_index       = 0;
  double closest_nucleus_distance = 1e100;
  
  for(int nucleus=0; nucleus<getNumberAtoms(); nucleus++)
    {      
      if( riA(e,nucleus) < closest_nucleus_distance )
	{
	  closest_nucleus_distance = riA(e,nucleus);
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
  rhs.usesPseudo.allocate(rhs.Natoms);

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
      rhs.usesPseudo(i) = false;
    }

  rhs.Zeff = rhs.Z;
  return strm;
}

void QMCMolecule::readGeometry(string runfile)
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

int QMCMolecule::readPseudoPotential(string runfile)
{
  if(Natoms == 0) return 0;

  ifstream input_file(runfile.c_str());

  if(!input_file)
    {
      cerr << "ERROR: Could not open " << runfile << "!" << endl;
      exit(1);
    }

  string temp_string;
  input_file >> temp_string;
  while((temp_string != "&pseudopotential") && (input_file.eof() != 1))
    {
      input_file >> temp_string;
    }

  //No pseudopotential to read?
  if(input_file.eof() == 1) return 0;

  int gridSizes[32] = {6,14,26,38,50,74,86,110,
		       146,170,194,230,266,302,350,434,
		       590,770,974,1202,1454,1730,2030,2354,
		       2702,3074,3470,3890,4334,4802,5294,5810};
  LDGRID functions[32] = {&LD0006,&LD0014,&LD0026,&LD0038,&LD0050,&LD0074,&LD0086,&LD0110,
			  &LD0146,&LD0170,&LD0194,&LD0230,&LD0266,&LD0302,&LD0350,&LD0434,
			  &LD0590,&LD0770,&LD0974,&LD1202,&LD1454,&LD1730,&LD2030,&LD2354,
			  &LD2702,&LD3074,&LD3470,&LD3890,&LD4334,&LD4802,&LD5294,&LD5810};

  grid.allocate(Natoms);
  gridWeights.allocate(Natoms);
  gridLegendre.allocate(Natoms);  
  pseudoTitle.allocate(Natoms,2);
  Vlocal.allocate(Natoms);
  Vnonlocal.allocate(Natoms);

  /*
    I'll follow the input pattern for GAMESS. You must have one line per atom,
    which must be in the same order as the geometry.
  */
  string item;
  stringstream itemSS;
  Array1D<string> ppTypes(Natoms);
  for(int i=0; i<Natoms; i++) 
    {
      input_file >> ws;
      getline(input_file,item);
      stringstream lineSS(item);
      lineSS >> item;
      ppTypes[i] = item;
      lineSS >> ws;

      pseudoTitle(i,0) = item;
      pseudoTitle(i,1) = "";

      if(lineSS.eof()){
	// No type (e.g. NONE or GEN) was listed, this means it's all exactly
	// the same as a previous entry.
	int matched = -1;
	for(int j=0; j<i; j++){
	  if(ppTypes(j) == item){
	    matched = j;
	    usesPseudo(i) = usesPseudo(j);
	    Zeff(i) = Zeff(j);
	    Vlocal(i) = Vlocal(j);
	    Vnonlocal(i) = Vnonlocal(j);
	    grid(i) = grid(j);
	    gridWeights(i) = gridWeights(j);
	    gridLegendre(i) = gridLegendre(j);
	    break;
	  }
	}
	if(matched >= 0){
	  continue;
	} else {
	  clog << "Error: Pseudopotential " << i << " didn't list a type, and didn't match a type.\n";
	  exit(1);
	}
      }

      lineSS >> item;
      pseudoTitle(i,1) = item;
      StringManipulation::toAllUpper(item);
      if(item != "GEN" && item != "NONE"){
	//In GAMESS, the other available options correspond to things like
	//SBKJC or HW
	clog << "Error: pseudopotential for " << pseudoTitle(i,0)
	     << " has unknown type " << item << ".\n";
	exit(1);
      }

      int elec_removed = 0;
      int numl         = 0;
      int vloc_size    = 0;
      if(item == "NONE")
	{
	  Vlocal(i).allocate(0,0);
	  Vnonlocal(i).allocate(0);
	  continue;
	} 

      usesPseudo(i) = true;

      lineSS >> item;
      elec_removed = atoi(item.c_str());
      Zeff(i) = Z(i) - elec_removed;

      lineSS >> item;
      numl = atoi(item.c_str());

      input_file >> item;
      vloc_size = atoi(item.c_str());      
      //There are vloc_size gaussians used to define the local potential operator
      Vlocal(i).allocate(vloc_size,3);
      for(int vloc=0; vloc<vloc_size; vloc++)
	{
	  input_file >> item;
	  (Vlocal(i))(vloc,0) = atof(item.c_str());//coeff
	  input_file >> item;
	  (Vlocal(i))(vloc,1) = atof(item.c_str());//r^n
	  input_file >> item;
	  (Vlocal(i))(vloc,2) = atof(item.c_str());//zeta
	}

      //This pseudopotential defines numl projection operators
      Vnonlocal(i).allocate(numl);
      for(int l=0; l<numl; l++)
	{
	  input_file >> item;
	  int vnloc_size = atoi(item.c_str());      
	  //There are vnloc_size gausssians used to define the projection operator for this l
	  (Vnonlocal(i))(l).allocate(vnloc_size,3);
	  for(int vnloc=0; vnloc<vnloc_size; vnloc++)
	    {
	      input_file >> item;
	      (Vnonlocal(i))(l)(vnloc,0) = atof(item.c_str());//coeff
	      input_file >> item;
	      (Vnonlocal(i))(l)(vnloc,1) = atof(item.c_str());//r^n
	      input_file >> item;
	      (Vnonlocal(i))(l)(vnloc,2) = atof(item.c_str());//exponent
	    }
	}

      int grSize = gridSizes[gridLevel];
      
      // The Lebedev-Laikov code wants to start indexing at N=1 for some
      // unknown reason, so make the arrays 1 bigger,
      // and the first elements are meaningless. The issue might relate to
      // Fortran vs C array indexing starting at 0 or 1.
      Array1D<double> x(grSize+1);
      Array1D<double> y(grSize+1);
      Array1D<double> z(grSize+1);
      Array1D<double> w(grSize+1);
      x = 0.0;
      y = 0.0;
      w = 0.0;
      z = 0.0;
      
      grid(i).allocate(grSize,3);
      gridWeights(i).allocate(grSize);
      gridLegendre(i).allocate(numl,grSize);
      
      (functions[gridLevel])(x.array(),y.array(),z.array(),w.array());
      
      for(int gp=0; gp<grSize; gp++)
	{
	  (grid(i))(gp,0)       = x[gp+1];
	  (grid(i))(gp,1)       = y[gp+1];
	  (grid(i))(gp,2)       = z[gp+1];
	  (gridWeights(i))(gp)  = w[gp+1];
	  for(int l=0; l<numl; l++)
	    (gridLegendre(i))(l,gp) = MathFunctions::legendre(l,z[gp+1]);
	}
    }

  int usePseudo = 0;
  for(int i=0; i<Natoms; i++){
    if(usesPseudo(i)){
      usePseudo = 1;
      cout << "On atom " << i << ": replacing " << (Z(i)-Zeff(i)) << " electrons with the " << pseudoTitle(i,0) << " ecp"
	   << " using " << grid(i).dim1() << " grid points." << endl;
      /*
      cout << "Zeff(" << ppTypes(i) << ") = " << Zeff(i) << " gridSize = " << grid(i).dim1() << endl;
      cout << "Vlocal(" << ppTypes(i) << ") = \n" << Vlocal(i);
      cout << "Vnonlocal(" << ppTypes(i) << ") = " << Vnonlocal(i) << endl;
      */
    }
  }
  return usePseudo;
}

Array2D<double> QMCMolecule::getGrid(int nuc, double r, bool translate)
{
  /*
    The grid points are on a spherical shell centered at
    the nucleus. So scale the unit shell by the radius, then
    translate to the nucleus.
  */
  Array2D<double> gNuc = grid(nuc); 
  gNuc *= r;

  //No translation
  if(!translate) return gNuc;

  double x0 = Atom_Positions(nuc,0);
  double y0 = Atom_Positions(nuc,1);
  double z0 = Atom_Positions(nuc,2);
  for(int gr=0; gr<gNuc.dim1(); gr++){
    gNuc(gr,0) = gNuc(gr,0) + x0;
    gNuc(gr,1) = gNuc(gr,1) + y0;
    gNuc(gr,2) = gNuc(gr,2) + z0;
  }
  return gNuc;
}

double QMCMolecule::evaluatePotential(Array2D<double> & V, double r)
{
  double v = 0.0;
  int ngauss = V.dim1();
  for(int g=0; g<ngauss; g++)
    {
      double coeff = V(g,0);
      double n     = V(g,1);
      double zeta  = V(g,2);
      /*
	By some convention, 2 powers of r are implicit, so we
	need to subract 2 to get what we want. It might come from
	the volume integral differential element r^2 sin phi...
      */
      v += coeff * pow(r,n-2) * exp(-1.0*zeta*r*r);
    }
  return v;
}

ostream& operator <<(ostream& strm, QMCMolecule& rhs)
{
  strm << "&geometry" << endl;
  strm.unsetf(ios::fixed);
  strm.unsetf(ios::scientific);
  for(int i=0; i<rhs.Natoms; i++)
    {
      strm << setw(5) << left << rhs.Atom_Labels(i);
      strm << setw(5) << rhs.Z(i);
      for(int j=0; j<3; j++)
	{
	  if(rhs.Atom_Positions(i,j) < 0.0) strm << " -";
	  else                              strm << "  ";
	  strm << setw(20) << left << fabs(rhs.Atom_Positions(i,j));
	}
      strm << endl;
    }
  strm << "&" << endl;
  strm << right;

  bool printPseudo = false;
  for(int i=0; i<rhs.Natoms; i++)
    if(rhs.usesPseudo(i)) printPseudo = true;
  if(!printPseudo) return strm;

  strm << "&pseudopotential" << endl;
  for(int i=0; i<rhs.Natoms; i++){
    strm << rhs.pseudoTitle(i,0);
    strm << " " << rhs.pseudoTitle(i,1);

    if(!rhs.usesPseudo(i) ||
       rhs.pseudoTitle(i,1) == ""){
      strm << endl;
      continue;
    }

    //Num elec removed
    strm << " " << (rhs.Z(i) - rhs.Zeff(i));

    int numl = rhs.Vnonlocal(i).dim1();
    strm << " " << numl;
    strm << endl;

    int numG = rhs.Vlocal(i).dim1();
    strm << setw(5) << numG << endl;
    for(int ng=0; ng<numG; ng++)
      strm << setw(20) << (rhs.Vlocal(i))(ng,0)
	   << setw(5)  << (rhs.Vlocal(i))(ng,1)
	   << setw(20) << (rhs.Vlocal(i))(ng,2) << endl;
    
    for(int l=0; l<numl; l++){
      numG = (rhs.Vnonlocal(i))(l).dim1();
      strm << setw(5) << numG << endl;
      for(int ng=0; ng<numG; ng++)
	strm << setw(20) << ((rhs.Vnonlocal(i))(l))(ng,0)
	     << setw(5)  << ((rhs.Vnonlocal(i))(l))(ng,1)
	     << setw(20) << ((rhs.Vnonlocal(i))(l))(ng,2) << endl;     
    }
  }
  strm << "&" << endl;
  return strm;
}



