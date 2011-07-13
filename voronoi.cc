//
//
//

#include "voronoi.h"

#include <cmath>

using namespace lpmd;

  void SkewStart(const int n, double x, double y, double z, Vector *centers);
  void Replicate(OrthogonalCell & unitcell, ParticleSet & unitset, unsigned long nx, unsigned long ny, unsigned long nz);
  void Rotate(OrthogonalCell & unitcell, ParticleSet & unitset, Vector rotate);
  void ReplicateRotate(//
  double num, double grains, OrthogonalCell unitcell, ParticleSet unitset, Vector & cellcenter, Vector &CellColor,//
  unsigned long na, unsigned long nb, unsigned long nc, Vector rotate, BasicParticleSet & atoms, BasicCell & celda);


VoronoiGenerator::VoronoiGenerator(std::string args): Plugin("voronoi","2.0") 
{
 ParamList & params = (*this);
 //
 DefineKeyword("symbol","Ar");
 DefineKeyword("type","sc");
 DefineKeyword("a","4.08"); // Default: Gold lattice constant
 DefineKeyword("grains","2");
 // hasta aqui los valores por omision
 ProcessArguments(args); 
 spc = std::string(params["symbol"]);
 type = (*this)["type"];
 a = double(params["a"]);
 grains = int(params["grains"]);
}

VoronoiGenerator::~VoronoiGenerator() { }

void VoronoiGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = voronoi                                                  \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin is used to build a no-structured crytal using the Voronoi     \n";
 std::cout << "      method.                                                                  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol  : Set the atomic species to generate in the cell, like as Cu,    \n";
 std::cout << "                Ar, Fe, etc.                                                   \n";
 std::cout << "      type    : Set the kind of crystaline graine (fcc, bcc, hcp and sc)       \n";
 std::cout << "      a       : Crystal constant of the cell in [A].                           \n";
 std::cout << "      grains  : Number of grains to put in the final cell.                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Using the plugin to generate a crystal                                       \n";
 std::cout << " input module voronoi symbol=Ar type=fcc a=3.61 grains=10                      \n";
 std::cout << " Explanation :                                                                 \n";
 std::cout << "      With this we generate a simulation cell with 10 crystal grains of Argon  \n";
 std::cout << "      atoms unifromly distributed (using the skewstart method), each one of    \n";
 std::cout << "      this grains is a perfect FCC crystal. The final number of atoms of the   \n";
 std::cout << "      cell, depends of the number of choiced grains and the size of the simula-\n";
 std::cout << "      tion cell (more grains, more small will be, and less atomes will have).  \n"; 
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void VoronoiGenerator::Generate(lpmd::Configuration & conf) const
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & celda = conf.Cell();
 OrthogonalCell basecell(a,a,a);
 ParticleSet ps;
 double rmin=0;         // Minimum separation between atoms

 bool create_atoms = (atoms.Size() == 0);
 // NOW WE FIX THE BASIC CELL
 //-- Simple cubic lattices --//
 if (type=="sc")
 {
  rmin=0.9*a;
  const lpmd::Vector t=0*e1;
  if(create_atoms) ps.Append(Atom(spc,t));
 }
 //-- Face-centered cubic lattices --//
 else if (type=="fcc")
 {
  rmin=0.9*a/sqrt(2);
  lpmd::Vector t[4];
  t[0]=0*e1; t[1]=0.5*a*(e2+e3); t[2]=0.5*a*(e1+e3); t[3]=0.5*a*(e1+e2);
 
  if(create_atoms)
  {
   for (int i=0; i<4; ++i) ps.Append(Atom(spc,t[i]));
  }
  
 }
 //-- Body-centered cubic lattices --//
 else if (type=="bcc")
 {
  rmin=0.9*sqrt(3)*a/2.0;
  lpmd::Vector t[2];
  t[0]=0*e1; t[1]=0.5*a*(e1+e2+e3);
 
  if(create_atoms)
  {
   for (int i=0; i<2; ++i) ps.Append(Atom(spc,t[i]));
  }
 }

 unsigned long nx,ny,nz;
 double V=celda.Volume();
 double x=celda[0].Module();
 double y=celda[1].Module();
 double z=celda[2].Module();
 nx=int(2.6*pow(V/(double)grains,1.0/3.0)/a); // V/N = (nx*a^3)
 ny=nx; nz=nx;
 lpmd::Vector *centers=new Vector [grains];
 lpmd::Vector *CellColor=new Vector [grains];

 std::cout<<"\nRUNNING VORONOI PLUGIN:\n"<<std::endl;
 std::cout << " -> Cell volume: "<<V<<" cubic angstroms."<<std::endl;
 std::cout << " -> Creating "<<grains<<" grains. Estimated average grain diameter: "<<pow(V/grains,1.0/3.0)<<" angstroms."<<std::endl;
 std::cout << " -> Unitary "<< type <<" cell replicated "<<nx<<" times in each axis..."<<std::endl;

 // CHOOSE CENTERS AND REPLICATE UNITARY CELLS
 SkewStart(grains, x, y, z, centers);
 for (int i=0; i<grains; ++i)
 {
  Vector rotate=2*M_PI*drand48()*e1+2*M_PI*drand48()*e2+2*M_PI*drand48()*e3;
  ReplicateRotate(i, grains, basecell, ps, centers[i], CellColor[i], nx, ny, nz, rotate, atoms, celda);
 } 
 
 std::cout << " -> Replicated unitary cells were rotated and moved into the cell with Skew-Start method..."<<std::endl;

 //-------------------- 1ST ELIMINATION -------------------//
 // OUTSIDE ELIMINATION: Eliminate the atoms out of the cell
 std::cout << " -> Eliminating atoms out of the cell..."<<std::endl;
 for (long i=0;i<atoms.Size();++i)
 {
  bool kill=false;
  Vector pos = atoms[i].Position();
  if (pos[0]<0 || pos[0]>x) {atoms.Delete(i); kill=true;}
  else if (pos[1]<0 || pos[1]>y) {atoms.Delete(i); kill=true;}
  else if (pos[2]<0 || pos[2]>z) {atoms.Delete(i); kill=true;}
  
  if (kill) i--;
 }

 //-------------------- 2ND ELIMINATION -------------------//
 // INTERSECTIONS ELIMINATION
 std::cout<< " -> Separating grains..."<<std::endl;
 for (long int i=0; i<atoms.Size(); ++i)
 {
  bool eliminated=false;
  for(int n=0; n<grains; ++n)
  {
   for(int m=0; m<grains; ++m)
   {
    if(n!=m)
    {
     Vector sep=centers[m]-centers[n];
     Vector pos=atoms[i].Position()-centers[n];
     Vector atmclr=ColorHandler::ColorOfAtom(atoms[i]);
     if ( Dot(pos,sep/sep.Module())>0.5*sep.Module() && atmclr==CellColor[n])
     {
      atoms.Delete(i);
      eliminated=true; i--;
     }
    }
    if (eliminated) break;
   }
   if (eliminated) break;
  }
 }
 //-------------------- 3D ELIMINATION -------------------//
 // PAIRS ELIMINATION: NOW THAT THE CELL IS FULL OF CELLS FILLED WITH ATOMS, WE ELIMINATE THE CLOSEST ATOMS
 std::cout << " -> Eliminating closest atoms..."<<std::endl;
 for (long i=0; i<atoms.Size(); ++i)
 {
  for (long j=i+1; j<atoms.Size(); ++j)
  {
   double dis=celda.Displacement(atoms[i].Position(), atoms[j].Position()).Module();
   if(dis<rmin)
   {
    atoms.Delete(i);
    i--; break;
   }
  }
 }

 // Updating positions (impose periodic boudary conditions)
// for (unsigned long i=0; i<atoms.Size(); ++i) atoms[i].Position()=atoms[i].Position()+1.5*(x*e1+y*e2+z*e3);

 delete [] CellColor;
 delete [] centers;
 std::cout << " -> A cell of "<<atoms.Size()<<" atoms was created."<<std::endl;
 std::cout<<"\nREADY.\n"<<std::endl;
 
}

//------------------------------------------------------------------------------------------------------------//
	
void SkewStart(int n, double x, double y, double z, Vector *centers)
{
 OrthogonalCell celda(x, y, z);
 ParticleSet atomos(n);
 int h, k, l;
 double dx, dy, dz;
 h = int(pow(double(n), 2.0/3.0));
 k = int(pow(double(n), 1.0/3.0));
 l = 1;
 dx = h / double(n);
 dy = k / double(n);
 dz = l / double(n);
 for (long i=0;i<n;++i)
 {
  atomos[i].Position()=celda.FittedInside(celda.Cartesian(Vector(dx*double(i)+0.5, dy*double(i)+0.5, dz*double(i)+0.5)));
  centers[i]=atomos[i].Position();
 }
}


void Replicate(OrthogonalCell & unitcell,ParticleSet & unitset,unsigned long nx,unsigned long ny,unsigned long nz)
{
 int N=unitset.Size();
 for (unsigned long i=1; i<=nx; ++i)
  for (int j=0; j<N; j++)
   unitset.Append(Atom(unitset[j].Z(),unitset[j].Position()+unitcell[0]*i));

 N=unitset.Size();
 for (unsigned long i=1; i<=ny; ++i)
  for (int j=0; j<N; j++)
   unitset.Append(Atom(unitset[j].Z(),unitset[j].Position()+unitcell[1]*i));
 
 N=unitset.Size();
 for (unsigned long i=1; i<=nz; ++i)
  for (int j=0; j<N; j++)
   unitset.Append(Atom(unitset[j].Z(),unitset[j].Position()+unitcell[2]*i));

 unitcell[0]=unitcell[0]*nx;
 unitcell[1]=unitcell[1]*ny;
 unitcell[2]=unitcell[2]*nz;

}

void Rotate(OrthogonalCell & unitcell, ParticleSet & unitset, Vector rotate)
{
 // Euler Rotation Matrix
 double rotmat[3][3];
 // Eulerian Angles
 double phi=rotate[0], psi=rotate[1], theta=rotate[2];
// double phi=M_PI/4.0, psi=0, theta=M_PI/2;
 // Eulerian Matrix
 rotmat[0][0] = (cos(phi)*cos(psi)-cos(theta)*sin(phi)*sin(psi));
 rotmat[0][1] = (cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi));
 rotmat[0][2] = (sin(theta)*sin(psi));
 rotmat[1][0] = (-cos(theta)*cos(psi)*sin(phi)-cos(phi)*sin(psi));
 rotmat[1][1] = ( cos(theta)*cos(phi)*cos(psi)-sin(phi)*sin(psi));
 rotmat[1][2] = (cos(psi)*sin(theta));
 rotmat[2][0] = (sin(theta)*sin(phi));
 rotmat[2][1] = (-cos(phi)*sin(theta));
 rotmat[2][2] = ( cos(theta));

 for (long int i=0; i<3; ++i)
 {
  double v0=unitcell[i][0], v1=unitcell[i][1], v2=unitcell[i][2];
  Vector newvec;
  for (int j=0; j<3; j++) newvec[j] = v0*rotmat[j][0]+v1*rotmat[j][1]+v2*rotmat[j][2];
  unitcell[i] = newvec;
 }

 long int N=unitset.Size();
 for (long int i=0; i<N; ++i)
 {
  Vector newvec;
  double v0=unitset[i].Position()[0], v1=unitset[i].Position()[1], v2=unitset[i].Position()[2];
  for (int j=0; j<3; j++) newvec[j] = v0*rotmat[j][0]+v1*rotmat[j][1]+v2*rotmat[j][2];
  unitset[i].Position() = newvec;
 }

}




void ReplicateRotate(//
  double num, double grains, OrthogonalCell unitcell, ParticleSet unitset, Vector & cellcenter, Vector &CellColor,//
  unsigned long na, unsigned long nb, unsigned long nc,//
  Vector rotate, BasicParticleSet & atoms, BasicCell & celda)
{
 ParticleSet us=unitset;
 OrthogonalCell uc=unitcell;
 Replicate(uc, us, na,nb,nc);
 Rotate(uc, us, rotate);

 Vector centro=0.5*(uc[0]+uc[1]+uc[2]);
 unsigned long N = us.Size();
 CellColor = ColorFromScalar(num/grains);
 for(unsigned long i=0;i<N;++i)
 {
  Vector vct=us[i].Position()+(cellcenter-centro);
  Atom at(us[i].Z(),vct);
  atoms.Append(Atom(us[i].Z(),vct));
  ColorHandler::ColorOfAtom(atoms[atoms.Size()-1]) = CellColor;
 }



 //-------------------- 1ST ELIMINATION -------------------//
 // OUTSIDE ELIMINATION: Eliminate the atoms out of the cell
 std::cout << " -> Eliminating atoms out of the cell of the grain number "<< num+1<<"..."<<std::endl;
 double x=celda[0].Module();
 double y=celda[1].Module();
 double z=celda[2].Module();
 for (long i=0;i<atoms.Size();++i)
 {
  bool kill=false;
  Vector pos = atoms[i].Position();
  if (pos[0]<0 || pos[0]>x) {atoms.Delete(i); kill=true;}
  else if (pos[1]<0 || pos[1]>y) {atoms.Delete(i); kill=true;}
  else if (pos[2]<0 || pos[2]>z) {atoms.Delete(i); kill=true;}
  
  if (kill) i--;
 }


}



// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new VoronoiGenerator(args); }
void destroy(Plugin * m) { delete m; }

