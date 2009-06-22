//
//
//
#include "voronoi.h"

#include <lpmd/particleset.h>
#include <lpmd/configuration.h>
#include <lpmd/orthogonalcell.h>
//#include "plugincommon.h"

//#include <algorithm>
//#include <functional>

#include <cmath>

using namespace lpmd;
double rmin;					// Minimum separation between atoms

void SkewStart(int n, double x, double y, double z, Vector *centers)
{
 // FIXME: int() corrige los warning, pero hay que chequear si es lo correcto o no
 OrthogonalCell celda(x, y, z);
 ParticleSet atomos;
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
  atomos[i].Position() = celda.ScaleByCell(Vector(dx*double(i), dy*double(i), dz*double(i)));
 }
 for (long i=0;i<n;++i)
 {
//  centers[i]=x*cell[i].Position()[0]*e1+y*cell[i].Position()[1]*e2+z*cell[i].Position()[2]*e3;
  centers[i]=atomos[i].Position();
 }
}

VoronoiGenerator::VoronoiGenerator(std::string args): Plugin("voronoi","2.0") 
{
 ParamList & params = (*this);
 //
 DefineKeyword("symbol","Ar");
 DefineKeyword("type","sc");
 DefineKeyword("a","3.61");
 DefineKeyword("grains","2");
 // hasta aqui los valores por omision
 ProcessArguments(args); 
// spc = ElemNum((*this)["symbol"]);
 spc = std::string(params["symbol"]);
 type = (*this)["type"];
 a = double(params["a"]);
 grains = int(params["grains"]);
}

VoronoiGenerator::~VoronoiGenerator() { }

void VoronoiGenerator::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This module is used to generate a nanostructured crystal.                \n";
 std::cout << " General Options   >>                                                         	\n";
 std::cout << "      symbol  : Specifies the atomic symbol of the species to generate like    \n";
 std::cout << "                Cu, Ar, Fe, etc.                                               \n";
 std::cout << "      type    : Specifies the cell type (sc, bcc, fcc).                        \n";
 std::cout << "      a       : Lattice constant.                                              \n";
 std::cout << "      grains  : Number of grains to put in the simulation cell                 \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Module calling :                                                              \n";
 std::cout << " input module=voronoi symbol=Ar type=fcc a=3.61 grains=10                      \n";
 std::cout << " Explanation :                                                                 \n";
 std::cout << "      This generates a simulation cell that contains 10 argon grains uniformly \n";
 std::cout << "      distributed (using skewstart method), each one of them being a fcc       \n";
 std::cout << "      cutted crystal. The number of resulting atoms depends on the number of		\n";
 std::cout << "      grains you choose and the size of the simulation cell (the more grains   \n";
 std::cout << "      you put, the smaller they become, and less atoms you have).               \n";
}

void VoronoiGenerator::Generate(lpmd::Configuration & conf) const
{
 OrthogonalCell basecell(a,a,a);  // cell es la celda de simulacion
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & celda = conf.Cell();
 bool create_atoms = (atoms.Size() == 0);
 // NOW WE FIX THE BASIC CELL
 //-- Simple cubic lattices --//
 if (type=="sc")
 {
  rmin=0.9*a;
  const lpmd::Vector t=0*e1;
  if(create_atoms) atoms.Append(Atom(spc,t));
 }
 //-- Face-centered cubic lattices --//
 else if (type=="fcc")
 {
  rmin=0.9*a/sqrt(2);
  lpmd::Vector t[4];
  t[0]=0*e1; t[1]=0.5*a*(e2+e3); t[2]=0.5*a*(e1+e3); t[3]=0.5*a*(e1+e2);
 
  if(create_atoms)
  {
   for (int i=0; i<4; i++) atoms.Append(Atom(spc,t[i]));
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
   for (int i=0; i<2; i++) atoms.Append(Atom(spc,t[i]));
  }
 }

 unsigned long nx,ny,nz;
 double V=celda.Volume();
// double x=celda[0].Module();
// double y=celda[1].Module();
//| double z=celda[2].Module();
 nx=int(2.2*pow(V/(double)grains,1.0/3.0)/a); ny=nx; nz=nx;
// lpmd::Vector *centers=new Vector [grains];
// lpmd::Vector *CellColor=new Vector [grains];

 std::cout<<"\nRUNNING VORONOI PLUGIN:"<<std::endl;

 // CHOOSE CENTERS AND REPLICATE CELLS
// SkewStart(grains, x, y, z, centers);
 for (int i=0; i<grains; i++)
 {
  Vector rotate=2*M_PI*drand48()*e1+2*M_PI*drand48()*e2+2*M_PI*drand48()*e3;
  // We put a new replicated cell in "cellcenter"
//  ReplicateRotate(basecell, centers[i], CellColor[i], nx, ny, nz, rotate, sc);
 }
/*
 // OUTSIDE ELIMINATION: ELIMINATE OUTSIDE ATOMS
 std::cout << "Elimination of atoms out of the cell..."<<std::endl;
 for (long i=0;i<atoms.Size();i++)
 {
  bool kill=false;
  Vector pos = atoms[i].Position();
  if (pos[0]<0 || pos[0]>x) {atoms.Delete(i); kill=true;}
  else if (pos[1]<0 || pos[1]>y) {atoms.Delete(i); kill=true;}
  else if (pos[2]<0 || pos[2]>z) {atoms.Delete(i); kill=true;}
  
  if (kill) i--;
 }


 // ELIMINATION BY CUTTING PLANE
 std::cout<< "Separating grains..."<<std::endl;
 for (unsigned long int i=0; i<atoms.Size(); i++)
 {
  bool eliminated=false;
  for(int n=0; n<Ncell; n++)
  {
   for(int m=0; m<Ncell; m++)
   {
    if(n!=m){
     Vector sep=centers[m]-centers[n];
     Vector pos=sc[i].Position()-centers[n];
     Vector atmclr=sc[i].Color();
     //FIXME : Sobrecarga de operador == en vector, ¿necesario?
//     if ( Dot(pos,sep/sep.Module())>0.5*sep.Module() && atmclr==CellColor[n])
     {
      sc.Destroy(&sc[i]); eliminated=true; i--;
     }
    }
    if (eliminated) break;
   }
   if (eliminated) break;
  }
 }


 // PAIRS ELIMINATION: NOW THAT THE CELL IS FULL OF CELLS FILLED WITH ATOMS, WE ELIMINATE THE CLOSEST ATOMS
 std::cout << "Eliminating closest atoms..."<<std::endl;
 for (unsigned long i=0; i<atoms.Size(); i++)
 {
  for (unsigned long j=i+1; j<atoms.Size(); j++)
  {
   if(sc.Distance(i,j)<rmin)
   {
    sc.Destroy(&sc[i]);
    i=0; break;
   }
  }
 }
*/
 // Updating positions (impose periodic boudary conditions)
// for (unsigned long i=0; i<atoms.Size(); i++) atoms[i].Position()=atoms[i].Position()+1.5*(x*e1+y*e2+z*e3);

// delete [] CellColor;
// delete [] centers;
 std::cout<<"Ready."<<std::endl;
 
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new VoronoiGenerator(args); }
void destroy(Plugin * m) { delete m; }

