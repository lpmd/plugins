//
//
//

#include "voronoi.h"
#include "plugincommon.h"
#include <algorithm>
#include <functional>

#include <lpmd/atom.h>
#include <lpmd/simulationcell.h>

#include <cmath>

using namespace lpmd;
double rmin;					// Minimum separation between atoms

void SkewStart(int n, double x, double y, double z, Vector *centers)
{
 // FIXME: int() corrige los warning, pero hay que chequear si es lo correcto o no
 SimulationCell cell(int(x), int(y), int(z), true, true, true);
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
  cell.Create(new Atom(18, Vector()));
  cell.SetFracPosition(i, Vector(dx*double(i)+0.25, dy*double(i)+0.25, dz*double(i)+0.25));
 }
 for (long i=0;i<n;++i)
 {
  centers[i]=x*cell[i].Position()[0]*e1+y*cell[i].Position()[1]*e2+z*cell[i].Position()[2]*e3;
 }
}

VoronoiGenerator::VoronoiGenerator(std::string args): Module("voronoi") 
{
 AssignParameter("type","sc");
 ProcessArguments(args); 
 spc = ElemNum(GetString("symbol"));
 a = GetDouble("a");
 Ncell = GetInteger("cells");
 type = GetString("type");
}

VoronoiGenerator::~VoronoiGenerator() { }

void VoronoiGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = voronoi                                                	\n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This module is used to generate a nanostructured crystal.                \n";
 std::cout << " General Options   >>                                                         	\n";
 std::cout << "      symbol  : Specifies the atomic symbol of the species to generate like    \n";
 std::cout << "                Cu, Ar, Fe, etc.                                               \n";
 std::cout << "      type    : Specifies the cell type (sc, bcc, fcc).                        \n";
 std::cout << "      a       : Lattice constant.                                              \n";
 std::cout << "      grains  : Number of grains to put in the simulation cell                 \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Module calling :                                                              \n";
 std::cout << " input module=voronoi symbol=Ar type=fcc a=3.61 grains=10                      \n";
 std::cout << " Explanation :                                                                 \n";
 std::cout << "      This generates a simulation cell that contains 10 argon grains uniformly \n";
 std::cout << "      distributed (using skewstart method), each one of them being a fcc       \n";
 std::cout << "      cutted crystal. The number of resulting atoms depends on the number of		\n";
 std::cout << "      grains you choose and the size of the simulation cell (the more grains   \n";
 std::cout << "      you put, the smaller they become, and less atoms you have).               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string VoronoiGenerator::Keywords() const
{
 return "symbol a cells type";
}

void VoronoiGenerator::Generate(SimulationCell & sc) const
{
 SimulationCell basecell(1, 1, 1, false, false, false);  // cell es la celda de simulacion
 Vector base[3];
 base[0]=a*e1; 	base[1]=a*e2; 	base[2]=a*e3;
 // WE SET THE VECTORS OF THE CELL
 for (int i=0; i<3; i++) basecell.GetCell()[i] = base[i];
 // NOW WE FIX THE BASIC CELL
 //-- Simple cubic lattices --//
 if (type=="sc")
 {
  rmin=0.9*a;
  Atom atm;
  atm.SetSpc(spc); atm.SetPos(0*e1); basecell.Create(new Atom(atm));
 }
 //-- Face-centered cubic lattices --//
 else if (type=="fcc")
 {
  rmin=0.9*a/sqrt(2);
  Atom atm[4];
  Vector t[4];
  t[0]=0*e1; t[1]=0.5*a*(e2+e3); t[2]=0.5*a*(e1+e3); t[3]=0.5*a*(e1+e2);
  for (int i=0; i<4; i++){ atm[i].SetSpc(spc); atm[i].SetPos(t[i]); basecell.Create(new Atom(atm[i]));}
 }
 //-- Body-centered cubic lattices --//
 else if (type=="bcc")
 {
  rmin=0.9*sqrt(3)*a/2.0;
  Atom atm[2];
  Vector t[2];
  t[0]=0*e1; t[1]=0.5*a*(e1+e2+e3);
  for (int i=0; i<2; i++){ atm[i].SetSpc(spc); atm[i].SetPos(t[i]); basecell.Create(new Atom(atm[i]));}
 }

 unsigned long nx,ny,nz;
 const Cell & scell = sc.GetCell();
 double V=scell.Volume();
 double x=scell[0].Module();
 double y=scell[1].Module();
 double z=scell[2].Module();
 nx=int(2.2*pow(V/(double)Ncell,1.0/3.0)/a); ny=nx; nz=nx;
 Vector *centers=new Vector [Ncell];
 Vector *CellColor=new Vector [Ncell];

 randomize();
 // CHOOSE CENTERS AND REPLICATE CELLS
 SkewStart(Ncell, x, y, z, centers);
 for (int i=0; i<Ncell; i++)
 {
  Vector rotate=dazar(0,2*M_PI)*e1+dazar(0,2*M_PI)*e2+dazar(0,2*M_PI)*e3;
  // We put a new replicated cell in "cellcenter"
  ReplicateRotate(basecell, centers[i], CellColor[i], nx, ny, nz, rotate, sc);
 }

 // OUTSIDE ELIMINATION: ELIMINATE OUTSIDE ATOMS
 for (unsigned long i=0;i<sc.size();i++)
 {
  bool kill=false;
  Vector pos = sc[i].Position();
  if (pos[0]<0 || pos[0]>sc.GetCell()[0].Module()) {sc.Destroy(&sc[i]); kill=true;}
  else if (pos[1]<0 || pos[1]>sc.GetCell()[1].Module()) {sc.Destroy(&sc[i]); kill=true;}
  else if (pos[2]<0 || pos[2]>sc.GetCell()[2].Module()) {sc.Destroy(&sc[i]); kill=true;}
  
  if (kill) i--;
 }


 // ELIMINATION BY CUTTING PLANE
 for (unsigned long int i=0; i<sc.size(); i++)
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
 for (unsigned long i=0; i<sc.size(); i++)
 {
  for (unsigned long j=i+1; j<sc.size(); j++)
  {
   if(sc.Distance(i,j)<rmin)
   {
    sc.Destroy(&sc[i]);
    i=0; break;
   }
  }
 }

 // Updating positions (impose periodic boudary conditions)
 for (unsigned long i=0; i<sc.size(); i++) sc.SetPosition(i,sc[i].Position()+1.5*(x*e1+y*e2+z*e3));

 delete [] CellColor;
 delete [] centers;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VoronoiGenerator(args); }
void destroy(Module * m) { delete m; }
