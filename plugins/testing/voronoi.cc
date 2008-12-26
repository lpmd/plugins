//
//
//

#include "voronoi.h"
#include "plugincommon.h"

#include <lpmd/atom.h>
#include <lpmd/simulationcell.h>

#include <cmath>

using namespace lpmd;
double rmin;			// minimum separation between atoms

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
 std::cout << " Module Name        = voronoi                                                \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para generar una celda compuesta por celdas       \n";
 std::cout << " mas pequenas (nanoestructuras).                                               \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol  : Especifica el simbolo atomico de la especie a generar.         \n";
 std::cout << "      type    : Especifica el tipo de celda base (sc, bcc, fcc)                \n";
 std::cout << "      a       : Especifica el tamano de la celda type                          \n";
 std::cout << "      cells   : Especifica la cantidad de celdas que se replicaran             \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=voronoi symbol=Ar a=2 cells=10 type=fcc                          \n";
 std::cout << "      De esta forma podemos generar una celda inicial, que contiene 10 celdas  \n";
 std::cout << " pequenas de argon, hechas a partir de la replicacion de una celda fcc,        \n";
 std::cout << " eliminando los atomos que se encuentren a una distancia menor que cierta      \n";
 std::cout << " distancia (que depende de type).                                              \n";
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
 for (int i=0; i<3; i++) basecell.SetVector(i, base[i]);     	
 // NOW WE FIX THE BASIC CELL
 //-- Simple cubic lattices --//
 if (type=="sc")
 {
  rmin=a;
  Vector vec(0,0,0);
  Atom atm(18,vec);
  basecell.AppendAtom(atm);
 }
 //-- Face-centered cubic lattices --//
 else if (type=="fcc")
 {
  rmin=a/sqrt(2);
  Atom atm[4];
  Vector t[4];
  t[0]=0*e1; t[1]=0.5*a*(e2+e3); t[2]=0.5*a*(e1+e3); t[3]=0.5*a*(e1+e2);
  for (int i=0; i<4; i++){ atm[i].SetSpc(18); atm[i].SetPos(t[i]); basecell.AppendAtom(atm[i]);}
 }
 //-- Body-centered cubic lattices --//
 else if (type=="bcc")
 {
  rmin=0.1;//0.9*sqrt(3)*a/2.0;
  Atom atm[4];
  Vector t[4];
  t[0]=0*e1; t[1]=0.5*a*(-e1+e2+e3); t[2]=0.5*a*(e1-e2+e3); t[3]=0.5*a*(e1+e2-e3);
  for (int i=0; i<4; i++){ atm[i].SetSpc(18); atm[i].SetPos(t[i]); basecell.AppendAtom(atm[i]);}
//  for (int i=0; i<3; i++) std::cout << "el vect "<<i<<" de la basecell es "<<basecell.GetVector(i)<<"\n";
//  for (int i=0; i<4; i++){ std::cout<<"atomo "<<i<<" de la basecell esta en: "<<basecell.GetAtom(i).Position()<<"\n";}
 }

 unsigned long nx,ny,nz;
 Vector cellcenter;
 Vector *centers=new Vector [Ncell];

 randomize();
 for (int i=0; i<Ncell; i++)
 {
  // Random "cellcenter" and amount of replications "nx", "ny", "nz".
//  Vector cellcenter(50-4.3*i,50,50+4.3*i);
  cellcenter=(dazar(0,1)*sc.GetVector(0)+dazar(0,1)*sc.GetVector(1)+dazar(0,1)*sc.GetVector(2));
//  cellcenter=0.4*(sc.GetVector(0)+sc.GetVector(1)+sc.GetVector(2));
//  std::cout<<"centro de celda: "<<cellcenter<<"\n";
  centers[i]=cellcenter;
  nx=int(sc.GetVector(0).Mod()/a); ny=nx; nz=nx;
//  std::cout << "voy a replicar "<< nx << " por cada eje"<< std::endl;
  // We put a new replicated cell in "cellcenter"
//  std::cout << "antes, la celda grande tenia un total de "<<sc.Size() << "atomos"<<std::endl;
  ReplicateRotate(basecell, cellcenter, nx, ny, nz, sc);
//  std::cout << "ahora, la celda grande tiene un total de "<<sc.Size() << "atomos"<<std::endl;
//  for (int i=0; i<sc.Size(); i++) std::cout << "posicion atomo "<<i<<"="<<sc.GetAtom(i).Position()<<"\n";
 }
 // FIRST ELIMINATION: NOW THAT THE CELL IS FULL OF CELLS FILLED WITH ATOMS, WE ELIMINATE THE CLOSEST ATOMS
 for (int i=0; i<sc.Size(); i++)
 {
//  std::cout<<"F. Elim: Particula "<< i <<".\n";
  Vector posi=sc.GetAtom(i).Position();
  for (int j=0; j<sc.Size(); j++)
  {
   Vector posj=sc.GetAtom(j).Position();
   if (i!=j && (posi-posj).Mod()<rmin){ /*std::cout<<"eliminando part "<<j<<", ya que rmin="<<rmin<<", y posi-posj="<<fabs((posi-posj).Mod())<<", posi="<<posi<<", posj="<<posj<<", i="<<i<<",j="<<j<<"\n";*/ sc.DeleteAtom(j); j--;}
  }
 }
 // SECOND ELIMINATION: NOW THAT THE CLOSEST ATOMS WERE DELETED, WE CHOOSE THE SUBCELL WHERE THEY BELONG
 for (int j=0; j<sc.Size(); j++)
 {
//  std::cout<<"S. Elim: Particula "<< j <<".\n";
  bool eliminated=false;
  Vector posj=sc.GetAtom(j).Position();
  for(int n=0; n<Ncell; n++)
  {
   for(int m=0; m<Ncell; m++)
   {
    if(n!=m)
    {
     double jn=(posj-centers[n]).Mod();
     double jm=(posj-centers[m]).Mod();
     if (fabs(jn-jm)<1.0e-1){
//      std::cout<<"el atomo "<<j<<" esta a "<<jn<<" de la celda "<<n <<" y a "<<jm<<" de la "<<m<<". Eliminando.\n";
      sc.DeleteAtom(j); eliminated=true;
//      std::cout<<"atomo "<<j<<" eliminado.\n";
      j--;
     }
    }
    if (eliminated) break;
   }
   if (eliminated) break;
  }
 }
// std::cout << "Ahora la celda quedo con "<<sc.Size()<< " atomos.\n";
 delete [] centers;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VoronoiGenerator(args); }
void destroy(Module * m) { delete m; }
