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
  rmin=0.9*a;
  Atom atm;
  atm.SetSpc(spc); atm.SetPos(0*e1); basecell.AppendAtom(atm);
 }
 //-- Face-centered cubic lattices --//
 else if (type=="fcc")
 {
  rmin=0.9*a/sqrt(2);
  Atom atm[4];
  Vector t[4];
  t[0]=0*e1; t[1]=0.5*a*(e2+e3); t[2]=0.5*a*(e1+e3); t[3]=0.5*a*(e1+e2);
  for (int i=0; i<4; i++){ atm[i].SetSpc(spc); atm[i].SetPos(t[i]); basecell.AppendAtom(atm[i]);}
 }
 //-- Body-centered cubic lattices --//
 else if (type=="bcc")
 {
  rmin=0.9*sqrt(3)*a/2.0;
  Atom atm[2];
  Vector t[2];
  t[0]=0*e1; t[1]=0.5*a*(e1+e2+e3);
  for (int i=0; i<2; i++){ atm[i].SetSpc(spc); atm[i].SetPos(t[i]); basecell.AppendAtom(atm[i]);}
//  for (int i=0; i<3; i++) std::cerr << "el vect "<<i<<" de la basecell es "<<basecell.GetVector(i)<<"\n";
//  for (int i=0; i<4; i++){ std::cerr<<"atomo "<<i<<" de la basecell tiene color: "<<basecell.GetAtom(i).Color()<<"\n";}
//  for (int i=0; i<4; i++){ std::cerr<<"atomo "<<i<<" de la basecell esta en: "<<basecell.GetAtom(i).Position()<<"\n";}
 }

 unsigned long nx,ny,nz;
 Vector cellcenter;
 Vector *centers=new Vector [Ncell];

 randomize();
 for (int i=0; i<Ncell; i++)
 {
  // Random "cellcenter" and amount of replications "nx", "ny", "nz".
  cellcenter=(dazar(0,1)*sc.GetVector(0)+dazar(0,1)*sc.GetVector(1)+dazar(0,1)*sc.GetVector(2));
//  cellcenter=(25-4.3*i)*e1+25*e2+(25+4.3*i)*e3;
//  cellcenter=0.0*(sc.GetVector(0)+sc.GetVector(1)+sc.GetVector(2));
//  std::cerr<<"centro de celda: "<<cellcenter<<"\n";
  centers[i]=cellcenter;
  nx=int(sc.GetVector(0).Mod()/a); ny=nx; nz=nx;
//  std::cerr << "voy a replicar "<< nx << " por cada eje"<< std::endl;
  // We put a new replicated cell in "cellcenter"
//  std::cerr << "antes, la celda grande tenia un total de "<<sc.Size() << " atomos"<<std::endl;
  ReplicateRotate(basecell, cellcenter, nx, ny, nz, sc);
//  std::cerr << "ahora, la celda grande tiene un total de "<<sc.Size() << " atomos"<<std::endl;
 }
// for (int i=0; i<sc.Size(); i++) std::cerr << "posicion atomo "<<i<<"="<<sc.GetAtom(i).Position()<<"\n";

 // FIRST ELIMINATION: NOW THAT THE CELL IS FULL OF CELLS FILLED WITH ATOMS, WE ELIMINATE THE CLOSEST ATOMS
 for (int i=0; i<sc.Size(); i++)
 {
//  std::cerr<<"F. Elim: Particula "<< i <<".\n";
  bool eliminated=false;
  Vector posi=sc.GetAtom(i).Position();
  for (int j=0; j<sc.Size(); j++)
  {
   Vector posj=sc.GetAtom(j).Position();
   if (i!=j && (posi-posj).Mod()<rmin)
   {
    int coin=iazar(0,1);
//    std::cerr<<"para i<j, i="<<i<<",j="<<j<<", eliminamos part ";
    if (coin==0)     { sc.DeleteAtom(j); /*std::cerr<<j*/; j--;}
    else if (coin==1){ sc.DeleteAtom(i); /*std::cerr<<i*/; i--; eliminated=true;}
//    std::cerr<<", ya que rmin="<<rmin<<" y posi-posj="<<fabs((posi-posj).Mod())<<":\n posi="<<posi<<",\n posj="<<posj<<"\n";
   }
   if (eliminated) break;
  }
 }
// std::cerr << "despues de la primera eliminacion, la sc tiene un total de "<<sc.Size() << " atomos"<<std::endl;
 // SECOND ELIMINATION: NOW THAT THE CLOSEST ATOMS WERE DELETED, WE CHOOSE THE SUBCELL WHERE THEY BELONG
 for (int i=0; i<sc.Size(); i++)
 {
//  std::cerr<<"S. Elim: Particula "<< j <<".\n";
  bool eliminated=false;
  Vector posj=sc.GetAtom(i).Position();
  for(int n=0; n<Ncell; n++)
  {
   for(int m=0; m<Ncell; m++)
   {
    if(n!=m)
    {
     double in=(posj-centers[n]).Mod();
     double im=(posj-centers[m]).Mod();
     if (fabs(in-im)<1.0e-1)
     {
//      std::cerr<<"el atomo "<<j<<" esta a "<<jn<<" de la celda "<<n <<" y a "<<jm<<" de la "<<m<<". Eliminando.\n";
      sc.DeleteAtom(i); eliminated=true; i--;
//      std::cerr<<"atomo "<<j<<" eliminado.\n";
     }
    }
    if (eliminated) break;
   }
   if (eliminated) break;
  }
 }
// std::cerr << "Ahora la celda quedo con "<<sc.Size()<< " atomos.\n";
 delete [] centers;
 // Updating positions (impose periodic boudary conditions)
 for (long i=0; i<sc.Size(); i++) sc.SetPosition(i,sc.GetAtom(i).Position());
// for (long i=0; i<sc.Size(); i++) std::cerr << "con cbp, posicion atomo "<<i<<"="<<sc.GetAtom(i).Position()<<"\n";
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VoronoiGenerator(args); }
void destroy(Module * m) { delete m; }
