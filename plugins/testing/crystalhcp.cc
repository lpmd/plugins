//
//
//

#include "crystalhcp.h"

#include <lpmd/atom.h>
#include <lpmd/simulationcell.h>

using namespace lpmd;

HCPGenerator::HCPGenerator(std::string args): Module("crystalhcp")
{
 AssignParameter("nx", "1");
 AssignParameter("ny", "1");
 AssignParameter("nz", "1");
 AssignParameter("symbol", "H");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 spc = ElemNum(GetString("symbol"));
 nx = GetInteger("nx");
 ny = GetInteger("ny");
 nz = GetInteger("nz");
}

HCPGenerator::~HCPGenerator() { }

void HCPGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = crystalhcp                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para crear celdas del tipo hexagonal close-packed \n";
 std::cout << "      (HCP).                                                                   \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << "      nz            : Repeticiones en la direccion Z.                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystalhcp symbol=Fe nx=3 ny=3 nz=3                                     \n";
 std::cout << " input crystalhcp Fe 3 3 3                                                     \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo HCP en la simulacion. \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string HCPGenerator::Keywords() const
{
 return "symbol nx ny nz";
}

void HCPGenerator::Generate(SimulationCell & sc) const
{
 Vector p;
 long int cc = 0;
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 double az = 1.0/double(nz);
 std::cerr << "-> Generating HCP cell, c/a ratio is " << (nx*sc.GetVector(2).Mod())/(nz*sc.GetVector(0).Mod()) << '\n';
 for (long k=0;k<nz;++k)
 {
  for (long j=0;j<ny;++j)
  {
   for (long i=0;i<nx;++i)
   {
    p = Vector((double(i)+(2.0/3.0))*ax, (double(j)+(1.0/3.0))*ay, (double(k)+(3.0/4.0))*az);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+(1.0/3.0))*ax, (double(j)+(2.0/3.0))*ay, (double(k)+(1.0/4.0))*az);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
   }
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new HCPGenerator(args); }
void destroy(Module * m) { delete m; }

