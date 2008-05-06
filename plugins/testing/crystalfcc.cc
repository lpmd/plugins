//
//
//

#include "crystalfcc.h"

#include <lpmd/atom.h>
#include <lpmd/simulationcell.h>

using namespace lpmd;

FCCGenerator::FCCGenerator(std::string args): Module("crystalfcc")
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

FCCGenerator::~FCCGenerator() { }

void FCCGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = crystalfcc                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para crear celdas del tipo cubica centrada en las \n";
 std::cout << "      caras (FCC).                                                             \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << "      nz            : Repeticiones en la direccion Z.                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystalfcc symbol=Ar nx=3 ny=3 nz=3                                     \n";
 std::cout << " input crystalfcc Ar 3 3 3                                                   \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo FCC en la simulacion. \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string FCCGenerator::Keywords() const
{
 return "symbol nx ny nz";
}

void FCCGenerator::Generate(SimulationCell & sc) const
{
 Vector p;
 long int cc = 0;
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 double az = 1.0/double(nz);
 for (long k=0;k<nz;++k)
 {
  for (long j=0;j<ny;++j)
  {
   for (long i=0;i<nx;++i)
   {
    p = Vector((double(i)+0.5)*ax, (double(j)+0.5)*ay, (double(k)+0.5)*az);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+0.5)*ax, double(j)*ay, (double(k)+1.0)*az);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+1.0)*ax, double(j)*ay, (double(k)+0.5)*az);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+1.0)*ax, (double(j)+0.5)*ay, (double(k)+1.0)*az);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
   }
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new FCCGenerator(args); }
void destroy(Module * m) { delete m; }


