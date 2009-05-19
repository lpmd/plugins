//
//
//

#include "crystalsc.h"

#include <lpmd/atom.h>
#include <lpmd/configuration.h>

using namespace lpmd;

SCGenerator::SCGenerator(std::string args): Module("crystalsc")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("nx", "1");
 DefineKeyword("ny", "1");
 DefineKeyword("nz", "1");
 DefineKeyword("symbol", "H");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 spc = ElemNum(GetString("symbol"));
 nx = GetInteger("nx");
 ny = GetInteger("ny");
 nz = GetInteger("nz");
}

SCGenerator::~SCGenerator() { }

void SCGenerator::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para crear celdas del tipo cubica simple (SC).    \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << "      nz            : Repeticiones en la direccion Z.                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystalfcc symbol=Ar nx=2 ny=2 nz=2                                     \n";
 std::cout << " input crystalfcc Ar 2 2 2                                                   \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo SC en la simulacion.  \n";
}


void SCGenerator::Generate(Configuration & conf) const
{
 long int cc = 0;
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 double az = 1.0/double(nz);
 for (long k=0;k<nz;++k)
  for (long j=0;j<ny;++j)
   for (long i=0;i<nx;++i)
   {
//    p = Vector(double(i)*ax, double(j)*ay, double(k)*az);
//    sc.Create(new Atom(spc));
//    sc.SetFracPosition(cc++, p);
    atoms[cc++].Position() = cell.ScaleByCell(Vector(double(i)*ax, double(j)*ay, double(k)*az));
   }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new SCGenerator(args); }
void destroy(Module * m) { delete m; }


