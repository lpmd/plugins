//
//
//

#include "crystalbcc.h"

#include <lpmd/atom.h>
#include <lpmd/configuration.h>

using namespace lpmd;

BCCGenerator::BCCGenerator(std::string args): Module("crystalbcc")
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

BCCGenerator::~BCCGenerator() { }

void BCCGenerator::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para crear celdas del tipo cubica centrada en el  \n";
 std::cout << "      cuerpo (BCC).                                                            \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << "      nz            : Repeticiones en la direccion Z.                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystalbcc symbol=Fe nx=3 ny=3 nz=3                                     \n";
 std::cout << " input crystalbcc Fe 3 3 3                                                   \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo BCC en la simulacion. \n";
}

void BCCGenerator::Generate(Configuration & conf) const
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
//    p = Vector((double(i)+0.5)*ax, (double(j)+0.5)*ay, (double(k)+0.5)*az);
//    sc.Create(new Atom(spc));
//    sc.SetFracPosition(cc++, p);
    atoms[cc++].Position() = cell.ScaleByCell(Vector(double(i)*ax, double(j)*ay, double(k)*az));
    atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+0.5)*ax, (double(j)+0.5)*ay, (double(k)+0.5)*az));
   }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new BCCGenerator(args); }
void destroy(Module * m) { delete m; }


