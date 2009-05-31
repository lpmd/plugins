//
//
//

#include "crystal3d.h"

#include <lpmd/atom.h>
#include <lpmd/configuration.h>

using namespace lpmd;

CrystalGenerator::CrystalGenerator(std::string args): Module("crystal3d")
{
 ParamList & params = (*this);
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("nx", "1");
 DefineKeyword("ny", "1");
 DefineKeyword("nz", "1");
 DefineKeyword("symbol", "H");
 DefineKeyword("type", "fcc");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 spc = ElemNum((*this)["symbol"]);
 type = (*this)["type"];
 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);
}

CrystalGenerator::~CrystalGenerator() { }

void CrystalGenerator::ShowHelp() const
{
 std::cout << " General Info      >>                                                              \n";
 std::cout << "      El modulo es utilizado para crear celdas del tipo BCC (base-centered cubic), \n";
 std::cout << "      FCC (face-centered cubic), HCP (Hexagonal close-packed) y SC (simple cubic). \n";
 std::cout << " General Options   >>                                                              \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.        \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                              \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                              \n";
 std::cout << "      nz            : Repeticiones en la direccion Z.                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                           \n";
 std::cout << " Utilizando el Modulo :                                                            \n";
 std::cout << " input crystal3d symbol=Ar type=fcc nx=3 ny=3 nz=3                                 \n";
 std::cout << " input crystal3d Ar fcc 3 3 3                                                      \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo FCC en la simulacion.     \n";
}

std::string CrystalGenerator::Keywords() const
{
 return "symbol type nx ny nz";
}

void CrystalGenerator::Generate(Configuration & conf) const
{
 long int cc = 0;
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 double az = 1.0/double(nz);
 bool create_atoms = (atoms.Size() == 0);
 
 if (type=="bcc")
 {
  for (long k=0;k<nz;++k)
   for (long j=0;j<ny;++j)
    for (long i=0;i<nx;++i)
    {
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector(double(i)*ax, double(j)*ay, double(k)*az));
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+0.5)*ax, (double(j)+0.5)*ay, (double(k)+0.5)*az));
    }
 }
 else if (type=="fcc")
 {
  for (long k=0;k<nz;++k)
   for (long j=0;j<ny;++j)
    for (long i=0;i<nx;++i)
    {
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+0.51)*ax, (double(j)+0.5)*ay, (double(k)+0.5)*az));
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+0.51)*ax, double(j)*ay, (double(k)+1.0)*az));
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+1.01)*ax, double(j)*ay, (double(k)+0.5)*az));
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+1.01)*ax, (double(j)+0.5)*ay, (double(k)+1.0)*az));
    }
 }
 else if (type=="hcp")
 {
  for (long k=0;k<nz;++k)
   for (long j=0;j<ny;++j)
    for (long i=0;i<nx;++i)
    {
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+(2.0/3.0))*ax, (double(j)+(1.0/3.0))*ay, (double(k)+(3.0/4.0))*az));
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector((double(i)+(1.0/3.0))*ax, (double(j)+(2.0/3.0))*ay, (double(k)+(1.0/4.0))*az));
    }
 }
 else if (type=="sc")
 {
  for (long k=0;k<nz;++k)
   for (long j=0;j<ny;++j)
    for (long i=0;i<nx;++i)
    {
     if (create_atoms) atoms.Append(Atom(spc));
     atoms[cc++].Position() = cell.ScaleByCell(Vector(double(i)*ax, double(j)*ay, double(k)*az));
    }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new CrystalGenerator(args); }
void destroy(Module * m) { delete m; }


