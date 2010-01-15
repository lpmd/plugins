//
//
//

#include "constantforce.h"

#include <lpmd/simulation.h>
#include <lpmd/session.h>

#include <iostream>

using namespace lpmd;

ConstantForcePotential::ConstantForcePotential(std::string args): Plugin("constantforce", "2.0")
{ 
 //
 DefineKeyword("force");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 force = Vector((*this)["force"].c_str());
}

ConstantForcePotential::~ConstantForcePotential() { }

void ConstantForcePotential::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo setea una fuerza constante, es utilizado para trabajo con un   \n";
 std::cout << " set de particulas que se quieran configurar con cierta fuerza.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      forcevector   : Seguido por 3 numeros reales que indican un vector para  \n";
 std::cout << "                      la fuerza.                                               \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use constantforce as gravity                                                  \n";
 std::cout << "     force <0.0,0.0,-9.8>                                                      \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " potential gravity Ar Ar                                                     \n\n";
 std::cout << "      De esta forma aplicamos una fuerza constante a la interaccion de los     \n";
 std::cout << " atomos de Argon.                                                              \n";
}

double ConstantForcePotential::energy(Configuration & con) {assert(&con != 0); return 0.0; }//icc 869

double ConstantForcePotential::AtomEnergy(lpmd::Configuration & con, long i) { return 0.0; }

void ConstantForcePotential::UpdateForces(Configuration & con) 
{ 
 lpmd::BasicParticleSet & atoms = con.Atoms();
 const double forcefactor = double(GlobalSession["forcefactor"]);
 long count = 0;
 for (long int i=0;i<atoms.Size();++i)
 {
  if ((atoms.Have(atoms[i], Tag("noconstantforce"))) && (atoms.GetTag(atoms[i], Tag("noconstantforce")) == "true")) continue;
  double mi = atoms[i].Mass();
  atoms[i].Acceleration() = atoms[i].Acceleration() + force*(forcefactor/mi);
  count++;
 }
 DebugStream() << "-> Applied constant force to " << count << " atoms" << '\n';
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new ConstantForcePotential(args); }
void destroy(Plugin * m) { delete m; }

