//
//
//

#include "osciforce.h"

#include <lpmd/simulation.h>
#include <lpmd/session.h>

#include <iostream>

using namespace lpmd;

OsciForcePotential::OsciForcePotential(std::string args): Plugin("constantforce", "2.0")
{ 
 //
 DefineKeyword("force");
 DefineKeyword("phase","0.0e0");
 DefineKeyword("n","10");
 DefineKeyword("counter","0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 force = Vector((*this)["force"].c_str());
 phase = double((*this)["phase"]);
 n = int((*this)["n"]);
 counter = int((*this)["counter"]);
}

OsciForcePotential::~OsciForcePotential() { }

void OsciForcePotential::ShowHelp() const
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

double OsciForcePotential::energy(Configuration & con) { return 0.0; }

void OsciForcePotential::UpdateForces(Configuration & con) 
{ 
 lpmd::BasicParticleSet & atoms = con.Atoms();
 const double forcefactor = double(GlobalSession["forcefactor"]);
 lpmd::Vector tmpforce = force*sin(2*M_PI*(double(counter/n)) + phase);
 for (long int i=0;i<atoms.Size();++i)
 {
  double mi = atoms[i].Mass();
  atoms[i].Acceleration() = atoms[i].Acceleration() + tmpforce*(forcefactor/mi);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new OsciForcePotential(args); }
void destroy(Plugin * m) { delete m; }

