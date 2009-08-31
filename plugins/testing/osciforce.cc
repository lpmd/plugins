//
//
//

#include "osciforce.h"

#include <lpmd/simulation.h>
#include <lpmd/session.h>

#include <iostream>

using namespace lpmd;

OsciForcePotential::OsciForcePotential(std::string args): Plugin("osciforce", "1.0")
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
 std::cout << "      El modulo setea una fuerza variable, es utilizado para trabajo con un    \n";
 std::cout << " set de particulas que se quieran configurar con ese tipo de fuerzas.          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      force        : Especifica el vector fuerza máximo en torno al cuál se    \n";
 std::cout << "                     oscila.                                                   \n";
 std::cout << "      phase        : Especifica una fase adicional (default=0)                 \n";
 std::cout << "      n            : número de subintervalos para la evaluación                \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use osciforce as sinforce                                                     \n";
 std::cout << "     force <0.0,0.0,-9.8>                                                      \n";
 std::cout << "     n 50                                                                      \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " potential sinforce Ar Ar                                                    \n\n";
 std::cout << "      De esta forma aplicamos una fuerza oscilatoria a la interaccion de los   \n";
 std::cout << " atomos de Argon.                                                              \n";
}

double OsciForcePotential::energy(Configuration & con) {assert(&con != 0); return 0.0; }//icc869

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
 counter++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new OsciForcePotential(args); }
void destroy(Plugin * m) { delete m; }

