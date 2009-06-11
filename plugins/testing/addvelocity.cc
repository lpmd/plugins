//
//
//

#include "addvelocity.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

ExtraVelModifier::ExtraVelModifier(std::string args): Plugin("extravel", "1.0")
{ 
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("velocity", "<0.0,0.0,0.0>");
 ProcessArguments(args); 
 start = int((*this)["start"]);
 end = int((*this)["end"]);
 each = int((*this)["each"]);
 velocity = Vector((*this)["velocity"].c_str());
}

ExtraVelModifier::~ExtraVelModifier() { }

void ExtraVelModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para ...                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use extravel                                                                  \n";
 std::cout << "     velocity <0.002,0.001,0.005>                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply extravel start=0 each=10 end=100                                      \n\n";
 std::cout << "      De esta forma aplicamos extravel entre 0 y 100 cada 10 steps.            \n";
}

void ExtraVelModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 DebugStream() << "-> Applying extravel with velocity ";
 FormattedWrite(DebugStream(), velocity);
 DebugStream() << '\n';
 if (!atoms.HaveAny(Tag("extravel"))) return;
 for (int i=0;i<atoms.Size();++i)
 {
  Atom at = atoms[i];
  if (atoms.Have(atoms[i], Tag("extravel")) && (bool(Parameter(atoms.GetTag(atoms[i], Tag("extravel")))))) 
     atoms[i].Velocity() += velocity;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new ExtraVelModifier(args); }
void destroy(Plugin * m) { delete m; }

