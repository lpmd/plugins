//
//
//

#include "quenchedmd.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

QuenchedMDModifier::QuenchedMDModifier(std::string args): Plugin("quenchedmd", "2.0")
{
 ProcessArguments(args);
 start = int((*this)["start"]);
 end = int((*this)["end"]);
 each = int((*this)["each"]);
}

QuenchedMDModifier::~QuenchedMDModifier() { }

void QuenchedMDModifier::ShowHelp() const
{
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use quenchedmd                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply quenchedmd start=0 each=1 end=1000                                    \n\n";
}

void QuenchedMDModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 //
 for (long int i=0;i<atoms.Size();++i)
 {
  if (Dot(atoms[i].Acceleration(), atoms[i].Velocity()) < 0.0) atoms[i].Velocity() = Vector(0.0, 0.0, 0.0);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new QuenchedMDModifier(args); }
void destroy(Plugin * m) { delete m; }

