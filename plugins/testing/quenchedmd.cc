//
//
//

#include "quenchedmd.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

QuenchedMDModifier::QuenchedMDModifier(std::string args): Module("quenchedmd")
{
 ProcessArguments(args);
 start = int((*this)["start"]);
 end = int((*this)["end"]);
 each = int((*this)["each"]);
}

QuenchedMDModifier::~QuenchedMDModifier() { }

void QuenchedMDModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = quenchedmd                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use quenchedmd                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply quenchedmd start=0 each=1 end=1000                                    \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string QuenchedMDModifier::Keywords() const
{
 return "start end each";
}

void QuenchedMDModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 //
 for (unsigned long int i=0;i<atoms.Size();++i)
 {
  if (Dot(atoms[i].Acceleration(), atoms[i].Velocity()) < 0.0) atoms[i].Velocity() = Vector(0.0, 0.0, 0.0);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new QuenchedMDModifier(args); }
void destroy(Module * m) { delete m; }

