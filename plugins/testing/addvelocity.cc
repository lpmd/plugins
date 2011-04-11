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
 DefineKeyword("debug","none");
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
 std::cout << "      The plugin is used to add a specific velocity to a set of atoms that     \n";
 std::cout << "      have the tag extravel setting in true.                                   \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Loading the plugin :                                                          \n";
 std::cout << " use addvelocity                                                               \n";
 std::cout << "     velocity <0.002,0.001,0.005>                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Apply the plugin :                                                            \n";
 std::cout << " apply extravel start=0 each=10 end=100                                      \n\n";
 std::cout << "      With this we apply a extra-velocity between steps 0 and 100 each 10,     \n";
 std::cout << "      to all atoms that have the tag extravel in true.                         \n";
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

