//
//
//

#include "setvelocity.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

SetVelocityModifier::SetVelocityModifier(std::string args): Plugin("setvelocity", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("velocity");
 DefineKeyword("filterby", "none");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 velocity = Vector(params["velocity"].c_str());
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

SetVelocityModifier::~SetVelocityModifier() { }

void SetVelocityModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << '\n';
}

void SetVelocityModifier::Apply(Simulation & sim)
{
 BasicParticleSet & atoms = sim.Atoms();
 DebugStream() << "-> Applying velocity " << velocity << " to " << atoms.Size() << " atoms\n";  
 for (long int i=0;i<atoms.Size();++i) atoms[i].Velocity() = velocity;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SetVelocityModifier(args); }
void destroy(Plugin * m) { delete m; }

