//
//
//

#include "setvelocity.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

SetVelocityModifier::SetVelocityModifier(std::string args): Module("setvelocity")
{
 ParamList & params = (*this);
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
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
Module * create(std::string args) { return new SetVelocityModifier(args); }
void destroy(Module * m) { delete m; }

