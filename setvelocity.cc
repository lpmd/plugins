//
//
//

#include "setvelocity.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

SetVelocityModifier::SetVelocityModifier(std::string args): Plugin("setvelocity", "2.1")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = setvelocity                                              \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to set a specific velocity to a set of atoms. Note   \n";
 std::cout << " the difference between add (addvelocity) and set (setvelocity) of velocity.   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      velocity      : Velocity vector that's going to be setted to each atom.  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Calling the module in a control file :                                        \n";
 std::cout << " use setvelocity                                                               \n";
 std::cout << "     velocity <0.002,0.001,0.005>                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Applying the plugin :                                                         \n";
 std::cout << " apply setvelocity start=0 each=10 end=100                                     \n";
 std::cout << "      With this we set the atom velocity between steps 0 and 100 each 10.      \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void SetVelocityModifier::Apply(Simulation & sim)
{
 BasicParticleSet & atoms = sim.Atoms();
 DebugStream() << "-> Setting velocity " << velocity << " to " << atoms.Size() << " atoms.\n";  
 for (long int i=0;i<atoms.Size();++i) atoms[i].Velocity() = velocity;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SetVelocityModifier(args); }
void destroy(Plugin * m) { delete m; }

