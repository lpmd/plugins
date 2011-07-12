//
//
//

#include "tempscaling.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

TempScalingModifier::TempScalingModifier(std::string args): Plugin("tempscaling", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("from", "300.0");
 DefineKeyword("to", "300.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 fromtemp = double(params["from"]);
 totemp = double(params["to"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

TempScalingModifier::~TempScalingModifier() { }

void TempScalingModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = tempscaling                                              \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      this plugin is used to scale the temperature of the system using the     \n";
 std::cout << " velocity rescaling on the atoms of the cell.                                  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      from          : Initial temperature.                                     \n";
 std::cout << "      to            : Final required temperature.                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use tempscaling                                                               \n";
 std::cout << "     from 84.0                                                                 \n";
 std::cout << "     to   10.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";
 std::cout << " apply tempscaling start=0 each=10 end=100                                     \n";
 std::cout << "      With this we reduce the temperature of the simulation from 84 to 10 K.   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void TempScalingModifier::Apply(Simulation & sim)
{
 double set_temp = ValueAtStep(sim.CurrentStep(), fromtemp, totemp);
 DebugStream() << "-> Rescaling temperature to T = " << set_temp << '\n';  
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 sim.SetTemperature(set_temp,(atoms.HaveAny(Tag("fixedvel"))||atoms.HaveAny(Tag("fixedpos"))));
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new TempScalingModifier(args); }
void destroy(Plugin * m) { delete m; }

