//
//
//

#include "nullmixer.h"

#include <iostream>

using namespace lpmd;

NullMixer::NullMixer(std::string args): Module("nullmixer")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 ProcessArguments(args);
}

NullMixer::~NullMixer() { }

void NullMixer::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << "  Aplicando el Modulo                                                          \n";
 std::cout << " apply nullmixer                                                               \n";
}

SimulationCell NullMixer::Apply(SimulationCell & sc1, SimulationCell & sc2)
{
 SimulationCell mixed;
 return mixed;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullMixer(args); }
void destroy(Module * m) { delete m; }

