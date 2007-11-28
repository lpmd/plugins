//
//
//

#include <iostream>

#include "nullpotential.h"

using namespace lpmd;

NullPotential::NullPotential(std::string args): Module("nullpotential") 
{ 
 ProcessArguments(args);
}

NullPotential::~NullPotential() { }

double NullPotential::energy(SimulationCell & sc) { return 0.0; }

void NullPotential::UpdateForces(SimulationCell & sc) { }

//
//
//
void NullPotential::SetParameter(std::string name) { }

void NullPotential::Show() const
{
 Module::Show();
}

std::string NullPotential::Keywords() const { return ""; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullPotential(args); }
void destroy(Module * m) { delete m; }


