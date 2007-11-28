//
//
//

#include <iostream>

#include "nullpairpotential.h"

using namespace lpmd;

NullPairPotential::NullPairPotential(std::string args): Module("nullpairpotential") 
{ 
 ProcessArguments(args);
}

double NullPairPotential::pairEnergy(const double & r) const { return 0.0; }

Vector NullPairPotential::pairForce(const Vector & r) const { return Vector(0.0); }

//
//
//
void NullPairPotential::SetParameter(std::string name) { }

void NullPairPotential::Show() const
{
 Module::Show();
}

std::string NullPairPotential::Keywords() const { return ""; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullPairPotential(args); }
void destroy(Module * m) { delete m; }


