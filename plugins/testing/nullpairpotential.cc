//
//
//

#include "nullpairpotential.h"

#include <iostream>

using namespace lpmd;

NullPairPotential::NullPairPotential(std::string args): Plugin("nullpairpotential", "2.0")
{ 
 ProcessArguments(args); 
}

void NullPairPotential::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo activa un potencial de pares nulo.                             \n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use nullpairpotential                                                         \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential nullpairpotential Ar Kr                                           \n\n";
 std::cout << "      De esta forma activa el potencial de pares nulos entre Ar y Kr.          \n";
}

double NullPairPotential::pairEnergy(const double & r) const { return 0.0; }

Vector NullPairPotential::pairForce(const Vector & r) const { return Vector(); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new NullPairPotential(args); }
void destroy(Plugin * m) { delete m; }


