//
//
//

#include "nullpotential.h"

#include <iostream>

using namespace lpmd;

NullPotential::NullPotential(std::string args): Module("nullpotential") { ProcessArguments(args); }

NullPotential::~NullPotential() { }

void NullPotential::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = nullpotential                                            \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo activa un potencial nulo.                                      \n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use nullpotential                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential nullpotential Ar Kr                                               \n\n";
 std::cout << "      De esta forma activa el potencial nulo entre Ar y Kr.                    \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string NullPotential::Keywords() const { return ""; }

double NullPotential::energy(SimulationCell & sc) { return 0.0; }

void NullPotential::UpdateForces(SimulationCell & sc) { }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullPotential(args); }
void destroy(Module * m) { delete m; }


