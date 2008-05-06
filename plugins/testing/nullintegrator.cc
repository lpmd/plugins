//
//
//

#include "nullintegrator.h"

#include <lpmd/potential.h>

#include <iostream>

using namespace lpmd;

NullIntegrator::NullIntegrator(std::string args): Module("nullintegrator")
{
 AssignParameter("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start_step = GetInteger("start");
}

NullIntegrator::~NullIntegrator() { }

void NullIntegrator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = nullintegrator                                           \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo realia integraciones fantasmas.                                \n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use nullintegrator                                                            \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " integrator nullintegrator start=15000                                       \n\n";
 std::cout << "      De esta forma activa el integrador nulo a partir del step 150000.        \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string NullIntegrator::Keywords() const { return "start"; }

void NullIntegrator::Advance(SimulationCell & sc, Potential & p) { p.UpdateForces(sc); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullIntegrator(args); }
void destroy(Module * m) { delete m; }


