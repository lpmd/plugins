//
//
//

#include "nullintegrator.h"

#include <lpmd/potential.h>

#include <iostream>

using namespace lpmd;

NullIntegrator::NullIntegrator(std::string args): Module("nullintegrator")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start_step = GetInteger("start");
}

NullIntegrator::~NullIntegrator() { }

void NullIntegrator::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo realia integraciones fantasmas.                                \n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use nullintegrator                                                            \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " integrator nullintegrator start=15000                                       \n\n";
 std::cout << "      De esta forma activa el integrador nulo a partir del step 150000.        \n";
}

void NullIntegrator::Advance(SimulationCell & sc, Potential & p) { p.UpdateForces(sc); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullIntegrator(args); }
void destroy(Module * m) { delete m; }


