//
//
//

#include "euler.h"

#include <lpmd/simulationcell.h>

#include <iostream>

#ifdef __ICC
#pragma warning (disable:869)
#endif

using namespace lpmd;

Euler::Euler(std::string args): Module("euler")
{
 AssignParameter("dt", "1.0");
 AssignParameter("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = GetDouble("dt");
 start_step = GetInteger("start");
}

Euler::~Euler() { }

void Euler::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = euler                                                    \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando el metodo de euler.      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n" ;
 std::cout << "                      integrador.                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use euler                                                                     \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator euler start=1000                                                 \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Euler::Keywords() const { return "dt start"; }

void Euler::Advance(SimulationCell & sc, long i)
{
 const Atom & now = sc[i];
 sc.SetPosition(i, now.Position() + now.Velocity()*dt);
 sc.SetVelocity(i, now.Velocity() + now.Acceleration()*dt);
}

#ifdef __ICC
#pragma warning (default:869)
#endif

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Euler(args); }
void destroy(Module * m) { delete m; }

