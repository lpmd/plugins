//
//
//

#include "leapfrog.h"

#include <lpmd/simulationcell.h>
#include <lpmd/potential.h>

#include <iostream>

using namespace lpmd;

Leapfrog::Leapfrog(std::string args): Module("leapfrog")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = GetDouble("dt");
 start = GetInteger("start");
}

Leapfrog::~Leapfrog() { }

void Leapfrog::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando el algoritmo leapfrog.   \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femtosegundos para el            \n";
 std::cout << "                      integrador.                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use leapfrog                                                                  \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator leapfrog start=1000                                                \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void Leapfrog::Initialize(SimulationCell & sc, Potential & p)
{
 UseOldCell(sc);
 SimulationCell & oldsc = OldCell();
 p.UpdateForces(oldsc);
 // Necesitamos las velocidades al tiempo -0.5*dt, no -dt
 for (unsigned long int i=0;i<oldsc.size();++i)
 {
  const Atom & old = oldsc[i];
  oldsc[i].Velocity() += old.Acceleration()*0.5*dt;
 } 
}

void Leapfrog::Advance(SimulationCell & sc, long i)
{
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];                      // was Atom now = sc[i];
 Vector vhalf = oldsc[i].Velocity() + now.Acceleration()*dt;
 sc[i].Position() += vhalf*dt;
 sc[i].Velocity() = 0.5*(vhalf+oldsc[i].Velocity());
 oldsc[i].Velocity() = vhalf;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Leapfrog(args); }
void destroy(Module * m) { delete m; }


