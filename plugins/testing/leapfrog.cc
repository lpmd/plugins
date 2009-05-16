//
//
//

#include "leapfrog.h"

#include <lpmd/simulation.h>
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

void Leapfrog::Initialize(Simulation & sim, Potential & p)
{
 UseOldConfig(sim);
 BasicParticleSet & oldatoms = OldConfig().Atoms();
 p.UpdateForces(OldConfig());
 // Necesitamos las velocidades al tiempo -0.5*dt, no -dt
 for (long int i=0;i<oldatoms.Size();++i)
 {
  const Atom & old = oldatoms[i];
  oldatoms[i].Velocity() += old.Acceleration()*0.5*dt;
 } 
}

void Leapfrog::Advance(Simulation & sim, long i)
{
 BasicParticleSet & atoms = sim.Atoms();
 BasicCell & cell = sim.Cell();
 BasicParticleSet & oldatoms = OldConfig().Atoms();
 const Atom & now = atoms[i];                      // was Atom now = sc[i];
 Vector vhalf = oldatoms[i].Velocity() + now.Acceleration()*dt;
 atoms[i].Position() = cell.FittedInside(atoms[i].Position() + vhalf*dt);
 atoms[i].Velocity() = 0.5*(vhalf+oldatoms[i].Velocity());
 oldatoms[i].Velocity() = vhalf;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Leapfrog(args); }
void destroy(Module * m) { delete m; }


