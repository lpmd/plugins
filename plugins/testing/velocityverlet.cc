//
//
//

#include "velocityverlet.h"

#include <lpmd/simulation.h>

#include <iostream>

using namespace lpmd;

VelocityVerlet::VelocityVerlet(std::string args): Plugin("velocityverlet", "2.0")
{
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = double((*this)["dt"]);
 start = int((*this)["start"]);
}

VelocityVerlet::~VelocityVerlet() { }

void VelocityVerlet::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando el metodo velocityverlet.\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n" ;
 std::cout << "                      integrador.                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use velocityverlet                                                            \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator velocityverlet start=1000                                        \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void VelocityVerlet::Initialize(Simulation & sim, Potential & p) { UseOldConfig(sim); }

void VelocityVerlet::AdvancePosition(Simulation & sim, long i)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms(); 
 lpmd::BasicCell & cell = sim.Cell();
 lpmd::BasicParticleSet & oldatoms = OldConfig().Atoms(); 
 Vector newpos = atoms[i].Position() + atoms[i].Velocity()*dt + 0.5*atoms[i].Acceleration()*dt*dt;
 oldatoms[i].Acceleration() = atoms[i].Acceleration();
 atoms[i].Position() = cell.FittedInside(newpos);
}

void VelocityVerlet::AdvanceVelocity(Simulation & sim, long i)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicParticleSet & oldatoms = OldConfig().Atoms();
 atoms[i].Velocity() = atoms[i].Velocity() + 0.5*dt*(oldatoms[i].Acceleration() + atoms[i].Acceleration());
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new VelocityVerlet(args); }
void destroy(Plugin * m) { delete m; }



