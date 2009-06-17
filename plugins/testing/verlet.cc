//
//
//

#include "verlet.h"

#include <lpmd/simulation.h>

#include <iostream>

using namespace lpmd;

Verlet::Verlet(std::string args): Plugin("verlet", "2.0")
{
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = double((*this)["dt"]);
 start = int((*this)["start"]);
}

Verlet::~Verlet() { }

void Verlet::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando el metodo de verlet.     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n" ;
 std::cout << "                      integrador.                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use verlet                                                                    \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator verlet start=1000                                                \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void Verlet::Initialize(Simulation & sim, Potential & p) { UseOldConfig(sim); }

void Verlet::Advance(Simulation & sim, long i)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicParticleSet & oldatoms = OldConfig().Atoms();
 lpmd::BasicCell & cell = sim.Cell();
 Vector oldpos = oldatoms[i].Position();
 Vector newpos = 2.0*atoms[i].Position() - oldatoms[i].Position() + atoms[i].Acceleration()*dt*dt;
 oldatoms[i].Position() = atoms[i].Position();
 oldatoms[i].Velocity() = atoms[i].Velocity();
 atoms[i].Position() = cell.FittedInside(newpos);
 Vector newvel = cell.Displacement(oldpos, atoms[i].Position())/(2.0*dt);
 atoms[i].Velocity() = newvel;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Verlet(args); }
void destroy(Plugin * m) { delete m; }


