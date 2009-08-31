//
//
//

#include "hardspheres.h"

#include <lpmd/simulation.h>

#include <iostream>

using namespace lpmd;

HardSpheres::HardSpheres(std::string args): Plugin("verlet", "2.0")
{
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int((*this)["start"]);
}

HardSpheres::~HardSpheres() { }

void HardSpheres::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando el metodo de esferas     \n";
 std::cout << "      duras, este metodo no requiere de potenciales para su utilizacion.       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol        : Especifica el radio del atomo asociado que sera          \n" ;
 std::cout << "                      utilizado en la simulacion.                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use hardspheres                                                               \n";
 std::cout << "     Ar 5.0                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator hardspheres start=1000                                           \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void HardSpheres::Initialize(Simulation & sim, Potential & p) { assert(&p != 0); UseOldConfig(sim); }//icc 869

void HardSpheres::Advance(Simulation & sim, long i)
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
Plugin * create(std::string args) { return new HardSpheres(args); }
void destroy(Plugin * m) { delete m; }


