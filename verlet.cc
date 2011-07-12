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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = verlet                                                   \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin is used to itnegrate the movement equation using the verlet   \n";
 std::cout < "  method.                                                                       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Set the integration time in femto-seconds.               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use verlet                                                                    \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";
 std::cout << " integrator verlet start=1000                                                  \n";
 std::cout << "      The plugin can be called at the begin (without start option) or during   \n";
 std::cout << " the simulation at any other time step, with this you can change the           \n";
 std::cout << " integration plugin during the simulation.                                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

}

void Verlet::Initialize(Simulation & sim, Potential & p) { assert(&p != 0); UseOldConfig(sim); }

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


