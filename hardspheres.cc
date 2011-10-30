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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = hardspheres                                              \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin is used to integrate the equations of movement using the hard \n";
 std::cout << "      spheres method. This method does not require potentials to be used.      \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      symbol        : Sets the atom's radius.                                  \n";
 std::cout << "      start         : Determines in which step the plugin begins to be applied.\n";
 std::cout << "      end           : Determines in which step the plugin ceases to be applied.\n";
 std::cout << "      each          : Determines how often (each how many time-steps) the      \n";
 std::cout << "                      plugin must be applied.                                  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use hardspheres                                                               \n";
 std::cout << "     Ar 5.0                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";
 std::cout << " integrator hardspheres start=1000                                           \n\n";
 std::cout << "      The plugin can be called at the beginning of the simulation (without the \n";
 std::cout << "      start option or setting start=0) or at any other time step (like         \n";
 std::cout << "      start=1000). This allows you to change integration method during the     \n";
 std::cout << "      simulation.                                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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


