//
//
//

#include "beeman.h"

#include <lpmd/simulation.h>
#include <lpmd/potential.h>

using namespace lpmd;

Beeman::Beeman(std::string args): Plugin("beeman", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = params["dt"];
 start = int(params["start"]);
}

Beeman::~Beeman() { }

void Beeman::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando el metodo de beeman.     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n" ;
 std::cout << "                      integrador.                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use beeman                                                                    \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator beeman start=1000                                                \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void Beeman::Initialize(Simulation & sim, Potential & p)
{
 BasicParticleSet & atoms = sim.Atoms();
 for (long int i=0;i<atoms.Size();++i) auxlist.push_back(Vector());
 UseOldConfig(sim);
 p.UpdateForces(OldConfig());
}

void Beeman::AdvancePosition(Simulation & sim, long i)
{
 BasicParticleSet & atoms = sim.Atoms();
 BasicParticleSet & oldatoms = OldConfig().Atoms();
 BasicCell & cell = sim.Cell();
 const Atom & now = atoms[i];
 const Atom & old = oldatoms[i];
 auxlist[i] = old.Acceleration();
 Vector newpos = now.Position() + now.Velocity()*dt + (2.0/3.0)*now.Acceleration()*dt*dt - (1.0/6.0)*old.Acceleration()*dt*dt;
 oldatoms[i].Acceleration() = now.Acceleration();
 atoms[i].Position() = cell.FittedInside(newpos);
}

void Beeman::AdvanceVelocity(Simulation & sim, long i)
{
 BasicParticleSet & atoms = sim.Atoms();
 BasicParticleSet & oldatoms = OldConfig().Atoms();
 const Atom & now = atoms[i];
 const Atom & old = oldatoms[i];
 atoms[i].Velocity() = now.Velocity()+(1.0/3.0)*now.Acceleration()*dt+(5.0/6.0)*old.Acceleration()*dt-(1.0/6.0)*auxlist[i]*dt;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Beeman(args); }
void destroy(Plugin * m) { delete m; }

