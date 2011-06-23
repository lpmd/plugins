//
//
//

#include "nosehoover.h"

#include <lpmd/simulation.h>
#include <lpmd/potential.h>
#include <lpmd/session.h>
#include <lpmd/properties.h>

using namespace lpmd;

NoseHoover::NoseHoover(std::string args): Plugin("nosehoover", "2.0")
{
 //
 DefineKeyword("fmass");
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 DefineKeyword("t", "300.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = double((*this)["dt"]);
 start = int((*this)["start"]);
 q = double((*this)["fmass"]);
 temp = double((*this)["t"]);
}

NoseHoover::~NoseHoover() { }

void NoseHoover::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para simular ensemble NVT.                        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n";
 std::cout << "                      integrador.                                              \n";
 std::cout << "      fmass         : Masa ficticia Q                                          \n";
 std::cout << "      t             : Temperatura fijada por el integrador                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use nosehoover                                                                \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << "     fmass 1.0                                                                 \n";
 std::cout << "     t  600.0                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator nosehoover start=1000                                            \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void NoseHoover::Initialize(Simulation & sim, Potential & p)
{
 for (long int i=0;i<sim.Atoms().Size();++i) auxlist.push_back(Vector());
 UseOldConfig(sim);
 Configuration & oldsc = OldConfig();
 p.UpdateForces(oldsc);
 friction = 0.0;
}

void NoseHoover::AdvancePosition(Simulation & sim, long i)
{
 const double kboltzmann = double(GlobalSession["kboltzmann"]);
 
 Configuration & oldsc = OldConfig();
 BasicCell & cell = sim.Cell();
 Atom now = sim.Atoms()[i];
 Atom old = oldsc.Atoms()[i];
 auxlist[i] = old.Acceleration();
 friction += ((3.0*sim.Atoms().Size()/q)*kboltzmann*(Temperature(sim.Atoms())-temp)*dt);
 const Vector newacc = now.Acceleration() - now.Velocity()*friction;
 sim.Atoms()[i].Acceleration() = newacc;
 Vector newpos = now.Position() + now.Velocity()*dt + (2.0/3.0)*newacc*dt*dt - (1.0/6.0)*old.Acceleration()*dt*dt;
 oldsc.Atoms()[i].Acceleration() =  newacc;
 sim.Atoms()[i].Position() = cell.FittedInside(newpos);
}

void NoseHoover::AdvanceVelocity(Simulation & sim, long i)
{
 Configuration & oldsc = OldConfig();
 Atom now = sim.Atoms()[i];
 Atom old = oldsc.Atoms()[i];
 sim.Atoms()[i].Velocity() = now.Velocity() + (1.0/3.0)*now.Acceleration()*dt + (5.0/6.0)*old.Acceleration()*dt - (1.0/6.0)*auxlist[i]*dt;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new NoseHoover(args); }
void destroy(Plugin * m) { delete m; }
