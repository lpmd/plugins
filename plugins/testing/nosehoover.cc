//
//
//

#include "nosehoover.h"

#include <lpmd/simulation.h>
#include <lpmd/potential.h>
#include <lpmd/session.h>

using namespace lpmd;

NoseHoover::NoseHoover(std::string args): Module("nosehoover")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
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

void NoseHoover::Initialize(Simulation & sc, Potential & p)
{
 for (long int i=0;i<sc.Atoms().Size();++i) auxlist.push_back(Vector());
#warning llamado a UseOldCell, que no se usa más?
 /*
 UseOldCell(sc);
 SimulationCell & oldsc = OldCell();
 p.UpdateForces(oldsc);
 friction = 0.0;
 */
}

void NoseHoover::AdvancePosition(Simulation & sim, long i)
{
 const double kboltzmann = double(Parameter(sim.GetTag(sim,"kboltzmann")));
#warning de nuevo uso de OldCell ... confusion!!
 /*
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 const Atom & old = oldsc[i];
 auxlist[i] = old.Acceleration();
 friction += ((3.0*sc.size()/q)*kboltzmann*(sc.Temperature()-temp)*dt);
 const Vector newacc = now.Acceleration() - now.Velocity()*friction;
 sc.SetAcceleration(i, newacc);
 Vector newpos = now.Position() + now.Velocity()*dt + (2.0/3.0)*newacc*dt*dt - (1.0/6.0)*old.Acceleration()*dt*dt;
 oldsc.SetAcceleration(i, newacc);
 sc.SetPosition(i, newpos);
 */
}

void NoseHoover::AdvanceVelocity(Simulation & sim, long i)
{
#warning nuevamente uso de OldCell ... more confusion!!!!
 /*
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 const Atom & old = oldsc[i];
 sc.SetVelocity(i, now.Velocity() + (1.0/3.0)*now.Acceleration()*dt + (5.0/6.0)*old.Acceleration()*dt - (1.0/6.0)*auxlist[i]*dt);
 */
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NoseHoover(args); }
void destroy(Module * m) { delete m; }
