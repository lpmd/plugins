//
//
//

#include "nosehoover.h"

#include <lpmd/simulationcell.h>
#include <lpmd/potential.h>
#include <lpmd/physunits.h>

using namespace lpmd;

NoseHoover::NoseHoover(std::string args): Module("nosehoover")
{
 AssignParameter("dt", "1.0");
 AssignParameter("start", "1");
 AssignParameter("fmass", "1000000.0");
 AssignParameter("t", "300.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = GetDouble("dt");
 start_step = GetInteger("start");
 q = GetDouble("fmass");
 temp = GetDouble("t");
}

NoseHoover::~NoseHoover() { }

void NoseHoover::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = nosehoover                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para simular ensemble NVT.                        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n";
 std::cout << "                      integrador.                                              \n";
 std::cout << "      fmass         : Masa ficticia Q                                          \n";
 std::cout << "      t             : Temperatura fijada por el integrador                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string NoseHoover::Keywords() const { return "dt fmass t start"; }

void NoseHoover::Initialize(SimulationCell & sc, Potential & p)
{
 for (long i=0;i<sc.Size();++i) auxlist.push_back(Vector());
 UseOldCell(sc);
 SimulationCell & oldsc = OldCell();
 p.UpdateForces(oldsc);
 friction = 0.0;
}

void NoseHoover::AdvancePosition(SimulationCell & sc, long i)
{
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 const Atom & old = oldsc[i];
 auxlist[i] = old.Acceleration();
 friction += ((3.0*sc.Size()/q)*KBOLTZMANN*(sc.Temperature()-temp)*dt);
 const Vector newacc = now.Acceleration() - now.Velocity()*friction;
 sc.SetAcceleration(i, newacc);
 Vector newpos = now.Position() + now.Velocity()*dt + (2.0/3.0)*newacc*dt*dt - (1.0/6.0)*old.Acceleration()*dt*dt;
 oldsc.SetAcceleration(i, newacc);
 sc.SetPosition(i, newpos);
}

void NoseHoover::AdvanceVelocity(SimulationCell & sc, long i)
{
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 const Atom & old = oldsc[i];
 sc.SetVelocity(i, now.Velocity() + (1.0/3.0)*now.Acceleration()*dt + (5.0/6.0)*old.Acceleration()*dt - (1.0/6.0)*auxlist[i]*dt);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NoseHoover(args); }
void destroy(Module * m) { delete m; }

