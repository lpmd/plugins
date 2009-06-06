//
//
//

#include "velocityverlet.h"

#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

VelocityVerlet::VelocityVerlet(std::string args): Module("velocityverlet")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = GetDouble("dt");
 start = GetInteger("start");
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

void VelocityVerlet::Initialize(SimulationCell & sc, Potential & p) { UseOldCell(sc); }

void VelocityVerlet::AdvancePosition(SimulationCell & sc, long i)
{
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 Vector newpos = now.Position() + now.Velocity()*dt + 0.5*now.Acceleration()*dt*dt;
 oldsc[i].Acceleration() = now.Acceleration();
 sc.SetPosition(i, newpos);
}

void VelocityVerlet::AdvanceVelocity(SimulationCell & sc, long i)
{
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 const Atom & old = oldsc[i];
 sc[i].Velocity() = now.Velocity() + 0.5*dt*(old.Acceleration() + now.Acceleration());
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VelocityVerlet(args); }
void destroy(Module * m) { delete m; }



