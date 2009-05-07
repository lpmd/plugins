//
//
//

#include "verlet.h"

#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

Verlet::Verlet(std::string args): Module("verlet")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("dt", "1.0");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = GetDouble("dt");
 start = GetInteger("start");
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

void Verlet::Initialize(SimulationCell & sc, Potential & p) { UseOldCell(sc); }

void Verlet::Advance(SimulationCell & sc, long i)
{
 SimulationCell & oldsc = OldCell();
 const Atom & now = sc[i];
 Vector oldpos = oldsc[i].Position();
 Vector newpos = 2.0*now.Position() - oldpos + now.Acceleration()*dt*dt;
 oldsc.SetPosition(i, now.Position());
 oldsc.SetVelocity(i, now.Velocity());
 sc.SetPosition(i, newpos);
 const Vector newvel = sc.Displacement(oldpos, sc[i].Position())/(2.0*dt);
 sc.SetVelocity(i, newvel);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Verlet(args); }
void destroy(Module * m) { delete m; }


