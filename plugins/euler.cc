//
//
//

#include <iostream>

#include "euler.h"

#ifdef __ICC
#pragma warning (disable:869)
#endif

using namespace lpmd;

Euler::Euler(std::string args): Module("euler")
{
 start_step = 1;
 ProcessArguments(args);
}

Euler::~Euler() { }

void Euler::SetParameter(std::string name)
{
 if (name == "dt")
 {
  AssignParameter("dt", GetNextWord());
  dt = GetDouble("dt");
 }
 if (name == "start")
 {
  AssignParameter("start", GetNextWord());
  start_step = GetInteger("start");
 }
}

void Euler::Show() const
{
 Module::Show();
 std::cout << "   dt = " << dt << '\n';
 std::cout << "start = " << start_step << '\n';
}

std::string Euler::Keywords() const { return "dt start"; }

void Euler::Advance(SimulationCell & sc)
{
 Atom now;
 Vector newpos, newvel;
 for (int i=0;i<sc.Size();++i)
 {
  now = sc.GetAtom(i);
  if ((now.IsTypeSet()) && (now.Type().GetBool("fixedpos") == true)) continue;
  else
  {
   newpos = now.Position() + now.Velocity()*dt;
   newvel = now.Velocity() + now.Acceleration()*dt;
   sc.SetPosition(i, newpos);
   sc.SetVelocity(i, newvel);
  }
 }
}

#ifdef __ICC
#pragma warning (default:869)
#endif

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Euler(args); }
void destroy(Module * m) { delete m; }

