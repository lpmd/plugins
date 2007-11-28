//
//
//

#include "verlet.h"

using namespace lpmd;

Verlet::Verlet(std::string args): Module("verlet")
{
 start_step = 1;
 ProcessArguments(args);
}

Verlet::~Verlet() { }

void Verlet::SetParameter(std::string name)
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

void Verlet::Show() const
{
 Module::Show();
 std::cout << "   dt = " << dt << '\n';
 std::cout << "start = " << start_step << '\n';
}

std::string Verlet::Keywords() const { return "dt start"; }

void Verlet::Initialize(SimulationCell & sc, Potential & p)
{
 UseOldCell(sc);
}

void Verlet::Advance(SimulationCell & sc)
{
 Vector oldpos, newpos, newvel;
 SimulationCell & oldsc = OldCell();
 for (long i=0;i<sc.Size();++i)
 {
  const Atom & now = sc.GetAtom(i);
  if ((now.IsTypeSet()) && (now.Type().GetBool("fixedpos") == true)) continue;
  else
  {
   oldpos = oldsc[i].Position();
   newpos = 2.0*now.Position() - oldpos + now.Acceleration()*dt*dt;

   oldsc.SetPosition(i, now.Position());
   oldsc.SetVelocity(i, now.Velocity());
   sc.SetPosition(i, newpos);

   newvel = sc.Displacement(oldpos, sc.GetAtom(i).Position())/(2.0*dt);

   sc.SetVelocity(i, newvel);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Verlet(args); }
void destroy(Module * m) { delete m; }


