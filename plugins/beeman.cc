//
//
//

#include "beeman.h"

using namespace lpmd;

Beeman::Beeman(std::string args): Module("beeman")
{
 start_step = 1;
 ProcessArguments(args);
}

Beeman::~Beeman() { }

void Beeman::SetParameter(std::string name)
{
 if (name == "dt")
 {
  (*this)["dt"] = GetNextWord();
  dt = GetDouble("dt");
 }
 if (name == "start")
 {
  AssignParameter("start", GetNextWord());
  start_step = GetInteger("start");
 }
}

void Beeman::Show() const
{
 Module::Show();
 std::cout << "   dt = " << dt << '\n';
 std::cout << "start = " << start_step << '\n';
}

std::string Beeman::Keywords() const { return "dt start"; }

void Beeman::Initialize(SimulationCell & sc, Potential & p)
{
 for (long i=0;i<sc.Size();++i) auxlist.push_back(Vector(0.0));
 UseOldCell(sc);
 SimulationCell & oldsc = OldCell();
 p.UpdateForces(oldsc);
}

void Beeman::AdvancePositions(SimulationCell & sc)
{
 Vector newpos;
 SimulationCell & oldsc = OldCell();
 for (long i=0;i<sc.Size();++i)
 {
  const Atom & now = sc.GetAtom(i);
  if ((now.IsTypeSet()) && (now.Type().GetBool("fixedpos") == true)) continue;
  else
  {
   const Atom & old = oldsc.GetAtom(i);
   auxlist[i] = old.Acceleration();
   newpos = now.Position() + now.Velocity()*dt + (2.0/3.0)*now.Acceleration()*dt*dt - (1.0/6.0)*old.Acceleration()*dt*dt;
   oldsc.SetAcceleration(i, now.Acceleration());
   sc.SetPosition(i, newpos);
  }
 }
}

void Beeman::AdvanceVelocities(SimulationCell & sc)
{
 Vector newvel;
 SimulationCell & oldsc = OldCell();
 for (long i=0;i<sc.Size();++i)
 {
  const Atom & now = sc.GetAtom(i);
  const Atom & old = oldsc.GetAtom(i);
  newvel = now.Velocity() + (1.0/3.0)*now.Acceleration()*dt + (5.0/6.0)*old.Acceleration()*dt - (1.0/6.0)*auxlist[i]*dt;
  sc.SetVelocity(i, newvel);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Beeman(args); }
void destroy(Module * m) { delete m; }



