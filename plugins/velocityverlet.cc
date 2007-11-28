//
//
//

#include "velocityverlet.h"

using namespace lpmd;

VelocityVerlet::VelocityVerlet(std::string args): Module("velocityverlet")
{
 start_step = 1;
 ProcessArguments(args);
}

VelocityVerlet::~VelocityVerlet() { }

void VelocityVerlet::SetParameter(std::string name)
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

void VelocityVerlet::Show() const
{
 Module::Show();
 std::cout << "   dt = " << dt << '\n';
 std::cout << "start = " << start_step << '\n';
}

std::string VelocityVerlet::Keywords() const { return "dt start"; }

void VelocityVerlet::Initialize(SimulationCell & sc, Potential & p)
{
 UseOldCell(sc);
}

void VelocityVerlet::AdvancePositions(SimulationCell & sc)
{
 Vector newpos;
 SimulationCell & oldsc = OldCell();
 for (long i=0;i<sc.Size();++i)
 {
  const Atom & now = sc.GetAtom(i);
  if ((now.IsTypeSet()) && (now.Type().GetBool("fixedpos") == true)) continue;
  else
  {
   newpos = now.Position() + now.Velocity()*dt + 0.5*now.Acceleration()*dt*dt;
   oldsc.SetAcceleration(i, now.Acceleration());
   sc.SetPosition(i, newpos);
  }
 }
}

void VelocityVerlet::AdvanceVelocities(SimulationCell & sc)
{
 Vector newvel;
 SimulationCell & oldsc = OldCell();
 for (long i=0;i<sc.Size();++i)
 {
  const Atom & now = sc.GetAtom(i);
  const Atom & old = oldsc.GetAtom(i);
  newvel = now.Velocity() + 0.5*dt*(old.Acceleration() + now.Acceleration());
  sc.SetVelocity(i, newvel);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VelocityVerlet(args); }
void destroy(Module * m) { delete m; }



