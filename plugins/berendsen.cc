//
//
//

#include <lpmd/util.h>

#include "berendsen.h"

using namespace lpmd;

BerendsenModifier::BerendsenModifier(std::string args): Module("berendsen")
{
 fromtemp = 300.0;
 totemp = 300.0;
 tau = 400.0;
 stop_thermostat = -1;
 ProcessArguments(args);
}

BerendsenModifier::~BerendsenModifier() { }

void BerendsenModifier::SetParameter(std::string name)
{
 if (name == "from")
 {
  AssignParameter("from", GetNextWord());
  fromtemp = GetDouble("from");
 }
 if (name == "to")
 {
  AssignParameter("to", GetNextWord());
  totemp = GetDouble("to");
 }
 if (name == "tau")
 {
  AssignParameter("tau", GetNextWord());
  tau = GetDouble("tau");
 }
 if (name == "start")
 {
  AssignParameter("start", GetNextWord());
  start_step = GetInteger("start");
 }
 if (name == "end")
 {
  AssignParameter("end", GetNextWord());
  end_step = GetInteger("end");
 }
 if (name == "each")
 {
  AssignParameter("each", GetNextWord());
  interval = GetInteger("each");
 }
}

void BerendsenModifier::Show() const
{
 Module::Show();
 std::cout << "   from   = " << fromtemp << '\n';
 std::cout << "   to     = " << totemp << '\n';
 std::cout << "   tau    = " << tau << '\n';
 std::cout << "   start  = " << start_step << '\n';
 std::cout << "   end    = " << end_step << '\n';
 std::cout << "   each   = " << interval << '\n';
}

std::string BerendsenModifier::Keywords() const 
{
 return "from to tau start end each";
}

void BerendsenModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 double timestep = md.GetIntegrator().Timestep();
 double set_temp = LeverRule(md.CurrentStep(), start_step, end_step, fromtemp, totemp);
 if (stop_thermostat == -1)
 {
  stop_thermostat = long(double(md.CurrentStep()) + tau/timestep);
  std::cerr << "-> Berendsen thermostat activated (stops in step " << stop_thermostat << "), rescaling temperature to T = " << set_temp << '\n';
  old_step = interval;
  interval = 1;
 }
 else
 {
  if (md.CurrentStep() <= stop_thermostat) 
  {
   sc.SetTemperature(set_temp, timestep, tau);
  }
  else
  {
   stop_thermostat = -1;
   std::cerr << "-> Berendsen thermostat stopped." << '\n';
   interval = old_step;
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new BerendsenModifier(args); }
void destroy(Module * m) { delete m; }


