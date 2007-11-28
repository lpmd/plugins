//
//
//

#include <lpmd/util.h>

#include "tempscaling.h"

using namespace lpmd;

TempScalingModifier::TempScalingModifier(std::string args): Module("tempscaling")
{
 fromtemp = 300.0;
 totemp = 300.0;
 ProcessArguments(args);
}

TempScalingModifier::~TempScalingModifier() { }

void TempScalingModifier::SetParameter(std::string name)
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

void TempScalingModifier::Show() const
{
 Module::Show();
 std::cout << "   from  = " << fromtemp << '\n';
 std::cout << "   to    = " << totemp << '\n';
 std::cout << "   start = " << start_step << '\n';
 std::cout << "   end   = " << end_step << '\n';
 std::cout << "   each  = " << interval << '\n';
}

std::string TempScalingModifier::Keywords() const
{
 return "from to start end each";
}

void TempScalingModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 double set_temp = LeverRule(md.CurrentStep(), start_step, end_step, fromtemp, totemp);
 std::cerr << "-> Rescaling temperature to T = " << set_temp << '\n';  
 sc.SetTemperature(set_temp);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new TempScalingModifier(args); }
void destroy(Module * m) { delete m; }


