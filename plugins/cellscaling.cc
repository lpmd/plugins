//
//
//

#include "cellscaling.h"

using namespace lpmd;

CellScalingModifier::CellScalingModifier(std::string args): Module("cellscaling")
{
 axis = -1;
 percent = 0.0;
 ProcessArguments(args);
}

CellScalingModifier::~CellScalingModifier() { }

void CellScalingModifier::SetParameter(std::string name)
{
 if (name == "percent")
 {
  AssignParameter("percent", GetNextWord());
  percent = GetDouble("percent");
 }
 if (name == "axis")
 {
  AssignParameter("axis", GetNextWord());
  std::string ax = GetString("axis");
  if ((ax == "x") || (ax == "X")) axis = 0;
  if ((ax == "y") || (ax == "Y")) axis = 1;
  if ((ax == "z") || (ax == "Z")) axis = 2;
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

void CellScalingModifier::Show() const
{
 Module::Show();
 std::cout << "   percent = " << percent << '\n';
 if (axis == -1) std::cout << "   axis  = all" << '\n';
 if (axis == 0) std::cout << "   axis  = X" << '\n';
 if (axis == 1) std::cout << "   axis  = Y" << '\n';
 if (axis == 2) std::cout << "   axis  = Z" << '\n';
 std::cout << "   start = " << start_step << '\n';
 std::cout << "   end   = " << end_step << '\n';
 std::cout << "   each  = " << interval << '\n';
}

std::string CellScalingModifier::Keywords() const
{
 return "percent axis start end each";
}

void CellScalingModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 if (axis == -1)
 {
  std::cerr << "-> Rescaling cell by " << percent << "%\n";  
  sc.RescalePercent(percent); 
 }
 else
 {
  std::cerr<< "-> Rescaling axis " << axis << " in a " << percent << "%\n";
  sc.RescalePercent(percent, axis);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new CellScalingModifier(args); }
void destroy(Module * m) { delete m; }


