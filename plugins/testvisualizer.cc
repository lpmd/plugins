//
//
//

#include "testvisualizer.h"

using namespace lpmd;

TestVisualizer::TestVisualizer(std::string args): Module("testvisualizer")
{
 ProcessArguments(args);
}

TestVisualizer::~TestVisualizer() { }

void TestVisualizer::SetParameter(std::string name)
{
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

void TestVisualizer::Show() const
{
 Module::Show();
 std::cout << "   start = " << start_step << '\n';
 std::cout << "   end   = " << end_step << '\n';
 std::cout << "   each  = " << interval << '\n';
}

std::string TestVisualizer::Keywords() const
{
 return "start end each";
}

void TestVisualizer::Apply(const MD & md)
{
 SimulationCell & sc = md.GetCell();
 std::cerr << "DEBUG Applying TestVisualizer!" << '\n';
 // 
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new TestVisualizer(args); }
void destroy(Module * m) { delete m; }


