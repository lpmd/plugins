//
//
//

#include "printatoms.h"

#include <lpmd/md.h>
#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

PrintAtomsVisualizer::PrintAtomsVisualizer(std::string args): Module("printatoms")
{
 AssignParameter("from", "-1");
 AssignParameter("to", "-1");
 ProcessArguments(args);
 from_at = GetInteger("from");
 to_at = GetInteger("to");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

PrintAtomsVisualizer::~PrintAtomsVisualizer() { }

void PrintAtomsVisualizer::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = printatoms                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use printatoms                                                                \n";
 std::cout << "     from 10                                                                   \n";
 std::cout << "     to 20                                                                     \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " visualize printatoms start=1 end=1000 each=50                               \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string PrintAtomsVisualizer::Keywords() const
{
 return "from to";
}

void PrintAtomsVisualizer::Apply(const MD & md) 
{ 
 SimulationCell & sc = md.GetCell();
 for (long i=0;i<sc.Size();++i) 
 {  
  if ((i >= from_at) && (i <= to_at)) 
  {
   std::cout << "-> Atom " << i << " : " << sc[i].Symb() << '\n'; 
   std::cout << "   Position     : " << sc[i].Position() << '\n';
   std::cout << "   Velocity     : " << sc[i].Velocity() << '\n';
   std::cout << "   Acceleration : " << sc[i].Acceleration() << '\n';
  }
 }
 std::cout << '\n';
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new PrintAtomsVisualizer(args); }
void destroy(Module * m) { delete m; }

