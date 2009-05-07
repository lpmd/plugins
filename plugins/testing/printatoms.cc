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
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("from", "-1");
 DefineKeyword("to", "-1");
 ProcessArguments(args);
 from_at = GetInteger("from");
 to_at = GetInteger("to");
 start = GetInteger("start");
 end = GetInteger("end");
 each = GetInteger("each");
}

PrintAtomsVisualizer::~PrintAtomsVisualizer() { }

void PrintAtomsVisualizer::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use printatoms                                                                \n";
 std::cout << "     from 10                                                                   \n";
 std::cout << "     to 20                                                                     \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " visualize printatoms start=1 end=1000 each=50                               \n\n";
}

void PrintAtomsVisualizer::Apply(const MD & md) 
{ 
 SimulationCell & sc = md.GetCell();
 for (unsigned long int i=0;i<sc.size();++i) 
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

