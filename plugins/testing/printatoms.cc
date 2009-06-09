//
//
//

#include "printatoms.h"

#include <lpmd/simulation.h>

#include <iostream>

using namespace lpmd;

PrintAtomsVisualizer::PrintAtomsVisualizer(std::string args): Plugin("printatoms", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("from", "-1");
 DefineKeyword("to", "-1");
 ProcessArguments(args);
 from_at = int(params["from"]);
 to_at = int(params["to"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
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

void PrintAtomsVisualizer::Apply(const Simulation & sim) 
{ 
 const BasicParticleSet & atoms = sim.Atoms();
 for (long int i=0;i<atoms.Size();++i) 
 {  
  if ((i >= from_at) && (i <= to_at)) 
  {
   std::cout << "-> Atom " << i << " : " << atoms[i].Symbol() << '\n'; 
   std::cout << "   Position     : " << atoms[i].Position() << '\n';
   std::cout << "   Velocity     : " << atoms[i].Velocity() << '\n';
   std::cout << "   Acceleration : " << atoms[i].Acceleration() << '\n';
  }
 }
 std::cout << '\n';
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new PrintAtomsVisualizer(args); }
void destroy(Plugin * m) { delete m; }

