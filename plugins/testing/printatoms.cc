//
//
//

#include "printatoms.h"

#include <lpmd/simulation.h>
#include <lpmd/colorhandler.h>

#include <iostream>

using namespace lpmd;

PrintAtomsVisualizer::PrintAtomsVisualizer(std::string args): Plugin("printatoms", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("tags", "x,y,z");
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 tags = StringSplit(params["tags"], ',');
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
 std::cout << "-> ";
 for (long int i=0;i<atoms.Size();++i) 
 {  
  const Vector & pos = atoms[i].Position();
  const Vector & vel = atoms[i].Velocity();
  const Vector & acc = atoms[i].Acceleration();
  for (int j=0;j<tags.Size();++j)
  {
   const std::string & tag = tags[j];
   std::cout << tag << "=";
   if (tag == "x") std::cout << pos[0];
   else if (tag == "y") std::cout << pos[1];
   else if (tag == "z") std::cout << pos[2];
   else if (tag == "vx") std::cout << vel[0];
   else if (tag == "vy") std::cout << vel[1];
   else if (tag == "vz") std::cout << vel[2];
   else if (tag == "ax") std::cout << acc[0];
   else if (tag == "ay") std::cout << acc[1];
   else if (tag == "az") std::cout << acc[2];
   else if ((tag == "r") && ColorHandler::HaveColor(atoms[i])) std::cout << ColorHandler::ColorOfAtom(atoms[i])[0];
   else if ((tag == "g") && ColorHandler::HaveColor(atoms[i])) std::cout << ColorHandler::ColorOfAtom(atoms[i])[1];
   else if ((tag == "b") && ColorHandler::HaveColor(atoms[i])) std::cout << ColorHandler::ColorOfAtom(atoms[i])[2];
   else if (atoms.Have(atoms[i], Tag(tag))) std::cout << atoms.GetTag(atoms[i], Tag(tag));
   std::cout << " ";
  }
  std::cout << '\n';
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new PrintAtomsVisualizer(args); }
void destroy(Plugin * m) { delete m; }

