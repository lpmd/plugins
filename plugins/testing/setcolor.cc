//
//
//

#include "setcolor.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

SetColorModifier::SetColorModifier(std::string args): Plugin("setcolor", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("color");
 DefineKeyword("filterby", "none");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 color = Color(params["color"].c_str());
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

SetColorModifier::~SetColorModifier() { }

void SetColorModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << '\n';
}

void SetColorModifier::Apply(Simulation & sim)
{
 BasicParticleSet & atoms = sim.Atoms();
 DebugStream() << "-> Applying color ";
 FormattedWrite(DebugStream(), color);
 DebugStream() << " to " << atoms.Size() << " atoms\n";  
 for (long int i=0;i<atoms.Size();++i) ColorHandler::ColorOfAtom(atoms[i]) = color;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SetColorModifier(args); }
void destroy(Plugin * m) { delete m; }

