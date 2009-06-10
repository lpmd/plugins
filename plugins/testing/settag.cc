//
//
//

#include "settag.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

SetTagModifier::SetTagModifier(std::string args): Plugin("settag", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("velocity");
 DefineKeyword("filterby", "none");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 tag = params["tag"];
 value = params["value"];
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

SetTagModifier::~SetTagModifier() { }

void SetTagModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << '\n';
}

void SetTagModifier::Apply(Simulation & sim)
{
 BasicParticleSet & atoms = sim.Atoms();
 DebugStream() << "-> Applying tag " << tag << " to " << atoms.Size() << " atoms\n";  
 for (long int i=0;i<atoms.Size();++i) atoms.SetTag(atoms[i],Tag(tag),value);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SetTagModifier(args); }
void destroy(Plugin * m) { delete m; }

