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
 DefineKeyword("tag");
 DefineKeyword("value", "true");
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
 if (tag == "") throw PluginError("settag", "The tag to be set was not specified.");
 DebugStream() << "-> Applying tag " << tag << " with value " << value << " to " << atoms.Size() << " atoms\n";  
 BasicParticleSet & orig_atoms = sim.OriginalAtoms();
 for (long int i=0;i<atoms.Size();++i)
 {
  const BasicAtom & at = atoms[i];
  orig_atoms.SetTag(at, Tag(tag), value);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SetTagModifier(args); }
void destroy(Plugin * m) { delete m; }

