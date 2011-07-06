//
//
//

#include "settag.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

SetTagModifier::SetTagModifier(std::string args): Plugin("settag", "2.1")
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = settag                                                   \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to assign a specific tag to a set of atoms.          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      The specified tag can be any word. However a couple of special words are \n";
 std::cout << " available in the code.                                                        \n";
 std::cout << "      fixedpos      : If an atom have this tag, then the atom have a fixed     \n";
 std::cout << "                      position during the simulation.                          \n";
 std::cout << "      fixedvel      : If an atom have this tag, then the velocity of these     \n";
 std::cout << "                      will be the same during all the simulation.              \n";
 std::cout << "      A real option for tag is value, with this you can set ot change the value\n";
 std::cout << " of the tag.                                                                   \n";
 std::cout << "      value         : True or False. In mostly of the cases is True.           \n"; 
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use settag as floor                                                           \n";
 std::cout << "     tag floor_atoms                                                           \n";
 std::cout << "     value true                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";
 std::cout << " apply floor over box x=0-37 y=0-37 z=27.0-37.0                                \n";
 std::cout << "      With this, we set a tag for a specific set of atoms. With this tag set we\n";
 std::cout << " can now apply different properties to this specific set of atoms.             \n"; 
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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

