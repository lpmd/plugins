//
//
//

#include "skewstart.h"

#include <lpmd/atom.h>
#include <lpmd/configuration.h>
#include <cmath>

using namespace lpmd;

SkewStartGenerator::SkewStartGenerator(std::string args): Plugin("skewstart", "2.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("atoms");
 DefineKeyword("symbol", "H");
 ProcessArguments(args); 
 spc = ElemNum(params["symbol"]);
 n = int(params["atoms"]);
}

SkewStartGenerator::~SkewStartGenerator() { }

void SkewStartGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = skewstart                                                \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin is used to generate a initial simulation cell using the       \n";
 std::cout << " skewstart method development by Keith Refson for the moldy software.          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      atoms         : Set the number of atoms in the simulation cell.          \n";
 std::cout << "      symbol        : Set the atomic symbol for the simulation cell. (only one)\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " #Using the plugin :                                                           \n";
 std::cout << " input module=skewstart atoms=108 symbol=Ar                                    \n";
 std::cout << "      With this we generate a simulation cell using  the skewstart method, of  \n";
 std::cout << " 108 Argon atoms.                                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

//
// Skew Start method, like Moldy
// Only works with atoms (we have no molecules yet!)
//
void SkewStartGenerator::Generate(Configuration & config) const
{
 BasicParticleSet & atoms = config.Atoms();
 BasicCell & cell = config.Cell();  
 int h, k, l;
 double dx, dy, dz;
 h = int(pow(double(n), 2.0/3.0));
 k = int(pow(double(n), 1.0/3.0));
 l = 1;
 dx = h / double(n);
 dy = k / double(n);
 dz = l / double(n);
 bool create_atoms = (atoms.Size() == 0);
 for (long i=0;i<n;++i)
 {
  if (create_atoms) atoms.Append(Atom(spc));
  atoms[i].Position() = cell.FittedInside(cell.Cartesian(Vector(dx*double(i), dy*double(i), dz*double(i))));
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SkewStartGenerator(args); }
void destroy(Plugin * m) { delete m; }
