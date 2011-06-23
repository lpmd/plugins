//
//
//

#include "displace.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

DisplaceModifier::DisplaceModifier(std::string args): Plugin("displace", "2.0")
{
 ParamList & params = (*this);
 DefineKeyword("x", "1.0");
 DefineKeyword("y", "0.0");
 DefineKeyword("z", "0.0");
 // 
 ProcessArguments(args);
 offset = Vector(double(params["x"]), double(params["y"]), double(params["z"]));
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

DisplaceModifier::~DisplaceModifier() { }

void DisplaceModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      x             : Coord. X del vector de desplazamiento                    \n";
 std::cout << "      y             : Coord. Y del vector de desplazamiento                    \n";
 std::cout << "      z             : Coord. Z del vector de desplazamiento                    \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare displace x=1.0 y=0.0 z=0.0                                            \n";
}

void DisplaceModifier::Apply(lpmd::Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicCell & cell = sim.Cell();
 for (long int i=0;i<atoms.Size();++i)
 {
  atoms[i].Position() = cell.FittedInside(atoms[i].Position()+offset);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new DisplaceModifier(args); }
void destroy(Plugin * m) { delete m; }

