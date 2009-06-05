//
//
//

#include "displace.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

DisplaceModifier::DisplaceModifier(std::string args): Module("displace")
{
 ParamList & params = (*this);
 AssignParameter("x", "1.0");
 AssignParameter("y", "0.0");
 AssignParameter("z", "0.0");
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = displace                                                 \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      x             : Coord. X del vector de desplazamiento                    \n";
 std::cout << "      y             : Coord. Y del vector de desplazamiento                    \n";
 std::cout << "      z             : Coord. Z del vector de desplazamiento                    \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare displace x=1.0 y=0.0 z=0.0                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string DisplaceModifier::Keywords() const
{
 return "x y z start end each";
}

void DisplaceModifier::Apply(lpmd::Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();

 for (long int i=0;i<atoms.Size();++i)
 {
  atoms[i].Position() = atoms[i].Position() + offset;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new DisplaceModifier(args); }
void destroy(Module * m) { delete m; }

