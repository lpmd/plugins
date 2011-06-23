//
//
//

#include "pinatom.h"

#include <lpmd/simulation.h>

#include <iostream>
#include <iomanip>

using namespace lpmd;

PinAtomModifier::PinAtomModifier(std::string args): Plugin("pinatom", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("index", "0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 index = int(params["index"]);
 initial_pos = NULL;
}

PinAtomModifier::~PinAtomModifier() { delete initial_pos; }

void PinAtomModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para fijar la posicion de un atomo en particular y\n";
 std::cout << "      conservando los desplazamientos respecto a ese atomo.                    \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example:                                                                      \n";
 std::cout << " use pinatom                                                                   \n";
 std::cout << "     index 42                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply pinatom start=0 each=300 end=1000                                     \n\n";
}

void PinAtomModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 if (initial_pos == NULL) 
 {
  initial_pos = new Vector(atoms[index].Position());
 }
 Vector correction = (*initial_pos) - atoms[index].Position();
 for (long i=0;i<atoms.Size();++i)
 {
  atoms[i].Position() = atoms[i].Position() + correction;  
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new PinAtomModifier(args); }
void destroy(Plugin * m) { delete m; }

