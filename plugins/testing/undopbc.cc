//
//
//

#include "undopbc.h"

#include <lpmd/simulation.h>

#include <iostream>
#include <iomanip>

using namespace lpmd;

UndoPBCModifier::UndoPBCModifier(std::string args): Plugin("cellscaling", "2.0"), oldpositions(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 first_apply = true;
}

UndoPBCModifier::~UndoPBCModifier() { delete [] oldpositions; }

void UndoPBCModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para posicionar los atomos sin condiciones        \n";
 std::cout << " periodicas de borde.                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example:                                                                      \n";
 std::cout << " use undopbc                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply undopbc start=0 each=300 end=1000                                     \n\n";
 std::cout << "      De esta forma modificamos las posiciones atomicas de nuestra muestra     \n";
 std::cout << " desde el paso 0 hasta el 1000 cada 300.                                       \n";
}

void UndoPBCModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicCell & cell = sim.Cell();
 if (first_apply)
 {
  oldpositions = new Vector[atoms.Size()];
  for (long int i=0;i<atoms.Size();++i) oldpositions[i] = atoms[i].Position(); 
  first_apply = false;
 }
 else
 {
  for (long int i=0;i<atoms.Size();++i)
  {
   atoms[i].Position() = oldpositions[i] + cell.Displacement(oldpositions[i], atoms[i].Position());
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new UndoPBCModifier(args); }
void destroy(Plugin * m) { delete m; }

