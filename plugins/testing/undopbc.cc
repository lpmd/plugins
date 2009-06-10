//
//
//

#include "undopbc.h"

#include <lpmd/simulation.h>

#include <iostream>
#include <iomanip>

using namespace lpmd;

UndoPBCModifier::UndoPBCModifier(std::string args): Plugin("cellscaling", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
}

UndoPBCModifier::~UndoPBCModifier() { }

void UndoPBCModifier::Show(std::ostream & os) const
{
 Module::Show(os);
}

void UndoPBCModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para posicionar los atomos sin condiciones        \n";
 std::cout << " periodicas de borde.                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      percent       : Indica el porcentaje en el que se escalara la celda,     \n";
 std::cout << "                      puede ser positivo(expandir) o negativo(comprimir).      \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use undopbc                                                                   \n";
 std::cout << "     constant true                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply undopbc start=0 each=300 end=1000                                     \n\n";
 std::cout << "      De esta forma modificamos las posiciones atomicas de nuestra muestra     \n";
 std::cout << " desde el paso 0 hasta el 1000 cada 300.                                       \n";
}

void UndoPBCModifier::Initialize(Simulation & sim, Potential & pot){old  = sim;}

void UndoPBCModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicCell & cell = sim.Cell();
// int N=2;

 Vector ** noperiodic = new Vector*[N];
 for (int t=0;t<N;++t) noperiodic[t] = new Vector[nat];
 for (int i=0;i<nat;++i) noperiodic[0][i] = part[i].Position();

// for (int t=1;t<N;++t)
  for (int i=0;i<nat;++i)
  {
   atoms[0].Position() = old.Atoms()[i].Position();
   atoms[1].Position() = sim.Atoms()[i].Position();
   noperiodic[t][i] = noperiodic[t-1][i] + cell.Displacement(part[0].Position(), part[1].Position());
  }


}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new UndoPBCModifier(args); }
void destroy(Plugin * m) { delete m; }


