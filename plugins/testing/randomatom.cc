//
//
//

#include "randomatom.h"

#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

RandomAtomModifier::RandomAtomModifier(std::string args): Plugin("randomatom", "1.0")
{
 ParamList & params = (*this);
 AssignParameter("type", "del");
 AssignParameter("value", "10");
 AssignParameter("symbol", "e");
 AssignParameter("density", "fixed");
 // 
 ProcessArguments(args);
 type = params["type"];
 value = double(params["val"]);
 symbol = params["sym"];
 density = params["den"];
}

RandomAtomModifier::~RandomAtomModifier() { }

void RandomAtomModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para elmininar/modificar átomos dentro de una     \n";
 std::cout << " muestra. Puede aplicarse al inicio en la instruccion \"prepare\"              \n";
 std::cout << " como durante la simulacion en la instruccion \"apply\".                       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      type          : Tipo de Accion del/rep                                   \n";
 std::cout << "      value         : Valor porcentual de atomos a eliminar/reemplazar         \n";
 std::cout << "      symbol        : simbolo atomico, en el caso de reemplazar.               \n";
 std::cout << "      density       : fixed/free Para fijar o no la densidad de la muestra.    \n";
 std::cout << "                      Toma importancia solo para el modo elminar.              \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare randomatom type=rep value=20 symbol=Cu                                \n";
}

void RandomAtomModifier::Apply(Simulation & conf)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 lpmd::Array<int> random;
 random.Clear();
 if(value<=0) throw PluginError("randomatom","Error in value assignation in plugin.");
 int tochange = (int)(atoms.Size()*value/100.0);
 if(type=="rep")
 {
  while(random.Size()<tochange)
  {
   int rnd = drand48()*atoms.Size();
   if (rnd!=0) random.AppendUnique(rnd);
  }
  for (int i=0;i<random.Size();++i)
  {
   atoms[random[i]].Symbol() = symbol;
  }
 }
 else if(type=="del")
 {
  double density = atoms.Size()/cell.Volume();
 }
 else
 {
  throw PluginError("randomatom","Bad definition in type assignation.");
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new RandomAtomModifier(args); }
void destroy(Plugin * m) { delete m; }

