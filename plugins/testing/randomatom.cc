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
 AssignParameter("type", "delete");
 AssignParameter("value", "10");
 AssignParameter("symbol", "e");
 AssignParameter("density", "fixed");
 // 
 ProcessArguments(args);
 type = params["type"];
 value = double(params["value"]);
 symbol = params["symbol"];
 density = params["density"];
}

RandomAtomModifier::~RandomAtomModifier() { }

void RandomAtomModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para elmininar/modificar átomos dentro de una     \n";
 std::cout << " muestra. Puede aplicarse al inicio en la instruccion \"prepare\"              \n";
 std::cout << " como durante la simulacion en la instruccion \"apply\".                       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      type          : Tipo de Accion delete/replace                            \n";
 std::cout << "      value         : Valor porcentual de atomos a eliminar/reemplazar         \n";
 std::cout << "      symbol        : simbolo atomico, en el caso de reemplazar.               \n";
 std::cout << "      density       : fixed/free Para fijar o no la densidad de la muestra.    \n";
 std::cout << "                      Toma importancia solo para el modo elminar.              \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare randomatom type=replace value=20 symbol=Cu                            \n";
}

void RandomAtomModifier::Apply(Simulation & conf)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 lpmd::Array<int> random;
 random.Clear();
 if(value<=0) throw PluginError("randomatom","Error in value assignation in plugin.");
 int tochange = (int)(atoms.Size()*value/100.0);
 if(type=="replace")
 {
  while(random.Size()<tochange)
  {
   int rnd = int(drand48()*atoms.Size());
   if (rnd!=0) random.AppendUnique(rnd);
  }
  for (int i=0;i<random.Size();++i)
  {
   atoms[random[i]].Symbol() = symbol;
  }
 }
 else if(type=="delete")
 {
  double mass = 0.0e0;
  for(int i=0;i<atoms.Size();++i)
  {
   mass += atoms[i].Mass();
  }
  double density = mass/cell.Volume();
  while(random.Size()<tochange)
  {
   int rnd = int(drand48()*atoms.Size());
   if (rnd!=0) random.AppendUnique(rnd);
  }
  //first delete atoms.
  for (int i=0;i<random.Size();++i)
  {
   atoms.Delete(random[i]);
  }
  //adjust cell parameters to original density.
  mass = 0.0e0;
  for(int i=0;i<atoms.Size();++i)
  {
   mass += atoms[i].Mass();
  }
  double newdens=0.5*density;
  while(fabs(newdens-density)>1E-2)
  {
   if(newdens > density)
   {
    for (int j=0;j<3;++j) cell[j] = cell[j] + cell[j]*0.01;
   }
   else if(newdens < density)
   {
    for (int j=0;j<3;++j) cell[j] = cell[j] - cell[j]*0.01;
   }
   newdens = mass/cell.Volume();
  }
 }
 else
 {
  throw PluginError("randomatom","Bad definition in type assignation.");
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new RandomAtomModifier(args); }
void destroy(Plugin * m) { delete m; }

