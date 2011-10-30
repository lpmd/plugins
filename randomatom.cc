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
 DefineKeyword("type", "delete");
 DefineKeyword("percent", "10");
 DefineKeyword("symbol", "e");
 DefineKeyword("density", "fixed");
 // 
 ProcessArguments(args);
 type = params["type"];
 percent = double(params["percent"]);
 symbol = params["symbol"];
 density = params["density"];
}

RandomAtomModifier::~RandomAtomModifier() { }

void RandomAtomModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = randomatom                                               \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This module is used to delete or replace a percentage of atoms of the    \n";
 std::cout << "      simulation cell (see also 'random' plugin).                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      type          : Sets the type of action to be done (delete / replace).   \n";
 std::cout << "      percent       : Sets the percentage of atoms that will be deleted or     \n";
 std::cout << "                      replaced.                                                \n";
 std::cout << "      symbol        : In the case 'replace', sets the atomic symbol of the     \n";
 std::cout << "                      atoms that will be replaced.                             \n";
 std::cout << "      density       : In the case 'delete', sets if the density is forced to   \n";
 std::cout << "                      remain constant or not (fixed / free).                   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use randomatom                                                                \n";
 std::cout << "     type replace                                                              \n";
 std::cout << "     percent 10                                                                \n";
 std::cout << "     symbol Cu                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";
 std::cout << " prepare randomatom type=replace percent=20 symbol=Cu                          \n";
 std::cout << "      The plugin is used to replace 10\% of the atoms of the simulation cell by\n";
 std::cout << "      copper (Cu) atoms.                                                       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void RandomAtomModifier::Apply(Simulation & conf)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 lpmd::Array<int> random;
 random.Clear();
 if(percent<=0) throw PluginError("randomatom","Error in percent assignation in plugin.");
 int tochange = (int)(atoms.Size()*percent/100.0);
 if(type=="replace")
 {
  while(random.Size()<tochange)
  {
   int rnd = int(drand48()*atoms.Size());
   if (rnd!=0) random.AppendUnique(rnd);
  }
  DebugStream() << "-> Replacing " << random.Size() << " atoms\n";
  for (int i=0;i<random.Size();++i)
  {
   atoms[random[i]].Z() = ElemNum(symbol);
   assert (atoms[random[i]].Symbol() == symbol);
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
  DebugStream() << "-> Deleting " << random.Size() << " atoms\n";
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

