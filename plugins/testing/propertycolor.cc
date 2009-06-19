//
//
//

#include "propertycolor.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>
#include <lpmd/error.h>
#include <lpmd/session.h>

#include <iostream>

using namespace lpmd;

PropertyColorModifier::PropertyColorModifier(std::string args): Plugin("propertycolor", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("property");
 DefineKeyword("min");
 DefineKeyword("max");
 DefineKeyword("filterby", "none");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 property = params["property"];
 vmin = double(params["min"]);
 vmax = double(params["max"]);
}

PropertyColorModifier::~PropertyColorModifier() { }

void PropertyColorModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << '\n';
}

void PropertyColorModifier::Apply(Simulation & sim)
{
 const double kin2ev = double(GlobalSession["kin2ev"]);
 const double kboltzmann = double(GlobalSession["kboltzmann"]);
 BasicParticleSet & atoms = sim.Atoms();
 DebugStream() << "-> Applying color according to " << property << '\n';

 for (long int i=0;i<atoms.Size();++i)
 {
  double v = 0.0;
  // 
  // Por ahora seleccionamos aqui de donde sale el escalar para asignar el color
  // Esto hasta tener un sistema mas dinamico donde no esten en el codigo todas las posibilidades
  //
  if (property == "temperature") v = (1.0/3.0)*(kin2ev/kboltzmann)*atoms[i].Mass()*atoms[i].Velocity().SquareModule();
  else if (property == "random") v = vmin + (vmax-vmin)*drand48();
  else throw PluginError("propertycolor", "Cannot color atoms by property \""+property+"\"");
  double vnorm = (v-vmin)/(vmax-vmin);
  if (vnorm < 0.0) vnorm = 0.0;
  if (vnorm > 1.0) vnorm = 1.0;
  const Color c = lpmd::ColorFromScalar(1.0-vnorm);
  const BasicAtom & at = atoms[i]; 
  ColorHandler::ColorOfAtom(at) = c;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new PropertyColorModifier(args); }
void destroy(Plugin * m) { delete m; }

