//
//
//

#include "shear.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>
#include <lpmd/plugin.h>

#include <iostream>

using namespace lpmd;

ShearModifier::ShearModifier(std::string args): Plugin("shear", "2.0")
{
 ParamList & param = (*this);
 AssignParameter("axis", "X");
 AssignParameter("normal", "Y");
 AssignParameter("strain", "0.01");
 // 
 ProcessArguments(args);
 if ((param["axis"] == "x") || (param["axis"] == "X")) shear_axis = 0;
 else if ((param["axis"] == "y") || (param["axis"] == "Y")) shear_axis = 1;
 else if ((param["axis"] == "z") || (param["axis"] == "Z")) shear_axis = 2;
 else throw PluginError("shear", "Invalid shear axis");
 if ((param["normal"] == "x") || (param["normal"] == "X")) perp_axis = 0;
 else if ((param["normal"] == "y") || (param["normal"] == "Y")) perp_axis = 1;
 else if ((param["normal"] == "z") || (param["normal"] == "Z")) perp_axis = 2;
 else throw PluginError("shear", "Invalid normal axis");
 if (shear_axis == perp_axis) throw PluginError("shear", "Shear axis must be orthogonal to normal axis");
 strain = double(param["strain"]);
 start = int(param["start"]);
 end = int(param["end"]);
 each = int(param["each"]);
}

ShearModifier::~ShearModifier() { }

void ShearModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      axis          : Eje en que se produce el cizalle                         \n";
 std::cout << "      normal        : Eje perpendicular al eje de cizalle                      \n";
 std::cout << "      strain        : Desplazamiento maximo a aplicar es strain*L(normal)      \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare shear axis=X normal=Y strain=0.01                                     \n";
}

void ShearModifier::Apply(Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicCell & cell = sim.Cell();
 // primero cambia la forma de la celda...
 Vector deformation(0.0, 0.0, 0.0);
 deformation[shear_axis] = strain*cell[perp_axis].Module(); 
 Vector newaxis = cell[perp_axis]+deformation;
 newaxis = newaxis*double(cell[perp_axis].Module())/double(newaxis.Module());
 cell[perp_axis] = newaxis;
 
 // luego desplaza los atomos
 Vector offset(0.0, 0.0, 0.0);
 for (long int i=0;i<atoms.Size();++i)
 {
  Vector pos = atoms[i].Position(); 
  offset[shear_axis] = pos[perp_axis]*strain; 
  atoms[i].Position() = pos+offset;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new ShearModifier(args); }
void destroy(Plugin * m) { delete m; }

