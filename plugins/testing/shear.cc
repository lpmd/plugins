//
//
//

#include "shear.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

ShearModifier::ShearModifier(std::string args): Module("shear")
{
 AssignParameter("axis", "X");
 AssignParameter("normal", "Y");
 AssignParameter("strain", "0.01");
 // 
 ProcessArguments(args);
 if ((GetString("axis") == "x") || (GetString("axis") == "X")) shear_axis = 0;
 else if ((GetString("axis") == "y") || (GetString("axis") == "Y")) shear_axis = 1;
 else if ((GetString("axis") == "z") || (GetString("axis") == "Z")) shear_axis = 2;
 else throw PluginError("shear", "Invalid shear axis");
 if ((GetString("normal") == "x") || (GetString("normal") == "X")) perp_axis = 0;
 else if ((GetString("normal") == "y") || (GetString("normal") == "Y")) perp_axis = 1;
 else if ((GetString("normal") == "z") || (GetString("normal") == "Z")) perp_axis = 2;
 else throw PluginError("shear", "Invalid normal axis");
 if (shear_axis == perp_axis) throw PluginError("shear", "Shear axis must be orthogonal to normal axis");
 strain = GetDouble("strain");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

ShearModifier::~ShearModifier() { }

void ShearModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = shear                                                    \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      axis          : Eje en que se produce el cizalle                         \n";
 std::cout << "      normal        : Eje perpendicular al eje de cizalle                      \n";
 std::cout << "      strain        : Desplazamiento maximo a aplicar es strain*L(normal)      \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare shear axis=X normal=Y strain=0.01                                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string ShearModifier::Keywords() const
{
 return "axis normal strain start end each";
}

void ShearModifier::Apply(SimulationCell & sc)
{
 // primero cambia la forma de la celda...
 Vector deformation(0.0, 0.0, 0.0);
 deformation[shear_axis] = strain*sc.GetVector(perp_axis).Module(); 
 Vector newaxis = sc.GetVector(perp_axis)+deformation;
 newaxis = newaxis*double(sc.GetVector(perp_axis).Module())/double(newaxis.Module());
 sc.SetVector(perp_axis, newaxis);
 
 // luego desplaza los atomos
 Vector offset(0.0, 0.0, 0.0);
 for (unsigned long int i=0;i<sc.size();++i)
 {
  Vector pos = sc[i].Position(); 
  offset[shear_axis] = pos[perp_axis]*strain; 
  sc.SetPosition(i, pos+offset);
 }
}

void ShearModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 Apply(sc);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ShearModifier(args); }
void destroy(Module * m) { delete m; }

