//
//
//

#include "color.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

ColorModifier::ColorModifier(std::string args): Module("tempscaling")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "2.0"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("type", "vel");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 type = GetString("type");
}

ColorModifier::~ColorModifier() { }

void ColorModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para modificar el color de los átomos utilizando  \n";
 std::cout << " distintas características.                                                    \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      type          : vel,ace -> Tipos soportados.                             \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use color                                                                     \n";
 std::cout << "     vel                                                                       \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply color start=0 each=10 end=100                                         \n\n";
 std::cout << "     Modifica el color segun velocidad del atomo en los instantes indicados.   \n";
}

void ColorModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 long N = sc.size();
 double max=0.0;
 double color = 0.0;
 double *info;
 if (N >= 1) info = new double[N];
 else throw PluginError ("color", "Error in atom numbers less that 1.");
 if (type=="vel")
 {
  for (long i = 0 ; i<N ; ++i)
  {
   info[i] = 0.0e0;
   double vel = (sc[i].Velocity()).Mod();
   if (vel >= max) max = vel;
   info[i] = vel;
  }
  for (long i = 0 ; i < N ; ++i)
  {
   color = fabs(info[i]/max);
   sc[i].SetColor(color);
  }
 }
 else if(type=="ace")
 {
  for (long i = 0 ; i<N ; ++i)
  {
   info[i] = 0.0e0;
   double ace = (sc[i].Acceleration()).Mod();
   if (ace >= max) max = ace;
   info[i] = ace;
  }
  for (long i = 0 ; i < N ; ++i)
  {
   color = fabs(info[i]/max);
   sc[i].SetColor(color);
  }
 }
 else throw PluginError("color", "Error in plugin color type undefined");
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ColorModifier(args); }
void destroy(Module * m) { delete m; }


