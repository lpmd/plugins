//
//
//

#include "berendsen.h"

#include <lpmd/simulationcell.h>
#include <lpmd/integrator.h>
#include <lpmd/md.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

BerendsenModifier::BerendsenModifier(std::string args): Module("berendsen")
{
 AssignParameter("from", "300.0");
 AssignParameter("to", "300.0");
 AssignParameter("tau", "400.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 fromtemp = GetDouble("from");
 totemp = GetDouble("to");
 tau = GetDouble("tau");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 // stop_thermostat no es un parametro sino un flag interno
 stop_thermostat = -1;
}

BerendsenModifier::~BerendsenModifier() { }

void BerendsenModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = berendsen                                                \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar la temperatura del sistema utilizando\n";
 std::cout << " el termostato de berendsen.                                                   \n";
 std::cout << "      Este metodo es mas utilizado que el de rescalamiento de velocidades ya   \n";
 std::cout << " que presenta un decenso menos brusco de la temperatura del sistema.           \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      from          : Temperatura inicial para el escalamiento.                \n";
 std::cout << "      to            : Temperatura final para el sistema.                       \n";
 std::cout << "      tau           : Intervalo del termostato.                                \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use berendsen                                                                 \n";
 std::cout << "     tau  400.0                                                                \n";
 std::cout << "     from 84.0                                                                 \n";
 std::cout << "     to   10.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply berendsen start=0 each=10 end=100                                     \n\n";
 std::cout << "      De esta forma aplicamos el termostato entre 0 y 100 cada 10 steps.       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string BerendsenModifier::Keywords() const 
{
 return "from to tau start end each";
}

void BerendsenModifier::Apply(SimulationCell & sc)
{
 ShowWarning("plugin berendsen", "Applying the berendsen modifier to a single configuration does not make much sense.");
}

void BerendsenModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 double timestep = md.GetIntegrator().Timestep();
 double set_temp = LeverRule(md.CurrentStep(), start_step, end_step, fromtemp, totemp);
 if (stop_thermostat == -1)
 {
  stop_thermostat = long(double(md.CurrentStep()) + tau/timestep);
  std::cerr << "-> Berendsen thermostat activated (stops in step " << stop_thermostat << "), rescaling temperature to T = " << set_temp << '\n';
  old_step = interval;
  interval = 1;
 }
 else
 {
  if (md.CurrentStep() <= stop_thermostat) 
  {
   sc.SetTemperature(set_temp, timestep, tau);
  }
  else
  {
   stop_thermostat = -1;
   std::cerr << "-> Berendsen thermostat stopped." << '\n';
   interval = old_step;
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new BerendsenModifier(args); }
void destroy(Module * m) { delete m; }


