//
//
//

#include "berendsen.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>
#include <lpmd/integrator.h>
#include <lpmd/properties.h>

#include <iostream>

using namespace lpmd;

BerendsenModifier::BerendsenModifier(std::string args): Plugin("berendsen", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("from", "300.0");
 DefineKeyword("to", "300.0");
 DefineKeyword("tau", "400.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 fromtemp = double(params["from"]);
 totemp = double(params["to"]);
 tau = double(params["tau"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 // stop_thermostat no es un parametro sino un flag interno
 stop_thermostat = -1;
}

BerendsenModifier::~BerendsenModifier() { }

void BerendsenModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar la temperatura del sistema utilizando\n";
 std::cout << " el termostato de berendsen.                                                   \n";
 std::cout << "      Este metodo es mas utilizado que el de rescalamiento de velocidades ya   \n";
 std::cout << " que presenta un decenso menos brusco de la temperatura del sistema.           \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      from          : Temperatura inicial para el escalamiento.                \n";
 std::cout << "      to            : Temperatura final para el sistema.                       \n";
 std::cout << "      tau           : Intervalo del termostato.                                \n";
 std::cout << '\n';
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
}

void BerendsenModifier::Apply(lpmd::Simulation & sim)
{
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::Integrator & integrator = sim.Integrator();
 double timestep = integrator.Timestep();
 double set_temp = ValueAtStep(sim.CurrentStep(), fromtemp, totemp);
 if (stop_thermostat == -1)
 {
  stop_thermostat = long(double(sim.CurrentStep()) + tau/timestep);
  DebugStream() << "-> Berendsen thermostat activated (stops in step " << stop_thermostat << "), rescaling temperature to T = " << set_temp << '\n';
  old_step = each;
  each = 1;
 }
 else
 {
  if (sim.CurrentStep() <= stop_thermostat) 
  {
   Vector vel(0,0,0);
   double xi, ti=lpmd::Temperature(atoms);
   for (int i=0;i<atoms.Size();++i)
   {
    vel = atoms[i].Velocity();
    xi = sqrt(1.0 + (double(timestep)/tau)*(set_temp/ti - 1.0));
    vel = vel*xi;
    atoms[i].Velocity() = vel;
   }
  }
  else
  {
   stop_thermostat = -1;
   DebugStream() << "-> Berendsen thermostat stopped." << '\n';
   each = old_step;
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new BerendsenModifier(args); }
void destroy(Plugin * m) { delete m; }


