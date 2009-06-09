//
//
//

#include "tempscaling.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

TempScalingModifier::TempScalingModifier(std::string args): Plugin("tempscaling", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("from", "300.0");
 DefineKeyword("to", "300.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 fromtemp = double(params["from"]);
 totemp = double(params["to"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

TempScalingModifier::~TempScalingModifier() { }

void TempScalingModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar la temperatura del sistema utilizando\n";
 std::cout << " rescalamiento de velocidades.                                                 \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      from          : Temperatura inicial para el escalamiento.                \n";
 std::cout << "      to            : Temperatura final para el sistema.                       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use tempscaling                                                               \n";
 std::cout << "     from 84.0                                                                 \n";
 std::cout << "     to   10.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply tempscaling start=0 each=10 end=100                                     \n\n";
 std::cout << "      De esta forma aplicamos el termostato entre 0 y 100 cada 10 steps.       \n";
}

void TempScalingModifier::Apply(Simulation & sim)
{
 double set_temp = ValueAtStep(sim.CurrentStep(), fromtemp, totemp);
 DebugStream() << "-> Rescaling temperature to T = " << set_temp << '\n';  
 sim.SetTemperature(set_temp);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new TempScalingModifier(args); }
void destroy(Plugin * m) { delete m; }

