//
//
//

#include "temperature.h"

#include <lpmd/simulation.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

TemperatureModifier::TemperatureModifier(std::string args): Plugin("temperature", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("t", "300.0");
 DefineKeyword("filterby", "none");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 temp = double(params["t"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
}

TemperatureModifier::~TemperatureModifier() { }

void TemperatureModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para fijar la temperatura del sistema utilizando  \n";
 std::cout << " una distribucion uniforme de velocidades compatible con la temperatura.       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      temp          : Temperatura deseada para el sistema.                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare temperature t=600.0                                                   \n";
 std::cout << '\n';
 std::cout << "      De esta forma fijamos la temperatura inicial.                            \n";
}

void TemperatureModifier::Apply(Simulation & sim)
{
 DebugStream() << "-> Rescaling temperature to T = " << temp << '\n';  
 sim.SetTemperature(temp);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new TemperatureModifier(args); }
void destroy(Plugin * m) { delete m; }


