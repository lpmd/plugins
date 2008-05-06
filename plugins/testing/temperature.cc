//
//
//

#include "temperature.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

TemperatureModifier::TemperatureModifier(std::string args): Module("temperature")
{
 AssignParameter("t", "300.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 temp = GetDouble("t");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

TemperatureModifier::~TemperatureModifier() { }

void TemperatureModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = temperature                                              \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para fijar la temperatura del sistema utilizando  \n";
 std::cout << " una distribucion uniforme de velocidades compatible con la temperatura.       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      temp          : Temperatura deseada para el sistema.                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare temperature t=600.0                                                   \n";
 std::cout << "\n";
 std::cout << "      De esta forma fijamos la temperatura inicial.                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string TemperatureModifier::Keywords() const
{
 return "t start end each";
}

void TemperatureModifier::Apply(SimulationCell & sc)
{
 sc.InitVelocities();
 sc.SetTemperature(temp);
}

void TemperatureModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 std::cerr << "-> Rescaling temperature to T = " << temp << '\n';  
 sc.InitVelocities();
 sc.SetTemperature(temp);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new TemperatureModifier(args); }
void destroy(Module * m) { delete m; }


