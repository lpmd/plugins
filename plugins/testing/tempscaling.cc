//
//
//

#include "tempscaling.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

TempScalingModifier::TempScalingModifier(std::string args): Module("tempscaling")
{
 AssignParameter("from", "300.0");
 AssignParameter("to", "300.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 fromtemp = GetDouble("from");
 totemp = GetDouble("to");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

TempScalingModifier::~TempScalingModifier() { }

void TempScalingModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = tempscaling                                              \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar la temperatura del sistema utilizando\n";
 std::cout << " rescalamiento de velocidades.                                                 \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      from          : Temperatura inicial para el escalamiento.                \n";
 std::cout << "      to            : Temperatura final para el sistema.                       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use tempscaling                                                               \n";
 std::cout << "     from 84.0                                                                 \n";
 std::cout << "     to   10.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply tempscaling start=0 each=10 end=100                                     \n\n";
 std::cout << "      De esta forma aplicamos el termostato entre 0 y 100 cada 10 steps.       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string TempScalingModifier::Keywords() const
{
 return "from to start end each";
}

void TempScalingModifier::Apply(SimulationCell & sc)
{
 sc.SetTemperature(fromtemp);
}

void TempScalingModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 double set_temp = LeverRule(md.CurrentStep(), start_step, end_step, fromtemp, totemp);
 std::cerr << "-> Rescaling temperature to T = " << set_temp << '\n';  
 sc.SetTemperature(set_temp);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new TempScalingModifier(args); }
void destroy(Module * m) { delete m; }


