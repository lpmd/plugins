//
//
//

#include "cellscaling.h"

#include <lpmd/md.h>
#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

CellScalingModifier::CellScalingModifier(std::string args): Module("cellscaling")
{
 axis = -1;
 AssignParameter("percent", "0.0");
 AssignParameter("axis", "all");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 percent = GetDouble("percent");
 std::string ax = GetString("axis");
 if ((ax == "all") || (ax == "ALL")) axis = -1;
 if ((ax == "x") || (ax == "X")) axis = 0;
 if ((ax == "y") || (ax == "Y")) axis = 1;
 if ((ax == "z") || (ax == "Z")) axis = 2;
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

CellScalingModifier::~CellScalingModifier() { }

void CellScalingModifier::Show(std::ostream & os) const
{
 if (axis == -1) os << "   " << std::setw(10) << "axis" << " = " << "all\n";
 if (axis == 0) os << "   " << std::setw(10) << "axis" << " = " << "X\n";
 if (axis == 1) os << "   " << std::setw(10) << "axis" << " = " << "Y\n";
 if (axis == 2) os << "   " << std::setw(10) << "axis" << " = " << "Z\n";
 Module::Show(os);
}

void CellScalingModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = cellscaling                                              \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar el tamano de la celda, puede         \n";
 std::cout << " utilizarse para escalar un eje de la celda, o realizar un esacalado           \n";
 std::cout << " hidroestatico. los valores son entregados en forma porcentual.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      percent       : Indica el porcentaje en el que se escalara la celda,     \n";
 std::cout << "                      puede ser positivo(expandir) o negativo(comprimir).      \n";
 std::cout << "      axis          : Indica en que eje se va a comprimir la celda, en caso de \n";
 std::cout << "                      no indicarse se realiza compresion hidrostatica.         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cellscaling                                                               \n";
 std::cout << "     precent 20                                                                \n";
 std::cout << "     axis X                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply cellscaling start=0 each=300 end=1000                                 \n\n";
 std::cout << "      De esta forma escalamos la celda de simulacion en un 20% cada 300 pasos  \n";
 std::cout << " en el eje X, terminando el escalado enel paso 1000.                           \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

}

std::string CellScalingModifier::Keywords() const
{
 return "percent axis start end each";
}

void CellScalingModifier::Apply(SimulationCell & sc)
{
 if (axis == -1)
 {
  std::cerr << "-> Rescaling cell by " << percent << "%\n";  
  sc.RescalePercent(percent); 
 }
 else
 {
  std::cerr<< "-> Rescaling axis " << axis << " in a " << percent << "%\n";
  sc.RescalePercent(percent, axis);
 }
}

void CellScalingModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new CellScalingModifier(args); }
void destroy(Module * m) { delete m; }


