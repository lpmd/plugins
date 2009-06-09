//
//
//

#include "cellscaling.h"

#include <lpmd/simulation.h>

#include <iostream>
#include <iomanip>

using namespace lpmd;

CellScalingModifier::CellScalingModifier(std::string args): Plugin("cellscaling", "2.0")
{
 ParamList & params = (*this);
 //
 axis = -1;
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("percent", "0.0");
 DefineKeyword("axis", "all");
 AssignParameter("constant", "true");
 AssignParameter("first", "true");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 percent = double(params["percent"]);
 constant = bool(params["constant"]);
 first = bool(params["first"]);
 std::string ax = params["axis"];
 if ((ax == "all") || (ax == "ALL")) axis = -1;
 if ((ax == "x") || (ax == "X")) axis = 0;
 if ((ax == "y") || (ax == "Y")) axis = 1;
 if ((ax == "z") || (ax == "Z")) axis = 2;
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
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
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar el tamano de la celda, puede         \n";
 std::cout << " utilizarse para escalar un eje de la celda, o realizar un esacalado           \n";
 std::cout << " hidroestatico. los valores son entregados en forma porcentual.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      percent       : Indica el porcentaje en el que se escalara la celda,     \n";
 std::cout << "                      puede ser positivo(expandir) o negativo(comprimir).      \n";
 std::cout << "      axis          : Indica en que eje se va a comprimir la celda, en caso de \n";
 std::cout << "                      no indicarse se realiza compresion hidrostatica.         \n";
 std::cout << "      constant      : true/false, la variacion porcentual es modificada segun  \n";
 std::cout << "                      acomulativamente (constant=false) o no (constant=true).  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cellscaling                                                               \n";
 std::cout << "     precent 20                                                                \n";
 std::cout << "     axis X                                                                    \n";
 std::cout << "     constant true                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply cellscaling start=0 each=300 end=1000                                 \n\n";
 std::cout << "      De esta forma escalamos la celda de simulacion en un 20% cada 300 pasos  \n";
 std::cout << " en el eje X, terminando el escalado enel paso 1000.                           \n";
}

void CellScalingModifier::Apply(Simulation & sim)
{
 lpmd::BasicCell & cell = sim.Cell();
 if (constant == false)
 {
  if (axis == -1)
  {
   DebugStream() << "-> Rescaling cell by " << percent << "%\n";
   //sc.RescalePercent(percent);
   for (int i=0;i<3;++i) cell[i] = cell[i] + cell[i] * percent/100.0e0;
  }
  else
  {
   DebugStream() << "-> Rescaling axis " << axis << " in a " << percent << "%\n";
   //sc.RescalePercent(percent, axis);
   cell[axis] = cell[axis] + cell[axis]*percent/100.0e0;
  }
 }
 else if (constant == true)
 {
  if (first == true)
  {
   for(int i=0;i<3;++i) s[i] = cell[i]*percent/100;
   first=false;
  }
  if (axis == -1)
  {
   DebugStream() << "-> Rescaling constant by " << percent << "%\n";
   //sc.RescaleVector(s[0],s[1],s[2]);
   for (int i=0;i<3;++i) cell[i] = cell[i] + s[i];
  }
  else
  {
   DebugStream() << "-> Rescaling axis " << axis <<" constant by " << percent << "%\n";
   //sc.RescaleVector(s[axis],axis);
   cell[axis] = cell[axis] + s[axis];
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new CellScalingModifier(args); }
void destroy(Plugin * m) { delete m; }


