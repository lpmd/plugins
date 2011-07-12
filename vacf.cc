//
//
//

#include "vacf.h"
#include <lpmd/properties.h>

using namespace lpmd;

Vacf::Vacf(std::string args): Plugin("vacf", "2.1")
{
 ParamList & param = (*this);
 DefineKeyword("dt","1");
 //
 ProcessArguments(args);
 dt = param["dt"];
}

Vacf::~Vacf() { }

void VariableStep::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = vacf                                                     \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to evaluate velocity autocorrelation function of a   \n";
 std::cout << " set of atomic configurations from a file.                                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Set the time step in femto-seconds in case that the      \n";
 std::cout << "                      velocity are not specified in the input file.            \n";
 std::cout << "      output        : Set the output file-name.                                \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use varstep                                                                   \n";
 std::cout << "     output vacf.dat                                                           \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator varstep start=1000                                                 \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}


void Vacf::Evaluate(lpmd::ConfigurationSet & hist, Potential & pot)
{
 vacf(hist,pot,dt,CurrentValue());
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Vacf(args); }
void destroy(Plugin * m) { delete m; }
