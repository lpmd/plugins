//
//
//

#include "gdr.h"

#include <lpmd/simulation.h>
#include <lpmd/properties.h>

#include <sstream>

using namespace lpmd;

Gdr::Gdr(std::string args): Plugin("gdr", "2.0")
{
 ParamList & params=(*this);
 //
 DefineKeyword("rcut");
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("output");
 DefineKeyword("average", "false");
 DefineKeyword("bins", "200");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = params["rcut"];
 nb = int(params["bins"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 do_average = bool(params["average"]);
}

Gdr::~Gdr() { }

void Gdr::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular la funcion de autocorrelacion de pares    \n";
 std::cout << " utiliza las condiciones de borde periodicas de la celda.                      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Especifica el numero de divisiones entre 0 y rcut.       \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo de gdr.       \n";
 std::cout << "      output        : Fichero en el que se graba el RDF.                       \n";
 std::cout << "      average       : Setea si calculo o no el promedio de cada calculo.       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use gdr                                                                       \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     output filegdr.dat                                                        \n";
 std::cout << "     rcut 15.0                                                                 \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property gdr start=1 each=10 end=100                                        \n\n";
 std::cout << "      De esta forma calculamos la funcion de distribucion radial de pares en    \n";
 std::cout << " la simulacion entre 1 y 100 cada 10 pasos.                                    \n";
}

void Gdr::Evaluate(Configuration & con, Potential & pot)
{
 // fabs(rcut) < 1e-05 used to avoid comparing doubles
 if (nb == 0 || fabs(rcut) < 1e-05) throw PluginError("gdr", "Error in calculation: Cutoff or bins have wrong value.");
 //CurrentValue() = *(gdr(con,pot,nb,rcut));
 gdr(con,pot,nb,rcut,CurrentValue());
 #warning "GDR: memory leak?"
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Gdr(args); }
void destroy(Plugin * m) { delete m; }

