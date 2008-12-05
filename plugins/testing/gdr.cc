//
//
//

#include "gdr.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>

using namespace lpmd;

Gdr::Gdr(std::string args): Module("gdr")
{
 m = NULL;
 AssignParameter("average", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = GetDouble("rcut");
 nb = GetInteger("bins");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 outputfile = GetString("output");
 do_average = GetBool("average");
}

Gdr::~Gdr() { if (m != NULL) delete m; }

void Gdr::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = gdr                                                      \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular la funcion de autocorrelacion de pares    \n";
 std::cout << " utiliza las condiciones de borde periodicas de la celda.                      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Especifica el numero de divisiones entre 0 y rcut.       \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo de gdr.       \n";
 std::cout << "      output        : Fichero en el que se graba el RDF.                       \n";
 std::cout << "      average       : Setea si calculo o no el promedio de cada calculo.       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "      De esta forma clculamos la funcion de distribucion radial de pares en    \n";
 std::cout << " la simulacion entre 1 y 100 cada 10 pasos.                                    \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Gdr::Keywords() const { return "rcut bins start end each output average"; }

void Gdr::Evaluate(SimulationCell & simcell, Potential & pot)
{
 // fabs(rcut) < 1e-05 used to avoid comparing doubles
 if (nb == 0 || fabs(rcut) < 1e-05) throw PluginError("gdr", "Error in calculation: Cutoff or bins have wrong value.");
 if (m != NULL) delete m;
 m=gdr(simcell,pot,nb,rcut);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Gdr(args); }
void destroy(Module * m) { delete m; }


