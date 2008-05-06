//
//
//

#include "cell.h"

#include <lpmd/physunits.h>
#include <lpmd/simulationcell.h>
#include <lpmd/potential.h>

using namespace lpmd;

CellProp::CellProp(std::string args): Module("cell") { ProcessArguments(args); }

CellProp::~CellProp() { }

void CellProp::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = cell                                                     \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para calcular dimensiones de la celda y densidades\n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cell                                                                      \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " -----------------------------                                               \n\n";
 std::cout << "      De esta forma activa el calculo de las dimensiones utilizando cell       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string CellProp::Keywords() const { return ""; }

std::string CellProp::Provides() const { return "volume cell-a cell-b cell-c density particledensity volumeperatom "; }

double CellProp::GetProperty(const std::string & name)
{
 if (name == "volume") return volume;
 if (name == "volumeperatom") return (volume/double(nat));
 if (name == "cell-a") return a;
 if (name == "cell-b") return b;
 if (name == "cell-c") return c;
 if (name == "density") return dens;
 if (name == "particledensity") return partdens;
 throw UnknownProperty(name);
}

void CellProp::Evaluate(SimulationCell & sc, Potential & pot)
{
 nat = sc.Size();
 volume = sc.Volume();
 a = sc.GetVector(0).Mod();
 b = sc.GetVector(1).Mod();
 c = sc.GetVector(2).Mod();
 dens = sc.Density();
 partdens = sc.ParticleDensity();
}

const double & CellProp::Value() const { return volume; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new CellProp(args); }
void destroy(Module * m) { delete m; }

