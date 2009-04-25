//
//
//

#include "energy.h"

#include <lpmd/simulationcell.h>
#include <lpmd/potential.h>

using namespace lpmd;

Energy::Energy(std::string args): Module("energy") 
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 ProcessArguments(args); 
}

Energy::~Energy() { }

void Energy::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para calcular la energia cinetica y potencial     \n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use energy                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " -----------------------------                                               \n\n";
 std::cout << "      De esta forma activa el calculo de las energias utilizando energy        \n";
}

std::string Energy::Provides() const 
{ 
 return "kinetic-energy potential-energy total-energy momentum px py pz temperature"; 
}

double Energy::GetProperty(const std::string & name)
{
 if (name == "kinetic-energy") return ekin;
 if (name == "potential-energy") return epot;
 if (name == "total-energy") return etot;
 if (name == "temperature") return temp;
 if (name == "momentum") return pv.Mod();
 if (name == "px") return pv.Get(0);
 if (name == "py") return pv.Get(1);
 if (name == "pz") return pv.Get(2);
 throw UnknownProperty(name);
}

void Energy::Evaluate(SimulationCell & sc, Potential & pot)
{
 epot = pot.energy(sc);
 ekin = sc.KineticEnergy();
 etot = epot+ekin;
 temp = sc.Temperature();
 pv.Zero();
 for (unsigned long int i=0;i<sc.size();++i)
 {
  pv = pv + sc[i].Velocity()*sc[i].Mass();
 }
 pv = pv / double(sc.size());
}

const double & Energy::CurrentValue() const { return etot; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Energy(args); }
void destroy(Module * m) { delete m; }

