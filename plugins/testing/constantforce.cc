//
//
//

#include "constantforce.h"

#include <lpmd/simulationcell.h>
#include <lpmd/session.h>

#include <iostream>

using namespace lpmd;

ConstantForcePotential::ConstantForcePotential(std::string args): Module("constantforce") 
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("forcevector");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 fx = GetDouble("fx");
 fy = GetDouble("fy");
 fz = GetDouble("fz");
}

ConstantForcePotential::~ConstantForcePotential() { }

void ConstantForcePotential::SetParameter(std::string name) 
{ 
 if (name == "forcevector") 
 {
  AssignParameter("fx", GetNextWord());
  AssignParameter("fy", GetNextWord());
  AssignParameter("fz", GetNextWord());
 }
 else Module::SetParameter(name);
}

void ConstantForcePotential::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo setea una fuerza constante, es utilizado para trabajo con un   \n";
 std::cout << " set de particulas que se quieran configurar con cierta fuerza.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      forcevector   : Seguido por 3 numeros reales que indican un vector para  \n";
 std::cout << "                      la fuerza.                                               \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use constantforce as gravity                                                  \n";
 std::cout << "     forcevector 0.0 0.0 -9.8                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " potential gravity Ar Ar                                                     \n\n";
 std::cout << "      De esta forma aplicamos una fuerza constante a la interaccion de los     \n";
 std::cout << " atomos de Argon.                                                              \n";
}

double ConstantForcePotential::energy(SimulationCell & sc) { return 0.0; }

void ConstantForcePotential::UpdateForces(SimulationCell & sc) 
{ 
 const double forcefactor = GlobalSession.GetDouble("forcefactor");
 for (unsigned long int i=0;i<sc.size();++i)
 {
  const Atom & atom_i = sc[i];
  double mi = atom_i.Mass();
  Vector ff(fx, fy, fz);
  const Vector & acci = atom_i.Acceleration(); 
  sc.SetAcceleration(i, acci + ff*(forcefactor/mi));
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ConstantForcePotential(args); }
void destroy(Module * m) { delete m; }


