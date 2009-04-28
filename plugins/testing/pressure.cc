//
//
//

#include "pressure.h"

#include <lpmd/simulationcell.h>
#include <lpmd/session.h>

using namespace lpmd;

Pressure::Pressure(std::string args): Module("pressure") 
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 ProcessArguments(args); 
}

Pressure::~Pressure() { }

void Pressure::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para calcular la presion en la celda.             \n";
 std::cout << " General Options   >> No Requiere                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use pressure                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " -----------------------------                                               \n\n";
 std::cout << "      De esta forma activa el calculo de la presion utilizando pressure.       \n";
}

std::string Pressure::Provides() const { return "pressure virial-pressure kinetic-pressure sxx sxy sxz syx syy syz szx szy szz"; }

double Pressure::GetProperty(const std::string & name)
{
 if (name == "pressure") return press;
 if (name == "virial-pressure") return vpress;
 if (name == "kinetic-pressure") return kpress;
 if (name == "sxx") return s[0][0];
 if (name == "sxy") return s[0][1];
 if (name == "sxz") return s[0][2];
 if (name == "syx") return s[1][0];
 if (name == "syy") return s[1][1];
 if (name == "syz") return s[1][2];
 if (name == "szx") return s[2][0];
 if (name == "szy") return s[2][1];
 if (name == "szz") return s[2][2];
 throw UnknownProperty(name);
}

void Pressure::Evaluate(SimulationCell & sc, Potential & pot)
{
 const double pressfactor = GlobalSession.GetDouble("pressfactor");
 const double kin2ev = GlobalSession.GetDouble("kin2ev");
 double v = sc.Volume();
 double K=0.0e0;
 for (unsigned long int i=0;i<sc.size();++i)
 {
  const Atom & a = sc[i];
  double m = a.Mass();
  Vector vel = a.Velocity();
  K += m*vel.SquareModule();
 }
 K = 0.5*K;
 kpress = (pressfactor/v)*(2.0/3.0)*K*kin2ev;
 vpress = (pressfactor/v)*(1.0/3.0)*sc.Virial();
 press = kpress+vpress;
 for (int i=0;i<3;i++)
 {
  s[0][i] = (pressfactor/v)*sc.StressTensor(0,i);
  s[1][i] = (pressfactor/v)*sc.StressTensor(1,i);
  s[2][i] = (pressfactor/v)*sc.StressTensor(2,i);
 }
}

const double & Pressure::CurrentValue() const { return press; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Pressure(args); }
void destroy(Module * m) { delete m; }


