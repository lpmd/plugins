//
//
//

#include "buckingham.h"

#include <iostream>

using namespace lpmd;

Buckingham::Buckingham(std::string args): Module("buckingham")
{
 ProcessArguments(args);
 B1 = GetDouble("B1");
 B2 = GetDouble("B2");
 Ro = GetDouble("Ro");
 cutoff = GetDouble("cutoff");
 SetCutoff(cutoff);
}

void Buckingham::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = buckingham                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de buckingham SIN PARTE COULOMBIANA,   \n";
 std::cout << " para interaccion de pares.                                                    \n";
 std::cout << "      Se utiliza la pairpotential de la API para llevar a cabo el calculo.     \n";
 std::cout << "      Si se desea incluir la parte coulombiana, se podra hacer sumando otro    \n";
 std::cout << " potencial a la interaccion de las especies, en este caso ewald, que aun no    \n";
 std::cout << " esta disponible en version estable para lpmd.                                 \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      Ro            : Especifica el valor de rho para el potencial.            \n";
 std::cout << "      B1            : Especifica el valor para la constante B1 del potencial.  \n";
 std::cout << "      B2            : Especifica el valor para la constante B2 del potencial.  \n";
 std::cout << "      cutoff        : Radio de corte para el potencial.                        \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use buckingham as BK1                                                         \n";
 std::cout << "     Ro 1.0                                                                    \n";
 std::cout << "     B1 2.0                                                                    \n";
 std::cout << "     B2 1.0                                                                    \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential BK1 Ge O                                                          \n\n";
 std::cout << "      De esta forma seteamos el potencial de buckingham entre los atomos de Ge \n";
 std::cout << " y O con las constantes usadas en BK1.                                         \n";
 std::cout << "      Note la ventaja de definir mas de un tipo de potencial para un par de    \n";
 std::cout << " especies atomicas.                                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Buckingham::Keywords() const { return "B1 B2 Ro cutoff"; }

double Buckingham::pairEnergy(const double & r) const
{
 double T1 = B1*exp(-r/Ro);
 double r6 = r*r*r*r*r*r;
 double T2 = B2/r6;
 return T1-T2; 
}

Vector Buckingham::pairForce(const Vector & r) const 
{
 double rr2= r.Mod2();
 double rr = sqrt(rr2);
 double r8 = rr*rr*rr*rr*rr*rr*rr*rr;
 double f6 = 6*B2/r8;
 double pe = B1*exp(-rr/Ro)/(Ro*rr);
 double ff = f6 - pe;
 Vector fv = r;
 fv.Scale(ff);
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Buckingham(args); }
void destroy(Module * m) { delete m; }
