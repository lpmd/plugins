//
//
//

#include "harmonic.h"

#include <iostream>

using namespace lpmd;

Harmonic::Harmonic(std::string args): Module("harmonic") 
{ 
 ProcessArguments(args);
 k = GetDouble("k");
 a = GetDouble("a");
}

void Harmonic::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = harmonic                                                 \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial Harmonico para la interaccion de       \n";
 std::cout << " de pares.                                                                     \n";
 std::cout << "      Se utiliza pairpotential de la API para llevar a cabo el calculo.        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      k             : Especifica la constante del resorte interatomico.        \n";
 std::cout << "      a             : Especifica el valor para el largo interatomico.          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use harmonic as HP                                                            \n";
 std::cout << "     k 3.4                                                                     \n";
 std::cout << "     a 3                                                                       \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential HP Ar Ar                                                          \n\n";
 std::cout << "      De esta forma seteamos el potencial de Harmonic entre los atomos         \n";
 std::cout << " de Ar con las constantes usadas en HP.                                        \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Harmonic::Keywords() const { return "k a cutoff"; }

double Harmonic::pairEnergy(const double & r) const
{
 return 0.5e0*k*(fabs(r)-a)*(fabs(r)-a);
}

Vector Harmonic::pairForce(const Vector & r) const
{
 double rr = r.Mod();
 double ff = k*(rr-a)/rr;
 Vector fv = r;
 fv.Scale(ff);
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Harmonic(args); }
void destroy(Module * m) { delete m; }
