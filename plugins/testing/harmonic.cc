//
//
//

#include "harmonic.h"

#include <iostream>

using namespace lpmd;

Harmonic::Harmonic(std::string args): Plugin("harmonic", "2.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("k");
 DefineKeyword("a");
 ProcessArguments(args);
 k = params["k"];
 a = params["a"];
}

void Harmonic::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial Harmonico para la interaccion de       \n";
 std::cout << " de pares.                                                                     \n";
 std::cout << "      Se utiliza pairpotential de la API para llevar a cabo el calculo.        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      k             : Especifica la constante del resorte interatomico.        \n";
 std::cout << "      a             : Especifica el valor para el largo interatomico.          \n";
 std::cout << '\n';
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
}

double Harmonic::pairEnergy(const double & r) const
{
 return 0.5e0*k*(fabs(r)-a)*(fabs(r)-a);
}

Vector Harmonic::pairForce(const Vector & r) const
{
 double rr = r.Module();
 double ff = k*(rr-a)/rr;
 Vector fv = r*ff;
// fv.Scale(ff);
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Harmonic(args); }
void destroy(Plugin * m) { delete m; }
