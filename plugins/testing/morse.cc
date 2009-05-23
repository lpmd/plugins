//
//
//

#include "morse.h"

#include <iostream>

using namespace lpmd;

Morse::Morse(std::string args): Module("morse") 
{ 
 ParamList & params = (*this);
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("depth");
 DefineKeyword("a");
 DefineKeyword("re");
 ProcessArguments(args);
 depth = params["depth"];
 a = params["a"];
 re = params["re"];
}

void Morse::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de Morse para la interaccion de        \n";
 std::cout << " de pares.                                                                     \n";
 std::cout << "      Se utiliza pairpotential de la API para llevar a cabo el calculo.        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      depth         : Especifica la profundidad del pozo de potencial.         \n";
 std::cout << "      a             : Especifica el valor de a (ancho)para el potencial.       \n";
 std::cout << "      re            : Especifica valor de re (distancia eq) para el potencial. \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use morse as MO                                                               \n";
 std::cout << "     depth 3.4                                                                 \n";
 std::cout << "     a 3.0                                                                     \n";
 std::cout << "     re 1.4                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential MO Ar Ar                                                          \n\n";
 std::cout << "      De esta forma seteamos el potencial de Morse entre los atomos            \n";
 std::cout << " de Ar con las constantes usadas en MO.                                        \n";
}

double Morse::pairEnergy(const double & r) const
{
 return depth*(1.0-exp(-a*(r - re)))*(1.0-exp(-a*(r - re)));
}

Vector Morse::pairForce(const Vector & r) const
{
 double rr = r.Module();
 // FIXME: Checkear la expresion para la fuerza, signo correcto y esas cosas
 double ff = 2.0*a*(depth/rr)*(1.0-exp(-a*(rr-re)))*exp(-a*(rr-re));
 Vector fv = r * ff;
 //fv.Scale(ff);
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Morse(args); }
void destroy(Module * m) { delete m; }
