//
//
//

#include "simplebond.h"

#include <iostream>

using namespace lpmd;

SimpleBond::SimpleBond(std::string args): Plugin("simplebond", "1.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("epsilon");
 DefineKeyword("sigma");
 DefineKeyword("r0");
 ProcessArguments(args);
 epsilon = params["epsilon"];
 sigma = params["sigma"];
 r0 = params["r0"];
}

void SimpleBond::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial SimpleBond para la interaccion de       \n";
 std::cout << " de pares.                                                                     \n";
 std::cout << "      Se utiliza pairpotential de la API para llevar a cabo el calculo.        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      epsilon       : Especifica la energia de enlace.                         \n";
 std::cout << "      sigma         : Especifica el ancho del pozo.                            \n";
 std::cout << "      r0            : Especifica la longitud de enlace.                        \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use simplebond as sb                                                          \n";
 std::cout << "     sigma 0.1                                                                 \n";
 std::cout << "     r0    3.6                                                                 \n";
 std::cout << "     epsilon 0.001                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential sb Ar Ar                                                            \n\n";
 std::cout << "      De esta forma seteamos el potencial SimpleBond entre los atomos          \n";
 std::cout << " de Ar con las constantes usadas en sb.                                        \n";
}

double SimpleBond::pairEnergy(const double & r) const
{
 return -epsilon*(1.0-0.5*pow(r-r0,2.0))*exp(-pow(r-r0,2.0)/(2.0*sigma*sigma));
}

Vector SimpleBond::pairForce(const Vector & r) const
{
 double rr = r.Module();
 double ff = (1.0/rr)*epsilon*(rr-r0)*exp(-pow(rr-r0,2.0)/(2.0*sigma*sigma))*(1.0+(1.0-0.5*pow(rr-r0,2.0))/(sigma*sigma));
 Vector fv = r*ff;
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SimpleBond(args); }
void destroy(Plugin * m) { delete m; }

