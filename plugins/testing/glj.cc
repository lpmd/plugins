//
//
//

#include "glj.h"

#include <iostream>

using namespace lpmd;

GenericLennardJones::GenericLennardJones(std::string args): Plugin("glj", "1.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("sigma");
 DefineKeyword("epsilon");
 DefineKeyword("a", "4.0");
 DefineKeyword("b", "0.0");
 DefineKeyword("n", "12.0");
 DefineKeyword("m", "6.0");
 DefineKeyword("cutoff");
 ProcessArguments(args);
 sigma = params["sigma"];
 epsilon = params["epsilon"];
 a = params["a"];
 b = params["b"];
 n = params["n"];
 m = params["m"];
 cutoff = params["cutoff"];
 SetCutoff(cutoff);
}

void GenericLennardJones::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa un potencial generico de tipo LenardJones para inter-\n";
 std::cout << "accion de pares.                                                               \n";
 std::cout << "      Se utiliza la pairpotential de la API para llevar a cabo el calculo.     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      sigma         : Especifica el valor de sigma para el potencial.          \n";
 std::cout << "      epsilon       : Especifica el valor para epsilon del potencial.          \n";
 std::cout << "      a             : Especifica el valor para a del potencial.                \n";
 std::cout << "      b             : Especifica el valor para b del potencial.                \n";
 std::cout << "      n             : Especifica el valor para n del potencial.                \n";
 std::cout << "      m             : Especifica el valor para m del potencial.                \n";
 std::cout << "      cutoff        : Radio de corte para el potencial.                        \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use glj as LJ1                                                                \n";
 std::cout << "     sigma 3.4                                                                 \n";
 std::cout << "     epsilon 2.0                                                               \n";
 std::cout << "     a 4.0                                                                     \n";
 std::cout << "     b 0.0                                                                     \n";
 std::cout << "     n 12.0                                                                    \n";
 std::cout << "     m 6.0                                                                     \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential LJ1 Ar Ar                                                         \n\n";
 std::cout << "      De esta forma seteamos el potencial de lennardjones generico entre los   \n";
 std::cout << " atomos de Ar con las constantes seteadas en LJ1.                              \n";
}

double GenericLennardJones::pairEnergy(const double & r) const
{
 double rtmp=sigma/r;
 return a*epsilon*(pow(rtmp, n) - b*pow(rtmp, m)); 
}

Vector GenericLennardJones::pairForce(const Vector & r) const
{
 double rr2 = r.SquareModule();
 double rtmp=sigma/sqrt(rr2);
 double ff = -a*(epsilon/rr2)*(n*pow(rtmp, n) - m*b*pow(rtmp, m));
 Vector fv = r*ff;
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new GenericLennardJones(args); }
void destroy(Plugin * m) { delete m; }

