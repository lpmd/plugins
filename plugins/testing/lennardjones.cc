//
//
//

#include "lennardjones.h"

#include <iostream>

using namespace lpmd;

LennardJones::LennardJones(std::string args): Module("lennardjones") 
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("sigma");
 DefineKeyword("epsilon");
 DefineKeyword("cutoff");
 ProcessArguments(args);
 sigma = GetDouble("sigma");
 epsilon = GetDouble("epsilon");
 cutoff = GetDouble("cutoff");
 SetCutoff(cutoff);
}

void LennardJones::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de LenardJones para interaccion de     \n";
 std::cout << " pares.                                                                        \n";
 std::cout << "      Se utiliza la pairpotential de la API para llevar a cabo el calculo.     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      sigma         : Especifica el valor de sigma para el potencial.          \n";
 std::cout << "      epsilon       : Especifica el valor para epsilon del potencial.          \n";
 std::cout << "      cutoff        : Radio de corte para el potencial.                        \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use lennardjones as LJ1                                                       \n";
 std::cout << "     sigma 3.4                                                                 \n";
 std::cout << "     epsilon 2.0                                                               \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential LJ1 Ar Ar                                                         \n\n";
 std::cout << "      De esta forma seteamos el potencial de lennardjones entre los atomos de  \n";
 std::cout << " Ar con las constantes seteadas en LJ1.                                        \n";
}

double LennardJones::pairEnergy(const double & r) const
{
 double rtmp=sigma/r;
 double r6 = rtmp*rtmp*rtmp*rtmp*rtmp*rtmp;
 double r12 = r6*r6;
 return 4.0e0*epsilon*(r12 - r6); 
}

Vector LennardJones::pairForce(const Vector & r) const
{
 double rr2 = r.SquareModule();
 double r6 = pow(sigma*sigma / rr2, 3.0e0);
 double r12 = r6*r6;
 double ff = -48.0e0*(epsilon/rr2)*(r12 - 0.50e0*r6);
 Vector fv = r * ff;
// fv.Scale(ff);
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new LennardJones(args); }
void destroy(Module * m) { delete m; }


