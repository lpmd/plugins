//
//
//

#include "gupta.h"

#include <iostream>

using namespace lpmd;

Gupta::Gupta(std::string args): Module("gupta")
{
 ProcessArguments(args); 
 A = GetDouble("A");
 r0 = GetDouble("r0");
 p = GetDouble("p");
 B = GetDouble("B");
 qij = GetDouble("qij");
 rcut = GetDouble("cutoff");
}

void Gupta::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = gupta                                                    \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de Gupta para interaccion de           \n";
 std::cout << " atomos metalicos.                                                             \n";
 std::cout << "      Se utiliza la metalpotential de la API para llevar a cabo el calculo.    \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      A             : Especifica el valor para la constante A del potencial.   \n";
 std::cout << "      r0            : Especifica la distancia a primeros vecinos.              \n";
 std::cout << "      p             : Especifica el valor para la constante p del potencial.   \n";
 std::cout << "      B             : Especifica el valor para la constante B del potencial.   \n";
 std::cout << "      qij           : Especifica el valor para la constante qij del potencial. \n";
 std::cout << "      cutoff        : Radio de corte para el potencial.                        \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use gupta as Gup                                                              \n";
 std::cout << "     A 0.2061                                                                  \n";
 std::cout << "     r0 2.884                                                                  \n";
 std::cout << "     p 10.229                                                                  \n";
 std::cout << "     B 1.79                                                                    \n";
 std::cout << "     qij 4.036                                                                 \n";
 std::cout << "     cutoff 4.1                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential Gup Au Au                                                          \n\n";
 std::cout << "      De esta forma seteamos el potencial de gupta entre los atomos de Au      \n";
 std::cout << " Constantes obtenidas de ref. Fabrizio Cleri and Vittorio Rosato PRB 48, pag. 22 (1993).\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Gupta::Keywords() const { return " A r0 p B qij cutoff "; }

double Gupta::pairEnergy(const double &r) const
{
 return A*exp(-p*(r-r0)/r0);
}

double Gupta::rhoij(const double &r) const
{
 return exp(-2*qij*(r-r0)/r0);
}

double Gupta::F(const double &rhoi) const
{
 return -B*sqrt(rhoi);
}

Vector Gupta::PairForce(const Vector &rij) const
{
 Vector norm = rij;
 double rmod = rij.Mod();
 norm.Norm();
 return -(A*p/r0)*exp(-p*(rmod-r0)/r0)*(norm/rmod);
}

Vector Gupta::ManyBodies(const Vector &rij, const double &rhoi, const double &rhoj) const
{
 double tmp;
 double rmod = rij.Mod();
 tmp=(B*qij/r0)*((1/rhoi)+(1/rhoj))*exp(-2*qij*(rmod-r0)/r0)*(1.0/rmod);
 Vector ff = rij;
 ff.Norm();
 return tmp*ff;
}

double Gupta::deltarhoi(const double &rhobar, const int &N) const
{
 return 0;
}

double Gupta::deltaU1(const double &rhobar, const int &N) const
{
 return 0;
}

double Gupta::deltaU2(const double &rhobar, const int &N, const double &rhoi) const
{
 return 0.0e0;
}

double Gupta::VirialContribution(const double &r, const double &rhoi, const double &rhoj) const { return 0.0; }

double Gupta::VirialCorrection(const double &rhobar, const int &N, const double &rhoi) const
{
 return 0;
}

// Esto se inlcuye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) {return new Gupta(args);}
void destroy(Module * m) { delete m; }

