//
//
//

#include "gupta.h"

#include <iostream>

using namespace lpmd;

Gupta::Gupta(std::string args): Plugin("gupta", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("A");
 DefineKeyword("r0");
 DefineKeyword("p");
 DefineKeyword("B");
 DefineKeyword("qij");
 DefineKeyword("cutoff");
 ProcessArguments(args); 
 A = params["A"];
 r0 = params["r0"];
 p = params["p"];
 B = params["B"];
 qij = params["qij"];
 rcut = params["cutoff"];
}

void Gupta::ShowHelp() const
{
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
 std::cout << '\n';
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
}

double Gupta::pairEnergy(const double &r) const { return A*exp(-p*(r-r0)/r0); }

double Gupta::rhoij(const double &r) const { return exp(-2*qij*(r-r0)/r0); }

double Gupta::F(const double &rhoi) const { return -B*sqrt(rhoi); }

Vector Gupta::PairForce(const Vector &normrij, const double & rmod) const
{
 //FIXME : Chequear bien potencial.
// Vector norm = rij;
// double rmod = rij.Module();
// norm.Normalize();
 return -(A*p/r0)*exp(-p*(rmod-r0)/r0)*(normrij/rmod);
}

Vector Gupta::ManyBodies(const Vector &normrij, const double &invrhoi, const double &invrhoj, const double &rmod) const
{
 //FIXME : Chequear bien el potential
// double tmp;
// double rmod = rij.Module();
 double tmp=(B*qij/r0)*(invrhoi+invrhoj)*exp(-2*qij*(rmod-r0)/r0)*(1.0/rmod);
// Vector ff = rij;
// ff.Normalize();
 return tmp*normrij;
}

double Gupta::deltarhoi(const double &rhobar) const
{
 double tmp =r0/qij;
 return (2*M_PI*rhobar*tmp)*(rcut*rcut+2*rcut*tmp+2*tmp*tmp)*exp(-2*(rcut-r0)/tmp);
}

double Gupta::deltaU1(const double &rhobar, const int &N) const
{
 double tmp = r0/p;
 double f = 2*M_PI*N*rhobar*A*tmp;
 return f*(rcut*rcut+2*rcut*tmp+2*tmp*tmp)*exp(-(rcut-r0)/tmp);
}

double Gupta::deltaU2(const double &rhobar, const int &N, const double &rhoi) const
{
 assert(&rhobar != 0);//icc869
 assert(&N != 0);//icc 869
 assert(&rhoi !=0);//icc 869
 return 0.0e0;
}

double Gupta::VirialContribution(const double &r, const double &rhoi, const double &rhoj) const 
{ 
 assert(&r != 0); // icc 869
 assert(&rhoi != 0);//icc 869
 assert(&rhoj != 0);//icc 869
 return 0.0; 
}

double Gupta::VirialCorrection(const double &rhobar, const int &N, const double &rhoi) const
{
 double tmp = r0/p;
 double tmq = r0/qij;
 double dV1 = N*A*(rcut*rcut*rcut+3*rcut*rcut*tmp+6*rcut*tmp*tmp+6*tmp*tmp*tmp)*exp(-(rcut-r0)/tmp);
 double dV2 = (rcut*rcut*rcut+3*rcut*rcut*tmq+6*rcut*tmq*tmq+6*tmq*tmq*tmq)*exp(-2*(rcut-r0)/tmq)*(N*B/(2.0*sqrt(rhoi)));
 return -(2.0*M_PI*rhobar*(dV1-dV2));
}

// Esto se inlcuye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) {return new Gupta(args);}
void destroy(Plugin * m) { delete m; }

