//
//
//

#include "suttonchen.h"

#include <iostream>

using namespace lpmd;

SuttonChen::SuttonChen(std::string args): Plugin("suttonchen", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("e");
 DefineKeyword("a");
 DefineKeyword("n");
 DefineKeyword("m");
 DefineKeyword("c");
 DefineKeyword("cutoff");
 ProcessArguments(args); 
 e = params["e"];
 a = params["a"];
 n = params["n"];
 m = params["m"];
 c = params["c"];
 rcut = params["cutoff"];
}

void SuttonChen::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de SuttonChen para interaccion de      \n";
 std::cout << " atomo metalicos.                                                              \n";
 std::cout << "      Se utiliza la metalpotential de la API para llevar a cabo el calculo.    \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      e             : Especifica el valor de epsilon para el potencial.        \n";
 std::cout << "      a             : Especifica el valor para la constante a del potencial.   \n";
 std::cout << "      n             : Especifica el valor para la constante n del potencial.   \n";
 std::cout << "      m             : Especifica el valor para la constante m del potencial.   \n";
 std::cout << "      c             : Especifica el valor para la constante c del potencial.   \n";
 std::cout << "      cutoff        : Radio de corte para el potencial.                        \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use suttonchen as SC                                                          \n";
 std::cout << "     e 3.4                                                                     \n";
 std::cout << "     a 2.0                                                                     \n";
 std::cout << "     n 2.9                                                                     \n";
 std::cout << "     m 0.05                                                                    \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential SC Cu Cu                                                          \n\n";
 std::cout << "      De esta forma seteamos el potencial de suttonchen entre los atomos de    \n";
 std::cout << " Cu con las constantes seteadas en SC.                                         \n";
}

double SuttonChen::pairEnergy(const double &r) const { return e*pow((a/r),n); }

double SuttonChen::rhoij(const double &r) const { return pow((a/r),m); }

double SuttonChen::F(const double &rhoi) const { return -c*e*sqrt(rhoi); }

Vector SuttonChen::PairForce(const Vector &rij) const
{
 Vector norm = rij;
 double rmod = rij.Module();
 norm.Normalize();
 return -n*e*pow(a/rmod,n)*(norm/rmod);
}

Vector SuttonChen::ManyBodies(const Vector &rij, const double &rhoi, const double &rhoj) const
{
 double tmp;
 double rmod = rij.Module();
 tmp=(m/2)*c*e*((1/sqrt(rhoi))+(1/sqrt(rhoj)))*pow(a/rmod,m)*(1.0/rmod);
 Vector ff = rij;
 ff.Normalize();
 return tmp*ff;
}

double SuttonChen::deltarhoi(const double &rhobar) const
{
 return (4*M_PI*rhobar*a*a*a/(m-3))*pow(a/rcut,m-3);
}

double SuttonChen::deltaU1(const double &rhobar, const int &N) const
{
 double f = 2*M_PI*N*rhobar*e*a*a*a/(n-3);
 return f*pow(a/rcut,n-3);
}

double SuttonChen::deltaU2(const double &rhobar, const int &N, const double &rhoi) const { return 0.0e0; }

double SuttonChen::VirialCorrection(const double &rhobar, const int &N, const double &rhoi) const
{
 double dV1 = (n/(n-3))*pow(a/rcut, n-3);
 double dV2 = (2.0*m/(m-3))*pow(a/rcut, m-3)*c/(2.0*sqrt(rhoi));
 return -(2.0*M_PI*rhobar*N*e*pow(a, 3.0)*(dV1-dV2));
}

// Esto se inlcuye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) {return new SuttonChen(args);}
void destroy(Plugin * m) { delete m; }

