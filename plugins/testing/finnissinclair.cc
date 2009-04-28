//
//
//

#include "finnissinclair.h"

#include <iostream>

using namespace lpmd;

FinnisSinclair::FinnisSinclair(std::string args): Module("finnissinclair")
{
 ProcessArguments(args); 
 c0 = GetDouble("c0");
 c1 = GetDouble("c1");
 c2 = GetDouble("c2");
 A = GetDouble("A");
 B = GetDouble("B");
 c = GetDouble("c");
 d = GetDouble("d");
}

void FinnisSinclair::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = finnissinclair                                           \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de FinnisSinclair para interaccion de  \n";
 std::cout << " atomos metalicos.                                                             \n";
 std::cout << "      Se utiliza la metalpotential de la API para llevar a cabo el calculo.    \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      c0            : Especifica el valor de c0 para el potencial.             \n";
 std::cout << "      c1            : Especifica el valor de c1 para el potencial.             \n";
 std::cout << "      c2            : Especifica el valor de c2 para el potencial.             \n";
 std::cout << "      A             : Especifica el valor de A para el potencial.              \n";
 std::cout << "      B             : Especifica el valor de beta para el potencial.           \n";
 std::cout << "      c             : Especifica primer radio de corte (cutoff_1=c).           \n";
 std::cout << "      d             : Especifica segundo radio de corte (cutoff_2=d).          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use FinnisSinclair as FS                                                      \n";
 std::cout << "     c0 3.4                                                                    \n";
 std::cout << "     c1 2.0                                                                    \n";
 std::cout << "     c2 2.9                                                                    \n";
 std::cout << "     A 0.05                                                                    \n";
 std::cout << "     B 1.90                                                                    \n";
 std::cout << "     c 0.05                                                                    \n";
 std::cout << "     d 1.90                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential FS Cu Au                                                          \n\n";
 std::cout << "      De esta forma fijamos el potencial de FinnisSinclair entre los atomos de \n";
 std::cout << " cobre y oro con las constantes seteadas en FS.                                \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

}

std::string FinnisSinclair::Keywords() const { return " c0 c1 c2 A B c d "; }

double FinnisSinclair::pairEnergy(const double &r) const
{
	return (r-c)*(r-c)*(c0+c1*r+c2*r*r);
}

double FinnisSinclair::rhoij(const double &r) const
{
	return (r-d)*(r-d)*(1.0e0+B*(r-d)/d);
}

double FinnisSinclair::F(const double &rhoi) const
{
	return -A*sqrt(rhoi);
}

Vector FinnisSinclair::PairForce(const Vector &rij) const
{
	double r = rij.Module();
	double t1=2.0e0*(r-c)*(c0+c1*r+c2*r*r);
	double t2=(r-c)*(r-c)*(c1+2.0e0*c2*r);
	return (t1+t2)*(rij/r);
}

Vector FinnisSinclair::ManyBodies(const Vector &rij, const double &rhoi, const double &rhoj) const
{
	double r = rij.Module();
	double t1=(1.0/rhoi+1.0/rhoj)*A/2.0;
	double t2=2.0*(r-d)+3.0*B*(r-d)*(r-d)/d;
	return -t1*t2*rij/r;
}

// No longer-ranged corrections apply beyond cutoffs c and d (neither for deltarhoi, for deltaU1, for deltaU2 nor for virial)
double FinnisSinclair::deltarhoi(const double &rhobar, const int &N) const{	return 0.0;}
double FinnisSinclair::deltaU1(const double &rhobar, const int &N) const{	return 0.0;}
double FinnisSinclair::deltaU2(const double &rhobar, const int &N, const double &rhoi) const { return 0.0; }

// Esto se incluye para que el modulo pueda ser cargado dinámicamente
Module * create(std::string args) {return new FinnisSinclair(args);}
void destroy(Module * m) { delete m; }

