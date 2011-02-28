//Finnis-Sinclair Potential
//Finnis, M. W., and Sinclair, J. E., 1984, Philos. Mag. A, 50, 45. 3, 29, 30
//

#include "finnissinclair.h"
#include <iostream>

using namespace lpmd;

FinnisSinclair::FinnisSinclair(std::string args): Plugin("finnissinclair", "2.1")
{
 ParamList & params = (*this);
 ProcessArguments(args); 
 c0 = params["c0"];
 c1 = params["c1"];
 c2 = params["c2"];
 A = params["A"];
 B = params["B"];
 c = params["c"];
 d = params["d"];
}

void FinnisSinclair::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin incorporate the FinnisSinclair potential used frequently for \n";
 std::cout << " metallic atomic interaction. Based in embedded atom model.                    \n\n";
 std::cout << " V(r) = (r-c)^2 * (c0+c1*r+c2*r^2) ; rho(r) = (r-d)^2+(B/d)*(r-d)^3 ;          \n";
 std::cout << " F(rho) = -A/sqrt(rho)                                                         \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      c0            : Value of c0 for the potential.                           \n";
 std::cout << "      c1            : Value of c1 for the potential.                           \n";
 std::cout << "      c2            : Value of c2 for the potential.                           \n";
 std::cout << "      A             : Value of  A for the potential.                           \n";
 std::cout << "      B             : Value of  B for the potential.                           \n";
 std::cout << "      c             : First cutoff of the system (cutoff_1=c).                 \n";
 std::cout << "      d             : Second cutoff od the system (cutoff_2=d).                \n";
 std::cout << " Example          >>                                                           \n";
 std::cout << " #Loading the plugin                                                           \n";
 std::cout << " use FinnisSinclair as FS                                                      \n";
 std::cout << "     c0 3.4                                                                    \n";
 std::cout << "     c1 2.0                                                                    \n";
 std::cout << "     c2 2.9                                                                    \n";
 std::cout << "     A 0.05                                                                    \n";
 std::cout << "     B 1.90                                                                    \n";
 std::cout << "     c 0.05                                                                    \n";
 std::cout << "     d 1.90                                                                    \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #using the loaded plugin                                                      \n";
 std::cout << " potential FS Cu Au                                                          \n\n";
}

double FinnisSinclair::pairEnergy(const double &r) const
{
	if(r<c) return (r-c)*(r-c)*(c0+c1*r+c2*r*r);
    else return 0.0e0;
}

double FinnisSinclair::rhoij(const double &r) const
{
	if(r<d) return ((r-d)*(r-d))*(1.0e0+B*(r-d)/d);
    else return 0.0e0;
}

double FinnisSinclair::F(const double &rhoi) const
{
	return -A*sqrt(rhoi);
}

Vector FinnisSinclair::PairForce(const Vector &normrij, const double &r) const
{
    if(r<c)
    {
	 double t1=2.0e0*(r-c)*(c0+c1*r+c2*r*r);
	 double t2=(r-c)*(r-c)*(c1+2.0e0*c2*r);
	 return (t1+t2)*normrij;
    }
    else return Vector(0.0,0.0,0.0);
}

Vector FinnisSinclair::ManyBodies(const Vector &normrij, const double &rhoi, const double &rhoj, const double &r) const
{
    if(r<d)
    {
	double t1=((1.0e0/sqrt(rhoi))+(1.0e0/sqrt(rhoj)))*A/2.0e0;
	double t2=2.0e0*(r-d)+3.0e0*B*(r-d)*(r-d)/d;
	return -t1*t2*normrij;
    }
    else return Vector(0.0e0,0.0e0,0.0e0);
}

// No longer-ranged corrections apply beyond cutoffs c and d (neither for deltarhoi, for deltaU1, for deltaU2 nor for virial)
double FinnisSinclair::deltarhoi(const double &rhobar, const int &N) const
{
 assert(&rhobar != 0);//icc 869
 assert(&N != 0);//icc 869
 return 0.0;
}
double FinnisSinclair::deltaU1(const double &rhobar, const int &N) const
{
 assert(&rhobar != 0);//icc 869
 assert(&N != 0);//icc 869
 return 0.0;
}

// Esto se incluye para que el modulo pueda ser cargado dinámicamente
Plugin * create(std::string args) {return new FinnisSinclair(args);}
void destroy(Plugin * m) { delete m; }

