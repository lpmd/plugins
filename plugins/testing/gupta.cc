//Gupta Potential
//Cleri, F., and Rosato, F., 1993, Phys. Rev. B, 48, 22. 30
//

#include "gupta.h"
#include <iostream>

using namespace lpmd;

Gupta::Gupta(std::string args): Plugin("gupta", "2.1")
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
 std::cout << "      The plugin incorporate the Gupta potencial, used frequently              \n";
 std::cout << " for metallic atomic interaction. Based in embedded atom model.                \n\n";
 std::cout << " V(r) = A*exp(-p*(r-r0)/r0) ; F(rho) = -B*sqrt(rhoi) ;                         \n";
 std::cout << " F(rho) = exp(-2*qij*(r-r0)/r0)                                                \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      A             : Value of the A variable.                                 \n";
 std::cout << "      r0            : Value of r0, first neighbours distance.                  \n";
 std::cout << "      p             : Value of the p variable.                                 \n";
 std::cout << "      B             : Value of the B variable.                                 \n";
 std::cout << "      qij           : Value of the qij variable.                               \n";
 std::cout << "      cutoff        : Cutoff of the potential.                                 \n";
 std::cout << '\n';
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin                                                           \n";
 std::cout << " use gupta as Gup                                                              \n";
 std::cout << "     A 0.2061                                                                  \n";
 std::cout << "     r0 2.884                                                                  \n";
 std::cout << "     p 10.229                                                                  \n";
 std::cout << "     B 1.79                                                                    \n";
 std::cout << "     qij 4.036                                                                 \n";
 std::cout << "     cutoff 4.1                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Using the loaded plugin                                                      \n";
 std::cout << " potential Gup Au Au                                                          \n\n";
 std::cout << " #Real values from: Fabrizio Cleri and Vittorio Rosato PRB 48, pag. 22 (1993). \n";
}

double Gupta::pairEnergy(const double &r) const { return A*exp(-p*(r-r0)/r0); }

double Gupta::rhoij(const double &r) const { return exp(-2*qij*(r-r0)/r0); }

double Gupta::F(const double &rhoi) const { return -B*sqrt(rhoi); }

Vector Gupta::PairForce(const Vector &normrij, const double & rmod) const
{
 return -(A*p/r0)*exp(-p*(rmod-r0)/r0)*normrij;
}

Vector Gupta::ManyBodies(const Vector &normrij, const double &rhoi, const double &rhoj, const double &rmod) const
{
 double tmp=(B*qij/r0)*((1.0e0/sqrt(rhoi))+(1.0e0/sqrt(rhoj)))*exp(-2*qij*(rmod-r0)/r0);
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

