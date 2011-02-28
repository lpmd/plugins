//Sutton-Chen PseudoPotential plugin.
//Sutton, A. P., and Chen, J., 1990, Philos. Mag. Lett., 61, 139. 3, 30
//Rafii-Tabar, H., and Sutton, A. P., 1991, Philos. Mag. Lett., 63, 217. 3, 30, 36, 37
//Todd, B. D., and Lynden-Bell, R. M., 1993, Surf. Science, 281, 191. 3, 30

#include "suttonchen.h"
#include <iostream>

using namespace lpmd;

SuttonChen::SuttonChen(std::string args): Plugin("suttonchen", "2.1")
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
 e = double(params["e"]);
 a = double(params["a"]);
 n = double(params["n"]);
 m = double(params["m"]);
 c = double(params["c"]);
 rcut = double(params["cutoff"]);
 an = pow(a,n);
 am = pow(a,m);
}

void SuttonChen::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin incorporate the SuttonChen potencial, used frequently         \n";
 std::cout << " for metallic atomic interaction. Based in embedded atom model.                \n\n";
 std::cout << " V(r) = e*(a/r)^n ; F(rho) = -c*e*sqrt(rhoi) ; rho(r) = (a/r)^m                \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      e             : The value of epsilon in the potential.                   \n";
 std::cout << "      a             : the value of the A constant in the potencial.            \n";
 std::cout << "      n             : Value of n in the potential.                             \n";
 std::cout << "      m             : Value of m in the potential.                             \n";
 std::cout << "      c             : Value of the c constant of the potential.                \n";
 std::cout << "      cutoff        : Cutoff of the interatomic potential.                     \n";
 std::cout << '\n';
 std::cout << " Example          >>                                                           \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use suttonchen as sc                                                          \n";
 std::cout << "     e 3.4                                                                     \n";
 std::cout << "     a 2.0                                                                     \n";
 std::cout << "     n 2.9                                                                     \n";
 std::cout << "     m 0.05                                                                    \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Using the loaded plugin :                                                    \n";
 std::cout << " potential sc Cu Cu                                                          \n\n";
}

double SuttonChen::pairEnergy(const double &r) const { double ir = 1/r; return e*an*pow(ir,n); }

double SuttonChen::rhoij(const double &r) const { double ir = 1/r; return am*pow(ir,m); }

double SuttonChen::F(const double &rhoi) const { return -c*e*sqrt(rhoi); }

Vector SuttonChen::PairForce(const Vector &normrij, const double &rmod) const
{
 double ir = 1.0e0/rmod;
 return -n*e*an*pow(ir,n+1)*normrij;
}

Vector SuttonChen::ManyBodies(const Vector &normrij, const double &rhoi, const double &rhoj, const double &rmod) const
{
 double ir = 1.0e0/rmod;
 double tmp=0.5*m*c*e*((1/sqrt(rhoi))+1/(sqrt(rhoj)))*am*pow(ir,m+1);
 return tmp*normrij;
}

// Esto se inlcuye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) {return new SuttonChen(args);}
void destroy(Plugin * m) { delete m; }

