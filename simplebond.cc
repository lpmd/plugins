//
//
//

#include "simplebond.h"

#include <iostream>

using namespace lpmd;

SimpleBond::SimpleBond(std::string args): Plugin("simplebond", "2.1")
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = simplebond                                               \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to implement the SimpleBond interatomic pair         \n";
 std::cout << " potential.                                                                    \n";
 std::cout << " V(r) = -e*(1-0.5*(r-r0)^2)*exp(((r-r0)^2)/2*s^2)                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      epsilon       : Set the energy bond, epsilon (e).       [eV]             \n";
 std::cout << "      sigma         : Set the well widht value sigma (s).     [A]              \n";
 std::cout << "      r0            : Set the bond lenght (r0).               [A]              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use simplebond as sb                                                          \n";
 std::cout << "     sigma 0.1                                                                 \n";
 std::cout << "     r0    3.6                                                                 \n";
 std::cout << "     epsilon 0.001                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Applying the plugin :                                                         \n";
 std::cout << " potential sb Ar Ar                                                            \n";
 std::cout << "      With this, we set the SimpleBond potential between the Ar atoms types in \n";
 std::cout << " the simulation.                                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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

