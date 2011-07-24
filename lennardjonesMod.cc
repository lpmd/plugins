//
//
//

#include "lennardjonesMod.h"

#include <iostream>

using namespace lpmd;

LennardJonesMod::LennardJonesMod(std::string args): Plugin("lennardjonesMod", "2.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("sigma");
 DefineKeyword("epsilon");
 DefineKeyword("m", "12");
 DefineKeyword("n", "6"); 
 DefineKeyword("cohesive", "1");
 DefineKeyword("cutoff");
 ProcessArguments(args);
 sigma = params["sigma"];
 epsilon = params["epsilon"];
 m = params["m"];
 n = params["n"];
 cohesive = params["cohesive"];
 cutoff = params["cutoff"];
 SetCutoff(cutoff);
}

void LennardJonesMod::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = lennardjonesMod                                          \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The module implements a modified Lenard-Jones potential for pairs        \n";
 std::cout << "      interaction, with the form                                               \n";
 std::cout << "            V(r)=4*epsilon( (sigma/r)^m - (sigma/r)^n )                        \n"; 
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      sigma         : Determines de value of sigma for the potential.          \n";
 std::cout << "      epsilon       : Determines de value of  epsilon for the potential.       \n";
 std::cout << "      m             : Determines de value of the first exponent.               \n";
 std::cout << "      n             : EDetermines de value of the second exponent.             \n";
 std::cout << "      cohesive      : Cohesion parameter.                                      \n";
 std::cout << "      cutoff        : Cutoff for the potential.                                \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use lennardjonesMod as ljm                                                    \n";
 std::cout << "     sigma 3.4                                                                 \n";
 std::cout << "     epsilon 2.0                                                               \n";
 std::cout << "     cohesive 0.2                                                              \n";
 std::cout << "     m 4                                                                       \n";
 std::cout << "     n 2                                                                       \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";
 std::cout << " potential ljm Ar Ar                                                         \n\n";
 std::cout << "      The plugin implements the Lennard-Jones potential between argon (Ar)     \n";
 std::cout << "      atoms with little cohesion (cohesive=0 is a purely repulsive potential). \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

double LennardJonesMod::pairEnergy(const double & r) const
{
 double rtmp=sigma/r;
 double rn = pow(rtmp,n);
 double rm = pow(rtmp,m);
 return 4.0e0*epsilon*(rm - cohesive*rn); 
}

Vector LennardJonesMod::pairForce(const Vector & r) const
{
 double rr2 = r.SquareModule();
 double rn = pow(sigma*sigma / rr2, n/2.0);
 double rm = pow(sigma*sigma / rr2, m/2.0);
 double ff = -4.0e0*(epsilon/rr2)*(m*rm - cohesive*n*rn);
 Vector fv = r*ff;
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LennardJonesMod(args); }
void destroy(Plugin * m) { delete m; }

