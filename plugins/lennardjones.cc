//
//
//

#include <iostream>

#include "lennardjones.h"

using namespace lpmd;

LennardJones::LennardJones(std::string args): Module("lennardjones") 
{ 
 ProcessArguments(args);
}

double LennardJones::pairEnergy(const double & r) const
{
 // std::cerr << "Calculando energia LJ para r=" << r << std::endl; 
 if (r > cutoff) return 0.0;
 double rtmp=sigma/r;
 double r6 = rtmp*rtmp*rtmp*rtmp*rtmp*rtmp;
 double r12 = r6*r6;
 return 4.0e0*epsilon*(r12 - r6); 
}

Vector LennardJones::pairForce(const Vector & r) const
{
 // std::cerr << "Calculando fuerza LJ para vector r=" << r << std::endl; 
 double rr2 = r.Mod2();
 if (rr2 > cutoff*cutoff) return Vector(0.0);
 double r6 = pow(sigma*sigma / rr2, 3.0e0);
 double r12 = r6*r6;
 double ff = -48.0e0*(epsilon/rr2)*(r12 - 0.50e0*r6);
 Vector fv = r;
 fv.Scale(ff);
 return fv;
}

//
//
//
void LennardJones::SetParameter(std::string name)
{
 if (name == "sigma") 
 {
  AssignParameter("sigma", GetNextWord());
  sigma = GetDouble("sigma");
 }
 if (name == "epsilon")
 {
  AssignParameter("epsilon", GetNextWord());
  epsilon = GetDouble("epsilon");
 }
 if (name == "cutoff") 
 {
  AssignParameter("cutoff", GetNextWord());
  cutoff = GetDouble("cutoff");
 }
}

void LennardJones::Show() const
{
 Module::Show();
 std::cout << "   sigma = " << sigma << '\n';
 std::cout << "   epsilon = " << epsilon << '\n';
 std::cout << "   cutoff = " << cutoff << '\n';
}

std::string LennardJones::Keywords() const { return "sigma epsilon cutoff"; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new LennardJones(args); }
void destroy(Module * m) { delete m; }


