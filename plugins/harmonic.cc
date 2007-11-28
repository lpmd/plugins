//
//
//

#include <iostream>

#include "harmonic.h"

using namespace lpmd;

Harmonic::Harmonic(std::string args): Module("harmonic") 
{ 
 k = 0.0e0;
 a = 0.0e0;
 ProcessArguments(args);
}

double Harmonic::pairEnergy(const double & r) const
{
 return 0.5e0*k*(fabs(r)-a)*(fabs(r)-a);
}

Vector Harmonic::pairForce(const Vector & r) const
{
 double rr = r.Mod();
 double ff = k*(rr-a)/rr;
 Vector fv = r;
 fv.Scale(ff);
 return fv;
}

//
//
//
void Harmonic::SetParameter(std::string name)
{
 if (name == "k") 
 {
  AssignParameter("k", GetNextWord());
  k = GetDouble("k");
 }
 if (name == "a")
 {
  AssignParameter("a", GetNextWord());
  a = GetDouble("a");
 }
}

void Harmonic::Show() const
{
 Module::Show();
 std::cout << "   k = " << k << '\n';
 std::cout << "   a = " << a << '\n';
}

std::string Harmonic::Keywords() const { return "k a cutoff"; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Harmonic(args); }
void destroy(Module * m) { delete m; }


