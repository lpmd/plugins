//
//
//

#include <iostream>

#include "buckingham.h"

using namespace lpmd;

Buckingham::Buckingham(std::string args): Module("buckingham")
{
 cutoff = 0.0e0;
 ProcessArguments(args);
}

double Buckingham::pairEnergy(const double & r) const
{
 if (r > cutoff) return 0.0;
 double T1 = B1*exp(-r/Ro);
 double r6 = r*r*r*r*r*r;
 double T2 = B2/r6;
 return T1-T2; 
}

Vector Buckingham::pairForce(const Vector & r) const 
{
 double rr2= r.Mod2();
 if(rr2 > cutoff*cutoff) return Vector(0.0);
 double rr = sqrt(rr2);
 double r8 = rr*rr*rr*rr*rr*rr*rr*rr;
 double f6 = 6*B2/r8;
 double pe = B1*exp(-rr/Ro)/(Ro*rr);
 double ff = f6 - pe;
 Vector fv = r;
 fv.Scale(ff);
 return fv;
}

void Buckingham::SetParameter(std::string name)
{
 if (name == "B1") 
 {
  AssignParameter("B1", GetNextWord());
  B1 = GetDouble("B1");
 }
 if (name == "Ro")
 {
  AssignParameter("Ro", GetNextWord());
  Ro = GetDouble("Ro");
 }
 if (name == "B2") 
 {
  AssignParameter("B2", GetNextWord());
  B2 = GetDouble("B2");
 }
 if (name == "cutoff")
   {
     AssignParameter("cutoff", GetNextWord());
     cutoff = GetDouble("cutoff");
   }
}

void Buckingham::Show() const
{
 Module::Show();
 std::cout << "   B1 = " << B1 << '\n';
 std::cout << "   B2 = " << B2 << '\n';
 std::cout << "   Ro = " << Ro << '\n';
 std::cout << "   cutoff = " << cutoff << '\n';
}

std::string Buckingham::Keywords() const { return "B1 B2 Ro cutoff"; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Buckingham(args); }
void destroy(Module * m) { delete m; }
