//
//
//

#include <iostream>

#include "morse.h"

using namespace lpmd;

Morse::Morse(std::string args): Module("morse") 
{ 
 ProcessArguments(args);
}

double Morse::pairEnergy(const double & r) const
{
 return depth*(1.0-exp(-a*(r - re)));
}

Vector Morse::pairForce(const Vector & r) const
{
 double rr = r.Mod();
 // FIXME: Checkear la expresion para la fuerza, signo correcto y esas cosas
 double ff = 2.0*a*(depth/rr)*(1.0-exp(-a*(rr-re)))*exp(-a*(rr-re));
 Vector fv = r;
 fv.Scale(ff);
 return fv;
}

//
//
//
void Morse::SetParameter(std::string name)
{
 if (name == "depth") 
 {
  AssignParameter("depth", GetNextWord());
  depth = GetDouble("depth");
 }
 if (name == "a")
 {
  AssignParameter("a", GetNextWord());
  a = GetDouble("a");
 }
 if (name == "re") 
 {
  AssignParameter("re", GetNextWord());
  re = GetDouble("re");
 }
}

void Morse::Show() const
{
 Module::Show();
 std::cout << "   depth = " << depth << '\n';
 std::cout << "       a = " << a << '\n';
 std::cout << "      re = " << re << '\n';
}

std::string Morse::Keywords() const { return "depth a re"; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Morse(args); }
void destroy(Module * m) { delete m; }


