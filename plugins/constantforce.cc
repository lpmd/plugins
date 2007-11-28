//
//
//

#include <iostream>

#include <lpmd/physunits.h>

#include "constantforce.h"

using namespace lpmd;

ConstantForcePotential::ConstantForcePotential(std::string args): Module("constantforce") 
{ 
 fx = fy = fz = 0.0;
 ProcessArguments(args);
}

ConstantForcePotential::~ConstantForcePotential() { }

double ConstantForcePotential::energy(SimulationCell & sc) 
{ 
 return 0.0; 
}

void ConstantForcePotential::UpdateForces(SimulationCell & sc) 
{ 
 long n = sc.Size();
 for (long i=0;i<n;++i)
 {
  const Atom & atom_i = sc.GetAtom(i);
  double mi = atom_i.Mass();
  Vector ff(fx, fy, fz);
  const Vector & acci = atom_i.Acceleration(); 
  sc.SetAcceleration(i, acci + ff*(FORCEFACTOR/mi));
 }
}

//
//
//
void ConstantForcePotential::SetParameter(std::string name) 
{ 
 if (name == "forcevector") 
 {
  AssignParameter("fx", GetNextWord());
  AssignParameter("fy", GetNextWord());
  AssignParameter("fz", GetNextWord());
  fx = GetDouble("fx");
  fy = GetDouble("fy");
  fz = GetDouble("fz");
 }
}

void ConstantForcePotential::Show() const
{
 Module::Show();
 std::cout << "   fx = " << fx << '\n';
 std::cout << "   fy = " << fy << '\n';
 std::cout << "   fz = " << fz << '\n';
}

std::string ConstantForcePotential::Keywords() const { return "forcevector"; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ConstantForcePotential(args); }
void destroy(Module * m) { delete m; }


